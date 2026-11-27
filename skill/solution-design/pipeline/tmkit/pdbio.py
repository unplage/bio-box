import hashlib
import os
import tempfile
import numpy as np
from pathlib import Path
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO, Select
from Bio.PDB.Polypeptide import is_aa
from Bio.Align import PairwiseAligner


_CIF_CACHE = {}


def is_cif(path):
    return str(path).lower().endswith(('.cif', '.mmcif'))


def as_pdb(path, workdir=None):
    """Return a PDB path for a PDB or mmCIF input.

    Downstream code (fixed-column parsing, freesasa) expects PDB, so mmCIF
    models (Boltz-2 / AlphaFold3 / PDBx) are converted once and cached. B-factor
    (pLDDT) and chain/residue ids are preserved by the Biopython round trip.
    """
    path = Path(path)
    if not is_cif(path):
        return path
    key = (str(path.resolve()), path.stat().st_mtime)
    cached = _CIF_CACHE.get(key)
    if cached and Path(cached).exists():
        return Path(cached)
    base = Path(workdir or tempfile.gettempdir()) / 'tmkit_cif'
    base.mkdir(parents=True, exist_ok=True)
    out = base / (hashlib.md5(repr(key).encode()).hexdigest() + '.pdb')
    try:
        struct = MMCIFParser(QUIET=True, auth_chains=True,
                             auth_residues=True).get_structure(
                                 path.stem, str(path))
    except Exception:
        struct = MMCIFParser(QUIET=True, auth_chains=False,
                             auth_residues=False).get_structure(
                                 path.stem, str(path))
    io = PDBIO()
    io.set_structure(struct)
    io.save(str(out))
    _CIF_CACHE[key] = str(out)
    return out


AMINO3 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'MSE': 'M', 'SEC': 'C', 'PYL': 'O', 'HOH': 'O',
}


class ChainData:
    def __init__(self, chain_id, residues):
        self.chain_id = chain_id
        self.residues = residues
        self.raw_seq = ''.join(AMINO3.get(r.resname, 'X') for r in residues)


class Target:
    def __init__(self, structure, chains, reference_map=None, path=None):
        self.structure = structure
        self.chains = chains
        self.reference_map = reference_map or {}
        self.path = path


def parse_pdb(path, chain_ids=None):
    name = Path(path).stem
    pdb = as_pdb(path)
    model = PDBParser(QUIET=True).get_structure(name, str(pdb))[0]
    chains = {}
    for cid in chain_ids or [c.id for c in model]:
        try:
            ch = model[cid]
        except KeyError:
            print(f'[warn] chain {cid} not in {path}; skip')
            continue
        residues = [r for r in ch if r.id[0] == ' ' and is_aa(r)]
        chains[cid] = ChainData(cid, residues)
    return Target(model, chains, path=str(pdb))


def load_target(path, ecd_segments=None):
    target = parse_pdb(path)
    ecd_lookup = []
    for seg in ecd_segments or []:
        if isinstance(seg, (list, tuple)) and len(seg) == 2:
            ecd_lookup.append((int(seg[0]), int(seg[1])))
    return target, ecd_lookup


def _chains_are_homomeric(target, chain_ids):
    """True when chains are the same protein (residues agree at common
    numbers), even if their modelled ranges (termini) differ slightly."""
    if len(chain_ids) <= 1:
        return False
    maps = [{r.id[1]: r.resname for r in target.chains[cid].residues}
            for cid in chain_ids]
    common = set.intersection(*[set(m) for m in maps])
    if len(common) < 10:
        return False
    same = sum(1 for n in common if len({m[n] for m in maps}) == 1)
    return same / len(common) > 0.9


def harmonize_chains(target, chain_ids):
    """Restrict homomeric chains to their common residue numbers.

    Only homomers are intersected, so the same ordered positions can be
    compared/tied. Heteromers keep each chain's own residue set (their
    numbering is handled per chain by reference alignment).
    """
    homomer = _chains_are_homomeric(target, chain_ids)
    if homomer:
        sets = [set(r.id[1] for r in target.chains[cid].residues)
                for cid in chain_ids]
        common = set.intersection(*sets) if sets else set()
        for cid in chain_ids:
            chain = target.chains[cid]
            seqmap = {res.id[1]: res for res in chain.residues}
            kept = [seqmap[n] for n in sorted(common)]
            target.chains[cid] = ChainData(cid, kept)
    for cid in chain_ids:
        target.reference_map.pop(cid, None)
    return target


def aa1(resname):
    return AMINO3.get(resname.strip(), 'X')


def chain_sequence_from_pdb(pdb_path, chain_id):
    """[(pdb_resnum, aa1)] of CA atoms for one chain of a raw PDB file.

    Model files (ColabFold/AlphaFold) renumber residues 1..N per chain and
    never carry gaps; this function only needs the raw numbers for pairing.
    """
    pairs = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                if line[21] == chain_id:
                    pairs.append((int(line[22:26]), aa1(line[17:20])))
    return pairs


def map_model_to_reference(pdb_path, chain_id, ref_seq):
    """Map model PDB residue numbers -> reference positions (1-based) by
    global pairwise sequence alignment. Robust to renumbering and gaps.

    Returns {model_resnum: refn}. Order-based identity mapping when the
    model chain sequence equals the reference sequence.
    """
    pairs = chain_sequence_from_pdb(pdb_path, chain_id)
    if not pairs:
        return {}
    seq = ''.join(aa for _, aa in pairs)
    if len(seq) == len(ref_seq) and seq == ref_seq:
        return {n: i + 1 for i, (n, _) in enumerate(pairs)}
    aln = _aligner().align(seq, ref_seq)[0]
    mapping = {}
    for pi, ri in zip(*aln.indices):
        if pi < 0 or ri < 0:
            continue
        mapping[pairs[pi][0]] = ri + 1
    return mapping


def _aligner():
    """Semi-global pairwise aligner: free end gaps so partial chains (e.g.
    cryo-EM fragments) align to the correct reference sub-range instead of
    being forced to the reference termini."""
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5
    if hasattr(aligner, 'open_end_gap_score'):
        aligner.open_end_gap_score = 0
        aligner.extend_end_gap_score = 0
    else:
        aligner.end_open_gap_score = 0
        aligner.end_extend_gap_score = 0
    return aligner


def reference_numbering(target, chain_id, ref_seq):
    chain = target.chains[chain_id]
    if len(chain.raw_seq) != len(ref_seq):
        aln = _aligner().align(chain.raw_seq, ref_seq)[0]
        pdb_marks, ref_marks = aln.indices
        mapping = {}
        for pi, ri in zip(pdb_marks, ref_marks):
            if pi < 0 or ri < 0:
                continue
            mapping[chain.residues[pi].id[1]] = ri + 1
    else:
        # order-based mapping: i-th residue (sorted by number) -> reference pos i
        ordered = sorted(chain.residues, key=lambda r: r.id[1])
        mapping = {r.id[1]: i + 1 for i, r in enumerate(ordered)}
    target.reference_map[chain_id] = mapping
    return target.reference_map[chain_id]


def chain_ordinal(target, chain_id):
    return {r.id[1]: i + 1 for i, r in enumerate(target.chains[chain_id].residues)}


def mpnn_positions(target, chain_id):
    """Positions (1-based) of residues in the ProteinMPNN chain coordinate
    system: the PDB residue-number range (min..max), with X inserted at
    unmodeled positions. pos(pdb_res) = pdb_res - min_number + 1."""
    res = target.chains[chain_id].residues
    if not res:
        return {}
    first = min(r.id[1] for r in res)
    return {r.id[1]: r.id[1] - first + 1 for r in res}


def map_to_reference(target, chain_id):
    return target.reference_map.get(chain_id, {})


def ca_coords(target, chain_id, ref_numbers=None):
    mapping = target.reference_map.get(chain_id, {})
    coords, refs = [], []
    for r in target.chains[chain_id].residues:
        refn = mapping.get(r.id[1], r.id[1])
        if ref_numbers is not None and refn not in ref_numbers:
            continue
        try:
            ca = r['CA']
        except KeyError:
            continue
        coords.append(ca.get_coord())
        refs.append(refn)
    return (np.asarray(coords), refs) if coords else (np.zeros((0, 3)), [])


def full_length_sequence(ref_seq, designed_positions, designed_seq):
    design_map = dict(zip(designed_positions, designed_seq))
    return ''.join(design_map.get(i + 1, aa) for i, aa in enumerate(ref_seq))


class _ChainSelect(Select):
    def __init__(self, chains, omit_resnames=None, omit_resseq=None,
                 accept_resids=None):
        self.chains = chains
        self.omit_resnames = omit_resnames or set()
        self.omit_resseq = omit_resseq or set()
        self.accept_resids = accept_resids or {}

    def accept_chain(self, chain):
        return chain.id in self.chains

    def accept_residue(self, residue):
        parent = residue.get_parent()
        cid = parent.id if parent is not None else ''
        if self.accept_resids and residue.id[1] not in \
                self.accept_resids.get(cid, ()):
            return 0
        if residue.id[0] != ' ':
            return 0
        if residue.resname == 'HOH':
            return 0
        if residue.resname in self.omit_resnames:
            return 0
        if residue.id[1] in self.omit_resseq:
            return 0
        return 1 if AMINO3.get(residue.resname) else 0


def write_structure(target, path, chains=None, omit_resnames=None,
                    omit_resseq=None):
    """Write only the residues tracked in target.chains (harmonized), so the
    output chain order/ordinals exactly match chain_ordinal / ChainData."""
    io = PDBIO()
    chains = chains or list(target.chains)
    accept = {}
    for cid in chains:
        data = target.chains.get(cid)
        if data is not None:
            accept[cid] = {r.id[1] for r in data.residues}
    sel = _ChainSelect(chains, omit_resnames, omit_resseq, accept)
    io.set_structure(target.structure)
    io.save(str(path), sel)


def write_fasta(path, records):
    lines = []
    for name, seq in records:
        lines.append(f'>{name}')
        lines.append(seq)
    Path(path).write_text('\n'.join(lines) + '\n')


MAX_ASA = {
    'A': 121.0, 'R': 265.0, 'N': 187.0, 'D': 187.0, 'C': 148.0,
    'Q': 198.0, 'E': 223.0, 'G': 97.0, 'H': 216.0, 'I': 195.0,
    'L': 191.0, 'K': 236.0, 'M': 203.0, 'F': 228.0, 'P': 154.0,
    'S': 143.0, 'T': 163.0, 'W': 285.0, 'Y': 263.0, 'V': 165.0,
}


def residue_sasa_from_file(pdb_path):
    import freesasa
    structure = freesasa.Structure(str(as_pdb(pdb_path)))
    result = freesasa.calc(structure)
    sasa = {}
    for chain, residues in result.residueAreas().items():
        for num, area in residues.items():
            sasa[(chain, int(num))] = area.total
    return sasa


def residue_sasa_map(target, chains=None, workdir='/tmp/opencode'):
    Path(workdir).mkdir(parents=True, exist_ok=True)
    p = Path(workdir) / 'sasa_input.pdb'
    write_structure(target, p, chains=chains or list(target.chains))
    return residue_sasa_from_file(p)

def hetero_residues(structure, resnames=None, chain_ids=None):
    """Return ligand/cofactor residues from a Bio.PDB structure/model.

    Each entry is ``{'chain', 'resid', 'resname', 'coords'}`` with heavy-atom
    coordinates.  Waters are excluded.  ``resnames`` is case-insensitive.
    """
    wanted = {str(x).upper() for x in resnames} if resnames else None
    out = []
    for chain in structure:
        if chain_ids and chain.id not in chain_ids:
            continue
        for residue in chain:
            if residue.id[0] == ' ' or residue.resname == 'HOH':
                continue
            if wanted and residue.resname.upper() not in wanted:
                continue
            coords = [a.get_coord() for a in residue
                      if a.element != 'H']
            if not coords:
                continue
            out.append({
                'chain': chain.id, 'resid': residue.id[1],
                'resname': residue.resname, 'coords': np.asarray(coords, float),
            })
    return out


def relative_sasa_map(target, chains=None, workdir='/tmp/opencode'):
    """Per-residue SASA divided by the residue-type maximum (0-1)."""
    absolute = residue_sasa_map(target, chains=chains, workdir=workdir)
    out = {}
    for cid, data in (target.chains or {}).items():
        if chains and cid not in chains:
            continue
        by_num = {r.id[1]: r.resname for r in data.residues}
        for (chain, num), area in absolute.items():
            if chain != cid:
                continue
            aa = AMINO3.get(by_num.get(num, ''), 'X')
            out[(chain, num)] = area / MAX_ASA.get(aa, 1.0)
    return out
