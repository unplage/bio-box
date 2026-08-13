import os
import numpy as np
from pathlib import Path
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa
from Bio.Align import PairwiseAligner


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
    model = PDBParser(QUIET=True).get_structure(name, str(path))[0]
    chains = {}
    for cid in chain_ids or [c.id for c in model]:
        try:
            ch = model[cid]
        except KeyError:
            print(f'[warn] chain {cid} not in {path}; skip')
            continue
        residues = [r for r in ch if r.id[0] == ' ' and is_aa(r)]
        chains[cid] = ChainData(cid, residues)
    return Target(model, chains, path=str(path))


def load_target(path, ecd_segments=None):
    target = parse_pdb(path)
    ecd_lookup = []
    for seg in ecd_segments or []:
        if isinstance(seg, (list, tuple)) and len(seg) == 2:
            ecd_lookup.append((int(seg[0]), int(seg[1])))
    return target, ecd_lookup


def harmonize_chains(target, chain_ids):
    sets = [set(r.id[1] for r in target.chains[cid].residues) for cid in chain_ids]
    common = set.intersection(*sets) if sets else set()
    for cid in chain_ids:
        chain = target.chains[cid]
        seqmap = {res.id[1]: res for res in chain.residues}
        kept = [seqmap[n] for n in sorted(common)]
        target.chains[cid] = ChainData(cid, kept)
    for cid in chain_ids:
        target.reference_map.pop(cid, None)
    return target


def reference_numbering(target, chain_id, ref_seq):
    chain = target.chains[chain_id]
    if len(chain.raw_seq) != len(ref_seq):
        aligner = PairwiseAligner()
        aligner.mode = 'global'
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -2
        aligner.extend_gap_score = -0.5
        aln = aligner.align(chain.raw_seq, ref_seq)[0]
        pdb_marks, ref_marks = aln.indices
        mapping = {}
        ref_pos = 0
        for pi, ri in zip(pdb_marks, ref_marks):
            if pi < 0 or ri < 0:
                continue
            mapping[chain.residues[pi].id[1]] = ref_pos + 1
            ref_pos += 1
    else:
        # order-based mapping: i-th residue (sorted by number) -> reference pos i
        ordered = sorted(chain.residues, key=lambda r: r.id[1])
        mapping = {r.id[1]: i + 1 for i, r in enumerate(ordered)}
    target.reference_map[chain_id] = mapping
    return target.reference_map[chain_id]


def chain_ordinal(target, chain_id):
    return {r.id[1]: i + 1 for i, r in enumerate(target.chains[chain_id].residues)}


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
    def __init__(self, chains, omit_resnames=None, omit_resseq=None):
        self.chains = chains
        self.omit_resnames = omit_resnames or set()
        self.omit_resseq = omit_resseq or set()

    def accept_chain(self, chain):
        return chain.id in self.chains

    def accept_residue(self, residue):
        if residue.id[0] != ' ':
            return 0
        if residue.resname == 'HOH':
            return 0
        if residue.resname in self.omit_resnames:
            return 0
        if residue.id[1] in self.omit_resseq:
            return 0
        return 1 if AMINO3.get(residue.resname) else 0


def write_structure(target, path, chains=None, omit_resnames=None, omit_resseq=None):
    io = PDBIO()
    sel = _ChainSelect(chains or list(target.chains), omit_resnames, omit_resseq)
    io.set_structure(target.structure)
    io.save(str(path), sel)


def write_fasta(path, records):
    lines = []
    for name, seq in records:
        lines.append(f'>{name}')
        lines.append(seq)
    Path(path).write_text('\n'.join(lines) + '\n')


def residue_sasa_from_file(pdb_path):
    import freesasa
    structure = freesasa.Structure(str(pdb_path))
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