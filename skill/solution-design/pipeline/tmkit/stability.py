"""Stability-oriented, sequence + structure derived gates.

Complement the AF2 confidence gates with:
- residue solubility propensity (Goverde-style, clamped so over-polarization
  earns no extra credit),
- residual-hydrophobicity fraction over the redesigned positions,
- inter-chain interface contact detection (guards oligomer interfaces).
"""
import numpy as np

# aqueous solubility propensity, higher = more soluble. Polar/charged top,
# bulky hydrophobic bottom. Contributions are clamped (see POLAR_CAP).
SOLUBILITY = {
    'A': 0.50, 'R': 0.80, 'N': 0.85, 'D': 0.95, 'C': 0.30,
    'Q': 0.90, 'E': 0.95, 'G': 0.60, 'H': 0.70, 'I': 0.10,
    'L': 0.15, 'K': 0.85, 'M': 0.25, 'F': 0.10, 'P': 0.55,
    'S': 0.90, 'T': 0.80, 'W': 0.05, 'Y': 0.35, 'V': 0.15,
    'X': 0.40, '-': 0.50,
}
POLAR_CAP = 0.85          # beyond this, extra polarity is not rewarded
RESID_HYD_THR = 0.30      # propensity below this counts as residual hydrophobic


def sequence_solubility(seq, ref_seq, positions, aa1_map=None):
    """Score the redesigned subset of one chain.

    seq / ref_seq: full-length strings (reference-numbered).
    positions: reference positions (1-based) that were redesigned.
    Returns dict with soluble_score, residual_hydro_fraction, over_polar_fraction.
    """
    sol = []
    resid_hyd = 0
    over_polar = 0
    n = 0
    for p in positions:
        if p < 1 or p > len(seq):
            continue
        aa = seq[p - 1]
        if aa == '-' or aa == 'X':
            continue
        n += 1
        s = SOLUBILITY.get(aa, 0.40)
        sol.append(min(s, POLAR_CAP))
        if s <= RESID_HYD_THR:
            resid_hyd += 1
        if s >= POLAR_CAP:
            over_polar += 1
    return {
        'soluble_score': float(np.mean(sol)) if n else float('nan'),
        'residual_hydro': resid_hyd / n if n else float('nan'),
        'over_polar': over_polar / n if n else float('nan'),
    }


def designable_positions(mask, chain_id):
    return [refn for refn, info in mask[chain_id].items()
            if info['designable']]


def interface_contact_residues(target, chain_ids, dist=6.0):
    """Set of (chain, refn) residues whose any heavy atom contacts another chain.

    Leaves interface residues free only if allowed by the caller; used to
    hard-fix oligomer interfaces in the mask.
    """
    if len(chain_ids) < 2:
        return set()
    atoms = []
    for cid in chain_ids:
        for r in target.chains[cid].residues:
            for a in r.get_atoms():
                if a.element != 'H':
                    atoms.append((cid, r.id[1], a.get_coord()))

    contact = set()
    try:
        from scipy.spatial import cKDTree
        coords = np.asarray([a[2] for a in atoms])
        tree = cKDTree(coords)
        pairs = tree.query_pairs(dist, output_type='set')
        for i, j in pairs:
            a, b = atoms[i], atoms[j]
            if a[0] != b[0]:
                contact.add((a[0], a[1]))
                contact.add((b[0], b[1]))
        return contact
    except ImportError:
        pass

    # Biopython NeighborSearch (C KD-tree) fallback
    from Bio.PDB import NeighborSearch
    bio_atoms = []
    meta = {}
    for cid in chain_ids:
        for r in target.chains[cid].residues:
            for a in r.get_atoms():
                if a.element != 'H':
                    bio_atoms.append(a)
                    meta[id(a)] = (cid, r.id[1])
    ns = NeighborSearch(bio_atoms)
    for a in bio_atoms:
        ca, ra = meta[id(a)]
        for b in ns.search(a.get_coord(), dist):
            cb, rb = meta[id(b)]
            if ca != cb:
                contact.add((ca, ra))
                contact.add((cb, rb))
    return contact