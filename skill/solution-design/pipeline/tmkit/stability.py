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

# Charge / solubility heuristics. Acidic (-) and basic (+) pKa values, plus
# Eisenberg consensus hydrophobicity for helical hydrophobic moments.
PKA = {'D': 3.9, 'E': 4.2, 'C': 8.3, 'Y': 10.1, 'H': 6.0,
       'K': 10.5, 'R': 12.5}
ACIDIC = {'D', 'E', 'C', 'Y'}
BASIC = {'K', 'R', 'H'}
EISENBERG = {
    'A': 0.62, 'R': -2.53, 'N': -0.78, 'D': -0.90, 'C': 0.29,
    'Q': -0.85, 'E': -0.74, 'G': 0.48, 'H': -0.40, 'I': 1.38,
    'L': 1.06, 'K': -1.50, 'M': 0.64, 'F': 1.19, 'P': 0.12,
    'S': -0.18, 'T': -0.05, 'W': 0.81, 'Y': 0.26, 'V': 1.08,
}


def net_charge(seq, ph=7.0, termini=True):
    """Approximate net formal charge of a sequence at the given pH."""
    charge = 0.0
    for aa in seq:
        if aa in ACIDIC:
            charge -= 1.0 / (1.0 + 10.0 ** (PKA[aa] - ph))
        elif aa in BASIC:
            charge += 1.0 / (1.0 + 10.0 ** (ph - PKA[aa]))
    if termini:
        charge += 1.0 / (1.0 + 10.0 ** (ph - 9.6))       # N-terminus
        charge -= 1.0 / (1.0 + 10.0 ** (2.4 - ph))       # C-terminus
    return float(charge)


def isoelectric_point(seq):
    """Bisection estimate of the pH where net charge crosses zero."""
    lo, hi = 0.0, 14.0
    for _ in range(40):
        mid = (lo + hi) / 2.0
        if net_charge(seq, mid) > 0:
            lo = mid
        else:
            hi = mid
    return float((lo + hi) / 2.0)


def hydrophobic_moment(seq, tm_segments=None, window=11, angle=100.0):
    """Max per-residue Eisenberg hydrophobic moment over TM windows."""
    positions = [i for a, b in (tm_segments or [(1, len(seq))])
                 for i in range(int(a) - 1, int(b))]
    if not positions:
        positions = list(range(len(seq)))
    step = 1.0 if window <= 1 else 2.0 * np.pi * float(angle) / 360.0
    phase = [np.cos(step * i) for i in range(window)]
    quad = [np.sin(step * i) for i in range(window)]
    best = 0.0
    half = window // 2
    for center in positions:
        start = center - half
        vals = [EISENBERG.get(seq[i], 0.0)
                for i in range(start, start + window)
                if 0 <= i < len(seq)]
        if len(vals) < 3:
            continue
        c = sum(h * o for h, o in zip(vals, phase[:len(vals)]))
        s = sum(h * o for h, o in zip(vals, quad[:len(vals)]))
        best = max(best, float(np.hypot(c, s) / len(vals)))
    return best


def max_hydrophobic_window(seq, positions, window=9, hyd_cut=RESID_HYD_THR):
    """Largest fraction of hydrophobic residues in a designed-position run.

    A sequence-level proxy for exposed hydrophobic patches: consecutive
    redesigned positions are scanned in windows of `window`; a high value
    means hydrophobic residues cluster rather than being dispersed.
    """
    pos = sorted(p for p in positions if 1 <= p <= len(seq))
    if len(pos) < window:
        return float('nan')
    best = 0.0
    for i in range(0, len(pos) - window + 1):
        chunk = pos[i:i + window]
        hydro = sum(1 for p in chunk
                    if SOLUBILITY.get(seq[p - 1], 0.4) <= hyd_cut)
        best = max(best, hydro / float(window))
    return best


def sequence_properties(seq, positions, tm_segments=None):
    charge = net_charge(seq)
    moment = hydrophobic_moment(seq, tm_segments)
    patch = max_hydrophobic_window(seq, positions)
    return {
        'net_charge': charge,
        'abs_net_charge': abs(charge),
        'isoelectric_point': isoelectric_point(seq),
        'hydrophobic_moment': moment,
        'max_hydro_window': patch,
    }


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