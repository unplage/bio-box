import numpy as np


APOLAR = set('GPAVILMFYW')


def kabsch(A, B):
    A = np.asarray(A, float)
    B = np.asarray(B, float)
    ca = A.mean(0)
    cb = B.mean(0)
    A_c = A - ca
    B_c = B - cb
    H = A_c.T @ B_c
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    t = cb - R @ ca
    return R, t


def rmsd(A, B):
    A = np.asarray(A, float)
    B = np.asarray(B, float)
    n = min(len(A), len(B))
    if n == 0:
        return float('nan')
    A, B = A[:n], B[:n]
    R, t = kabsch(A, B)
    d2 = ((B - (A @ R.T + t)) ** 2).sum(1)
    return float(np.sqrt(d2.mean()))


def tm_score(A, B):
    A = np.asarray(A, float)
    B = np.asarray(B, float)
    n = min(len(A), len(B))
    if n == 0:
        return float('nan')
    A, B = A[:n], B[:n]
    d0 = 1.24 * (n - 15) ** (1.0 / 3.0) - 1.8
    R, t = kabsch(A, B)
    d2 = ((B - (A @ R.T + t)) ** 2).sum(1)
    score = (1.0 / (1.0 + (d2 / d0 ** 2))).sum() / n
    return float(score)


def surface_hydrophobic_fraction(sasa, residue_aa_map):
    surface = [aa for (res, aa) in residue_aa_map.items()
               if sasa.get(res, 0.0) > 40.0]
    if not surface:
        return float('nan')
    hydrophobic = [aa for aa in surface if aa in APOLAR]
    return len(hydrophobic) / len(surface)


def seq_identity(ref_seq, mut_seq):
    if not ref_seq:
        return float('nan')
    same = sum(1 for a, b in zip(ref_seq, mut_seq) if a == b and a != 'X')
    total = sum(1 for a, b in zip(ref_seq, mut_seq) if a != 'X' and b != 'X')
    return same / total if total else float('nan')


def mutation_calls(ref_seq, mut_seq):
    calls = []
    for i, (a, b) in enumerate(zip(ref_seq, mut_seq)):
        if a != b and a != '-' and b != '-' and a != 'X' and b != 'X':
            calls.append((i + 1, a, b))
    return calls


def parse_plddt(pdb_path):
    from collections import defaultdict
    values = defaultdict(list)
    with open(pdb_path) as f:
        for line in f:
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                ch = line[21]
                res = int(line[22:26])
                try:
                    v = float(line[60:66])
                except ValueError:
                    continue
                values[(ch, res)].append(v)
    return {res: float(np.mean(v)) for (ch, res), v in values.items()}


def is_transferable(ref_base, transferable_aa, inter_loop_residue,
                    contact_residues, to_exclude):
    """Classify each reference position as transferable.

    ref_base: dict refnum -> info dict with 'sasa' and 'avg_hyd'.
    transferable_aa: optional iterable of explicitly transferable refnums.
    contact_residues / to_exclude: iterables of refnums.
    """
    transferable_aa = set(transferable_aa or [])
    contacs = set(contact_residues or [])
    exclude = set(to_exclude or [])
    status = {}
    for refn, info in (ref_base or {}).items():
        sasa = info.get('sasa', 0.0)
        is_surface = (sasa > 30.0)
        is_contact = (refn in contacs) or (refn in inter_loop_residue or [])
        trans = (refn in transferable_aa) or (is_surface and not is_contact)
        if refn in exclude:
            trans = False
        status[refn] = {
            'transferable': trans,
            'surface': is_surface,
            'contact': is_contact,
            'sasa': float(sasa),
            'avg_hyd': float(info.get('avg_hyd', 0.0)),
        }
    return status


def tm_rmsd(native_pdb, model_pdb, species_list, ecd_range, epitope_range,
            n_d_dig, n_r_dig, n_d_dig2=None, n_r_dig2=None, plddt=None,
            hold=None, transferred=None, excluded=None, mode='verify'):
    """Aggregate structural scoring between a native and a designed model.

    species_list: ordered chain ids to compare.
    ecd_range / epitope_range: (start, end) reference numbers.
    If both pdbs use the same per-chain residue numbering (reference
    numbering), ecd_range/others are comparable directly.
    Returns dict with tm, rmsd_ecd, rmsd_others, plddt.
    """
    from .pdbio import parse_pdb, ca_coords
    native_target = parse_pdb(native_pdb, chain_ids=species_list)
    model_target = parse_pdb(model_pdb, chain_ids=species_list)

    tm_vals, n_id = [], 0
    for cid in species_list:
        A, _ = ca_coords(native_target, cid)
        B, _ = ca_coords(model_target, cid)
        if len(A) and len(B):
            tm_vals.append(tm_score(A, B))
            n_id += 1
    if not tm_vals:
        return {'tm': float('nan')}

    common_ref = set()
    if ecd_range:
        common_ref |= set(range(int(ecd_range[0]), int(ecd_range[1]) + 1))
    if epitope_range:
        common_ref |= set(range(int(epitope_range[0]),
                                int(epitope_range[1]) + 1))
    ecd_n = {r: ca_coords(native_target, species_list[0], [r])[0][0]
             for r in common_ref
             if len(ca_coords(native_target, species_list[0], [r])[0])}
    ecd_m = {r: ca_coords(model_target, species_list[0], [r])[0][0]
             for r in common_ref
             if len(ca_coords(model_target, species_list[0], [r])[0])}
    common = sorted(set(ecd_n) & set(ecd_m))
    ecd_rmsd = rmsd(np.asarray([ecd_n[r] for r in common]),
                    np.asarray([ecd_m[r] for r in common])) if common \
        else float('nan')

    others_penalized = []
    for cid in species_list:
        A, refs = ca_coords(native_target, cid)
        B, _ = ca_coords(model_target, cid)
        if len(A) != len(B):
            continue
        keep = [i for i, r in enumerate(refs) if r not in common_ref]
        if keep:
            others_penalized.append(
                rmsd(A[keep], B[keep]) if len(keep) else float('nan'))
    rmsd_others = float(np.mean([v for v in others_penalized
                                 if not isinstance(v, float) or not np.isnan(v)])) \
        if others_penalized else float('nan')

    plddt_mean = float('nan')
    if plddt:
        try:
            plddt_mean = float(np.mean(list(parse_plddt(model_pdb).values())))
        except Exception:
            pass

    return {
        'tm': float(np.mean(tm_vals)),
        'rmsd_ecd': ecd_rmsd,
        'rmsd_others': rmsd_others,
        'plddt': plddt_mean,
    }