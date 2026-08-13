import json
from .pdbio import chain_ordinal, map_to_reference, residue_sasa_map


def classify_region(refn, tm_segments, ecd_segments, icd_segments, nterm='ICD', cterm='ICD'):
    for a, b in tm_segments:
        if a <= refn <= b:
            return 'TMD'
    for a, b in ecd_segments:
        if a <= refn <= b:
            return 'ECD'
    for a, b in icd_segments:
        if a <= refn <= b:
            return 'ICD'
    if tm_segments and refn < tm_segments[0][0]:
        return nterm
    if tm_segments and refn > tm_segments[-1][1]:
        return cterm
    return 'LOOP'


def auto_loop_topology(tm_segments):
    ecd, icd = [], []
    for k in range(len(tm_segments) - 1):
        a = tm_segments[k][1] + 1
        b = tm_segments[k + 1][0] - 1
        if a <= b:
            if k % 2 == 0:
                ecd.append((a, b))
            else:
                icd.append((a, b))
    return ecd, icd


def build_mask(target, chain_ids, config):
    tm = [(int(a), int(b)) for a, b in config['tm_segments']]
    ecd = [(int(a), int(b)) for a, b in config.get('ecd_segments', [])]
    icd = [(int(a), int(b)) for a, b in config.get('icd_segments', [])]
    if not ecd and not icd:
        ecd, icd = auto_loop_topology(tm)
    include_regions = config.get('include_regions',
                                 config.get('include', 'designable'))
    if isinstance(include_regions, str):
        if include_regions in ('designable', 'tmd_only', 'default'):
            include_regions = ['TMD']
        elif include_regions == 'all':
            include_regions = ['TMD', 'ECD', 'ICD', 'LOOP']
        elif include_regions == 'tm_plus_icd':
            include_regions = ['TMD', 'ICD']
        else:
            include_regions = ['TMD']
    fixed_res = set(int(x) for x in config.get('fixed_residues', []))
    protect_pad = int(config.get('protect_pad', 3))
    protected = set()
    for f in fixed_res:
        for d in range(-protect_pad, protect_pad + 1):
            protected.add(f + d)
    surface_thr = float(config.get('surface_sasa_threshold', 40.0))
    surface_only = bool(config.get('design_tm_surface_only', True))
    junction_keep = int(config.get('tm_junction_keep', 2))

    sasa = {}
    if not surface_only or config.get('use_sasa', True):
        sasa = residue_sasa_map(target, chains=chain_ids,
                                workdir=config.get('workdir', '/tmp/opencode'))

    interface = set()
    fix_interface = bool(config.get('fix_interface_contacts', True))
    if fix_interface and len(chain_ids) > 1:
        from .stability import interface_contact_residues
        interface = interface_contact_residues(
            target, chain_ids,
            dist=float(config.get('interface_contact_dist', 6.0)))

    mask = {}
    for cid in chain_ids:
        mapping = map_to_reference(target, cid)
        mask[cid] = {}
        for r in target.chains[cid].residues:
            pdb_res = r.id[1]
            refn = mapping.get(pdb_res, pdb_res)
            region = classify_region(refn, tm, ecd, icd, config.get('nterm', 'ICD'),
                                     config.get('cterm', 'ICD'))
            is_interface = (cid, pdb_res) in interface
            designable = False
            reason = 'fixed'
            a = float('inf')
            if region in include_regions:
                key = (cid, pdb_res)
                a = sasa.get(key, float('inf'))
                junction = ecd if region == 'TMD' else []
                if surface_only and a <= surface_thr:
                    reason = 'core'
                elif _near_junction(refn, junction, junction_keep):
                    reason = 'junction'
                elif refn in protected:
                    reason = 'protected'
                elif is_interface:
                    reason = 'interface'
                else:
                    designable = True
                    reason = 'surface'
            mask[cid][refn] = {'region': region, 'designable': designable,
                               'reason': reason, 'pdb_residue': pdb_res,
                               'sasa': a, 'interface': is_interface}
    if config.get('symmetric_mask', True) and len(chain_ids) > 1:
        common = None
        for cid in chain_ids:
            s = {refn for refn, info in mask[cid].items() if info['designable']}
            common = s if common is None else (common & s)
        for cid in chain_ids:
            for refn, info in mask[cid].items():
                if info['designable'] and refn not in common:
                    info['designable'] = False
                    info['reason'] = 'asym_masked'
    return mask


def _near_junction(refn, ecd_segments, keep):
    for a, b in ecd_segments:
        if abs(refn - a) <= keep or abs(refn - b) <= keep:
            return True
    return False


def summary(mask):
    counts = {'designable': 0}
    reasons = {}
    regions = {}
    for cid, resmap in mask.items():
        for refn, info in resmap.items():
            regions[info['region']] = regions.get(info['region'], 0) + 1
            if info['designable']:
                counts['designable'] += 1
            reasons[info['reason']] = reasons.get(info['reason'], 0) + 1
    return {'regions': regions, 'designable': counts['designable'], 'reasons': reasons}


def resolve(target, chain_ids, mask):
    ord_map = []
    design = {}
    fixed = {}
    for cid in chain_ids:
        ordinal = chain_ordinal(target, cid)
        mapping = map_to_reference(target, cid)
        design[cid] = []
        fixed[cid] = []
        for refn in sorted(mask[cid]):
            info = mask[cid][refn]
            pdb_res = info['pdb_residue']
            if info['designable']:
                design[cid].append(ordinal[pdb_res])
            else:
                fixed[cid].append(ordinal[pdb_res])
    return design, fixed


def fixed_positions_dict(target, chain_ids, mask):
    _, fixed = resolve(target, chain_ids, mask)
    return fixed


def symmetry_ties_dict(target, chain_ids, mask):
    design, _ = resolve(target, chain_ids, mask)
    if not chain_ids:
        return {}
    first = chain_ids[0]
    if all(design[c] == design[chain_ids[0]] for c in chain_ids):
        positions = design[chain_ids[0]]
        item = {}
        for c in chain_ids:
            item[c] = [positions, [1.0] * len(positions)]
        return {'0': [item]}
    by_ref = {}
    for cid in chain_ids:
        ordinal = chain_ordinal(target, cid)
        for refn, info in mask[cid].items():
            if info['designable']:
                by_ref.setdefault(refn, {})[cid] = ordinal[info['pdb_residue']]
    items = []
    for refn in sorted(by_ref):
        d = by_ref[refn]
        if set(d) == set(chain_ids):
            item = {}
            for c in chain_ids:
                item[c] = [[d[c]], [1.0]]
            items.append(item)
    return {'0': items}