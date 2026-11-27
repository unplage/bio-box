"""Designability mask for soluble redesign of membrane proteins.

The mask is built from three layers:

1. topology (``tmkit.topology``) -- which residues are TMD/ECD/ICD/LOOP,
2. constraints (``tmkit.constraints``) -- ECD hard-fix, pockets, disulfides,
   PTMs, epitopes, pore, interface mode,
3. a design scope -- which regions may be designed and whether only the
   solvent-exposed TMD surface is eligible.

Nothing in this module is STEAP1-specific.  The same code handles a GPCR, a
single-pass receptor, a beta-barrel or a soluble control once the config
describes the topology and constraints.
"""
from .pdbio import map_to_reference, residue_sasa_map, relative_sasa_map, \
    mpnn_positions
from .topology import (MEMBRANE_CLASSES, auto_loop_topology, classify_region,
                       region_at, resolve_topology, resolve_topologies,
                       topology_sides)
from .constraints import build_constraints, symmetry_groups
from .membrane import membrane_geometries, membrane_exclusion, signed_distance

__all__ = [
    'MEMBRANE_CLASSES', 'auto_loop_topology', 'classify_region',
    'resolve_topology', 'topology_sides', 'build_mask', 'summary', 'resolve',
    'fixed_positions_dict', 'symmetry_ties_dict',
]

_SCOPE_ALIASES = {
    'designable': 'tm_surface',
    'tmd_only': 'tm_surface',
    'default': 'tm_surface',
    'tm_surface': 'tm_surface',
    'tm_core': 'tm_core',
    'tm_plus_icd': 'tm_plus_icd',
    'all_non_ecd': 'all_non_ecd',
    'surface': 'surface',
    'all': 'all',
    'scaffold': 'all_non_ecd',
}

_SCOPE_REGIONS = {
    'tm_surface': {'TMD'},
    'tm_core': {'TMD'},
    'tm_plus_icd': {'TMD', 'ICD'},
    'all_non_ecd': {'TMD', 'ICD', 'LOOP'},
    'surface': {'TMD', 'ECD', 'ICD', 'LOOP', 'SOLUBLE'},
    'all': {'TMD', 'ECD', 'ICD', 'LOOP', 'SOLUBLE'},
}


def resolve_design_scope(config):
    """Return ``(scope_name, regions, surface_only)`` for a config.

    ``include`` / ``include_regions`` remain supported for backwards
    compatibility; ``design_scope`` is the preferred general spelling.
    """
    raw = config.get('design_scope')
    if raw is None:
        raw = config.get('include_regions', config.get('include',
                                                       'designable'))
    if isinstance(raw, (list, tuple)):
        regions = {str(x).strip().upper() for x in raw}
        return 'explicit', regions, bool(config.get('design_tm_surface_only',
                                                     False))
    key = _SCOPE_ALIASES.get(str(raw).strip().lower())
    if key is None:
        raise ValueError(
            f'unknown design_scope/include={raw!r}; choose '
            f'{", ".join(sorted(_SCOPE_ALIASES))} or an explicit region list')
    regions = set(_SCOPE_REGIONS[key])
    if 'design_tm_surface_only' in config:
        surface_only = bool(config['design_tm_surface_only'])
    else:
        surface_only = key in ('tm_surface', 'surface')
    return key, regions, surface_only


def _auto_membrane(config, topology):
    """Materialize ``membrane_auto: true`` into an explicit membrane slab.

    A slab is the physically correct definition of a lipid-facing surface, but
    the PCA normal is only a heuristic; it is opt-in.  An explicit ``membrane``
    block or OPM/PPM structure takes precedence.
    """
    if config.get('membrane') is False:
        return config
    if config.get('membrane') is not None or config.get('membrane_opm_pdb') \
            or config.get('opm_pdb'):
        return config
    if topology['class'] not in MEMBRANE_CLASSES or \
            topology['class'] == 'peripheral':
        return config
    # A PCA-derived membrane normal is only a heuristic.  Do not silently
    # apply it: require `membrane_auto: true` or an explicit membrane/OPM
    # definition, and let the annotate step warn when it is absent.
    if config.get('membrane_auto', False):
        out = dict(config)
        out['membrane'] = {
            'normal': 'auto',
            'thickness': float(config.get('membrane_thickness', 30.0)),
            'headgroup': float(config.get('membrane_headgroup', 15.0)),
            'exclude_aromatic_belt': bool(
                config.get('exclude_aromatic_belt', True)),
        }
        return out
    return config


def _near_junction(refn, topology, keep, tm_end_keep=0):
    for segs in (topology.get('ecd_segments') or [],
                 topology.get('icd_segments') or []):
        for a, b in segs:
            if abs(refn - a) <= keep or abs(refn - b) <= keep:
                return True
    for a, b in topology.get('tm_segments') or []:
        if a <= refn <= b and tm_end_keep > 0:
            if refn - a < tm_end_keep or b - refn < tm_end_keep:
                return True
    return False


def build_mask(target, chain_ids, config):
    config = dict(config)
    ref_seqs = config.get('reference_sequences', {})
    ref_lengths = {c: len(ref_seqs[c]) for c in chain_ids if c in ref_seqs}
    topologies = resolve_topologies(config, ref_lengths, chain_ids)
    topology = topologies.get(chain_ids[0]) if chain_ids else \
        resolve_topology(config, ref_lengths)
    config = _auto_membrane(config, topology)

    constraints = build_constraints(target, chain_ids, config, topology)
    fixed_regions = constraints['fixed_regions']
    _, include_regions, surface_only = resolve_design_scope(config)
    relative_sasa = bool(config.get('relative_sasa', False))
    if relative_sasa:
        surface_thr = float(config.get('surface_rsa_threshold', 0.25))
    else:
        surface_thr = float(config.get('surface_sasa_threshold', 40.0))
    junction_keep = int(config.get('tm_junction_keep', 2))
    tm_end_keep = int(config.get('tm_end_keep', 0))
    interface_mode = constraints['interface_mode']

    sasa = {}
    if surface_only or config.get('use_sasa', True):
        sasa_map = relative_sasa_map if relative_sasa else residue_sasa_map
        sasa = sasa_map(target, chains=chain_ids,
                        workdir=config.get('workdir', '/tmp/opencode'))

    geoms = membrane_geometries(target, chain_ids, topologies, config)

    mask = {}
    for cid in chain_ids:
        mapping = map_to_reference(target, cid)
        chain_topology = topologies.get(cid, topology)
        mask[cid] = {}
        for r in target.chains[cid].residues:
            pdb_res = r.id[1]
            refn = mapping.get(pdb_res, pdb_res)
            region = region_at(refn, chain_topology)
            designable = region in include_regions
            reason = 'fixed'
            a = float('inf')
            mem_d = None

            if region == 'ECD' and constraints['fix_ecd']:
                designable = False
                reason = 'ecd'
            elif region == 'SIGNAL':
                designable = False
                reason = 'signal'
            elif region in fixed_regions:
                designable = False
                reason = region.lower()
            elif not designable:
                reason = 'fixed_region'
            else:
                fixed_reason = constraints['fixed_reasons'].get(cid, {}).get(
                    refn)
                if fixed_reason:
                    designable = False
                    reason = fixed_reason
                elif surface_only:
                    a = sasa.get((cid, pdb_res), float('inf'))
                    if a <= surface_thr:
                        designable = False
                        reason = 'core'
                if designable and _near_junction(refn, chain_topology,
                                                 junction_keep, tm_end_keep):
                    designable = False
                    reason = 'junction'
                if designable and (cid, pdb_res) in constraints[
                        'interface_core']:
                    designable = False
                    reason = ('interface' if interface_mode == 'fixed'
                              else 'interface_core')
                if designable and refn in constraints['pore_shell'].get(
                        cid, set()) and config.get('exclude_pore', True):
                    designable = False
                    reason = 'pore'
                if designable and surface_only:
                    a = sasa.get((cid, pdb_res), float('inf'))
                    reason = 'surface'
                elif designable:
                    reason = 'designable'

            geom = geoms.get(cid)
            if geom is not None:
                try:
                    coord = r['CA'].get_coord()
                except KeyError:
                    coord = None
                if coord is not None:
                    mem_d = signed_distance(coord, geom)
                    exclusion = membrane_exclusion(coord, geom, r.resname)
                    if exclusion is not None and designable:
                        designable = False
                        reason = exclusion

            mask[cid][refn] = {
                'region': region, 'designable': designable,
                'reason': reason, 'pdb_residue': pdb_res,
                'sasa': a, 'interface': (cid, pdb_res) in constraints[
                    'interface_counts'],
                'interface_core': (cid, pdb_res) in constraints[
                    'interface_core'],
                'membrane_d': mem_d,
            }

    groups = symmetry_groups(config, chain_ids)
    if config.get('symmetric_mask', bool(groups)) and groups:
        for group in groups:
            common = None
            for cid in group:
                s = {refn for refn, info in mask[cid].items()
                     if info['designable']}
                common = s if common is None else (common & s)
            for cid in group:
                for refn, info in mask[cid].items():
                    if info['designable'] and refn not in common:
                        info['designable'] = False
                        info['reason'] = 'asym_masked'
    return mask


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
    return {'regions': regions, 'designable': counts['designable'],
            'reasons': reasons}


def resolve(target, chain_ids, mask):
    design = {}
    fixed = {}
    for cid in chain_ids:
        positions = mpnn_positions(target, cid)
        design[cid] = []
        fixed[cid] = []
        for refn in sorted(mask[cid]):
            info = mask[cid][refn]
            pdb_res = info['pdb_residue']
            pos = positions.get(pdb_res)
            if pos is None:
                continue
            if info['designable']:
                design[cid].append(pos)
            else:
                fixed[cid].append(pos)
    return design, fixed


def fixed_positions_dict(target, chain_ids, mask):
    _, fixed = resolve(target, chain_ids, mask)
    return fixed


def symmetry_ties_dict(target, chain_ids, mask, config=None):
    """Return one ProteinMPNN tie group per equivalent residue position.

    Each returned dict ties a single reference position across a symmetry
    group, e.g. ``{"A": [[5], [1.0]], "B": [[5], [1.0]]}``.  The old code
    returned one group containing every designable position of every chain,
    which forced all designed sites to the same amino acid.  Ties are only
    emitted inside explicit ``symmetry.groups`` (or an exact homomer); a
    heteromer gets no implicit ties.
    """
    if not chain_ids:
        return []
    groups = symmetry_groups(config or {}, chain_ids)
    if not groups:
        return []
    items = []
    for group in groups:
        group = [c for c in group if c in chain_ids]
        if len(group) < 2:
            continue
        by_ref = {}
        for cid in group:
            positions = mpnn_positions(target, cid)
            for refn, info in mask[cid].items():
                if not info.get('designable'):
                    continue
                pos = positions.get(info['pdb_residue'])
                if pos is not None:
                    by_ref.setdefault(refn, {})[cid] = pos
        for refn in sorted(by_ref):
            d = by_ref[refn]
            if set(d) == set(group):
                items.append({c: [[d[c]], [1.0]] for c in group})
    return items
