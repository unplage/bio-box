"""Config normalization and validation for the generic membrane pipeline.

The pipeline accepts a deliberately small, design-oriented schema:

    membrane_class, tm_segments/ecd_segments/icd_segments, nterm/cterm,
    first_loop, design_scope, fixed_regions/fixed_residues, pockets,
    disulfides/ptms/epitopes, interface_mode, membrane, assembly,
    design_mode/backbone_dir, predict_tool.

This module fills safe defaults and reports warnings without silently changing
the biological interpretation.  Errors are raised by ``resolve_topology`` and
``build_mask``.
"""
from pathlib import Path

from .topology import MEMBRANE_CLASSES, normalize_membrane_class, \
    resolve_topology, resolve_topologies

_SCOPE_ALIASES = {
    'designable': 'tm_surface',
    'tmd_only': 'tm_surface',
    'default': 'tm_surface',
    'tm_surface': 'tm_surface',
    'surface': 'surface',
    'non_ecd': 'surface',
    'non-ecd': 'surface',
    'tm_core': 'tm_core',
    'tm_plus_icd': 'tm_plus_icd',
    'all_non_ecd': 'all_non_ecd',
    'all': 'all',
    'scaffold': 'all_non_ecd',
}


def normalize_config(cfg):
    """Return a copy of ``cfg`` with explicit, design-oriented defaults."""
    out = dict(cfg or {})
    chains = list(out.get('chains') or [])
    if not chains:
        raise ValueError('config must define chains')

    if not out.get('membrane_class') and not out.get('class'):
        out['class'] = 'alpha_helical_multi' if out.get('tm_segments') \
            else 'peripheral'
    out['membrane_class'] = normalize_membrane_class(out)

    out.setdefault('fix_ecd', True)
    out.setdefault('interface_mode',
                   'fixed' if out.get('fix_interface_contacts', True)
                   else 'none')
    out.setdefault('assembly', 'as_is')
    out.setdefault('design_mode', 'sequence')
    out.setdefault('strict_predict', True)
    out.setdefault('disulfides', 'auto')
    out.setdefault('ptms', False)
    out.setdefault('fix_epitopes', True)
    out.setdefault('fix_pore', True)
    out.setdefault('exclude_pore', True)
    out.setdefault('membrane_auto', False)
    out.setdefault('workdir', '/tmp/opencode')
    out.setdefault('interface_contact_dist', 6.0)
    out.setdefault('interface_core_contacts', 3)
    out.setdefault('pocket_shell_dist', 8.0)
    out.setdefault('fix_pocket_shell', True)
    out.setdefault('tm_junction_keep', 2)
    out.setdefault('tm_end_keep', 0)
    out.setdefault('surface_sasa_threshold', 40.0)
    out.setdefault('surface_rsa_threshold', 0.25)
    out.setdefault('relative_sasa', False)
    out.setdefault('cofactor_shell', 4.5)
    out.setdefault('symmetry', {})
    out.setdefault('interface_pairs', [])
    out.setdefault('cofactors', [])

    if not out.get('design_scope') and not out.get('include_regions') \
            and not out.get('include') \
            and out['membrane_class'] == 'peripheral':
        out['design_scope'] = 'surface'

    if not out.get('design_scope'):
        raw = out.get('include_regions', out.get('include', 'designable'))
        if isinstance(raw, (list, tuple)):
            out['design_scope'] = [str(x).upper() for x in raw]
        else:
            out['design_scope'] = _SCOPE_ALIASES.get(
                str(raw).strip().lower(), 'tm_surface')
    return out


def resolve_template_path(cfg, pdb_path):
    """Choose the structure to use.

    ``assembly: prebuilt`` + ``assembly_pdb`` (or ``assembly: /path.pdb``)
    lets the user supply a biological assembly without the pipeline guessing
    symmetry.  Automatic BIOMT expansion is intentionally not attempted: a
    wrong assembly is worse than an explicit warning.
    """
    assembly = cfg.get('assembly', 'as_is')
    if isinstance(assembly, str) and assembly.lower() not in (
            'as_is', 'asis', 'none', 'prebuilt', 'biological'):
        return assembly
    if cfg.get('assembly_pdb'):
        return cfg['assembly_pdb']
    return pdb_path


def validate_config(cfg, ref_lengths=None):
    """Return a list of warnings; raise for structural config errors."""
    cfg = normalize_config(cfg)
    chains = cfg['chains']
    topologies = resolve_topologies(cfg, ref_lengths, chains)

    warnings = []
    seen = set()

    def warn(message):
        if message not in seen:
            seen.add(message)
            warnings.append(message)

    for cid, topology in topologies.items():
        for message in topology.get('warnings') or []:
            warn(f'{cid}: {message}')

    refs = cfg.get('reference_sequences') or {}
    missing = [c for c in chains if c not in refs]
    if missing:
        raise ValueError(f'reference_sequences missing for chains {missing}')

    for cid, topology in topologies.items():
        cls = topology['class']
        if cls in MEMBRANE_CLASSES and cls != 'peripheral' \
                and not (cfg.get('membrane') or cfg.get('membrane_opm_pdb')
                         or cfg.get('opm_pdb') or cfg.get('membrane_auto')):
            warn(
                'no membrane slab configured; designability falls back to '
                'SASA. Provide `membrane`/`membrane_opm_pdb` or '
                '`membrane_auto: true` for a lipid-facing definition.')
        if cls == 'beta_barrel' and not cfg.get('pore_residues'):
            warn(
                'beta_barrel without pore_residues: interior/pore-facing '
                'strands may be designed. Provide `pore_residues` (and '
                'optionally `pore_shell`) or verify the mask manually.')

    scope = cfg.get('design_scope')
    if scope == 'all' and not cfg.get('fix_ecd', True):
        warn(
            'design_scope=all with fix_ecd=false: ECD residues may be '
            'redesigned; this violates the usual ECD-fixed principle.')

    assembly = str(cfg.get('assembly', 'as_is')).lower()
    if assembly == 'biological':
        warn(
            'assembly=biological is not expanded automatically; supply a '
            'prebuilt assembly PDB via `assembly_pdb` or `assembly: /path.pdb`.')
    if assembly == 'prebuilt' and not cfg.get('assembly_pdb') \
            and not cfg.get('pdb') and not cfg.get('template_pdb'):
        warn(
            'assembly=prebuilt without assembly_pdb: using the template PDB '
            'as supplied.')

    if str(cfg.get('design_mode', 'sequence')).lower() in (
            'scaffold', 'rf_diffusion', 'rfdiffusion',
            'partial_diffusion') and not cfg.get('backbone_dir') \
            and not cfg.get('rfdiffusion_script') \
            and not cfg.get('rfdiffusion_dir'):
        warn(
            'scaffold/RFdiffusion design mode without backbone_dir or '
            'rfdiffusion_script: the design step will fall back to the '
            'template backbone.')

    try:
        end_keep = int(cfg.get('tm_end_keep', 0))
    except (TypeError, ValueError):
        raise ValueError('tm_end_keep must be a non-negative integer')
    if end_keep < 0:
        raise ValueError('tm_end_keep must be a non-negative integer')

    try:
        sasa_threshold = float(cfg.get('surface_sasa_threshold', 40.0))
    except (TypeError, ValueError):
        raise ValueError('surface_sasa_threshold must be a number')
    if sasa_threshold < 0:
        raise ValueError('surface_sasa_threshold must be non-negative')

    symmetry = cfg.get('symmetry') or {}
    groups = symmetry.get('groups') if isinstance(symmetry, dict) else symmetry
    if groups:
        assigned = set()
        for group in groups:
            if isinstance(group, str):
                group = [group]
            if not isinstance(group, (list, tuple)) or len(group) < 2:
                raise ValueError(
                    'symmetry.groups entries must contain at least two chains')
            for cid in group:
                if cid not in chains:
                    raise ValueError(
                        f'symmetry group references unknown chain {cid!r}')
                if cid in assigned:
                    raise ValueError(
                        f'chain {cid!r} appears in more than one symmetry group')
                assigned.add(cid)

    for pair in cfg.get('interface_pairs') or []:
        if not isinstance(pair, (list, tuple)) or len(pair) != 2:
            raise ValueError(
                f'interface_pairs entries must be [chain_a, chain_b], got '
                f'{pair!r}')
        a, b = pair
        if a == b or a not in chains or b not in chains:
            raise ValueError(
                f'invalid interface pair {pair!r}; chains must differ and '
                f'exist in {chains}')

    for i, cofactor in enumerate(cfg.get('cofactors') or []):
        if isinstance(cofactor, str):
            continue
        if not isinstance(cofactor, dict):
            raise ValueError(
                f'cofactors[{i}] must be a mapping or a residue-name string')
        chain = cofactor.get('chain')
        if chain and chain not in chains:
            raise ValueError(
                f'cofactors[{i}] references unknown chain {chain!r}')
        if not (cofactor.get('resnames') or cofactor.get('resname')):
            warn(
                f'cofactors[{i}] has no resnames/resname; it will not be '
                f'fixed by the mask.')

    return warnings
