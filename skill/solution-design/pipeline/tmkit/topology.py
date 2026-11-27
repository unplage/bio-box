"""Membrane-protein topology model.

The pipeline is not STEAP1-specific: a target is described by a small number
of design-relevant axes instead of a protein-family label.

    membrane_class:
      alpha_helical_multi | single_pass | beta_barrel | tail_anchored | peripheral

For alpha-helical multi-pass and beta-barrel proteins, TM segments are ordered
along the sequence and inter-TM loops alternate between the two sides of the
membrane.  ``first_loop`` says which side the first inter-TM loop is on; the
N-terminus is on the opposite side.  Single-pass / tail-anchored targets have
no inter-TM loop, so their N/C sides must be stated explicitly (a warning is
emitted when the default type-I orientation is assumed).
"""
from collections import OrderedDict

MEMBRANE_CLASSES = (
    'alpha_helical_multi',
    'single_pass',
    'beta_barrel',
    'tail_anchored',
    'peripheral',
)

_CLASS_ALIASES = {
    'alpha_helical_multi': 'alpha_helical_multi',
    'alpha_helical': 'alpha_helical_multi',
    'alpha-helical': 'alpha_helical_multi',
    'multi_pass': 'alpha_helical_multi',
    'multipass': 'alpha_helical_multi',
    'tm_alpha': 'alpha_helical_multi',
    'single_pass': 'single_pass',
    'single-pass': 'single_pass',
    'single': 'single_pass',
    'bitopic': 'single_pass',
    'beta_barrel': 'beta_barrel',
    'beta-barrel': 'beta_barrel',
    'beta': 'beta_barrel',
    'omp': 'beta_barrel',
    'tail_anchored': 'tail_anchored',
    'tail-anchored': 'tail_anchored',
    'tail': 'tail_anchored',
    'peripheral': 'peripheral',
    'soluble': 'peripheral',
    'peripheral_membrane': 'peripheral',
}

_SIDES = ('ECD', 'ICD')


def normalize_membrane_class(config):
    raw = config.get('membrane_class', config.get('class', 'alpha_helical_multi'))
    key = str(raw).strip().lower()
    if key not in _CLASS_ALIASES:
        raise ValueError(
            f'unknown membrane_class={raw!r}; choose one of '
            f'{", ".join(MEMBRANE_CLASSES)}')
    return _CLASS_ALIASES[key]


def _segments(value, name):
    if value is None:
        return []
    if isinstance(value, dict):
        # Heteromer-specific segments are not meaningful for a global topology;
        # accept a single (start, end) pair or a list under a chain key.
        raise ValueError(f'{name} must be a list of [start, end] pairs')
    out = []
    for seg in value:
        if not isinstance(seg, (list, tuple)) or len(seg) != 2:
            raise ValueError(f'{name} entries must be [start, end], got {seg!r}')
        a, b = int(seg[0]), int(seg[1])
        if a > b:
            a, b = b, a
        out.append((a, b))
    return out


def _normalize_side(value, default=None, name='side'):
    if value is None:
        return default
    side = str(value).strip().upper()
    if side not in _SIDES:
        raise ValueError(f'{name} must be ECD or ICD, got {value!r}')
    return side


def auto_loop_topology(tm_segments, first_loop='ECD'):
    """Assign inter-TM loops alternately, starting with ``first_loop``.

    This is the topology rule for both alpha-helical multi-pass proteins and
    beta-barrels: consecutive loops lie on opposite sides of the membrane.
    """
    ecd, icd = [], []
    first = 'ECD' if str(first_loop).upper().startswith('E') else 'ICD'
    for k in range(len(tm_segments) - 1):
        a = tm_segments[k][1] + 1
        b = tm_segments[k + 1][0] - 1
        if a <= b:
            kind = first if k % 2 == 0 else \
                ('ICD' if first == 'ECD' else 'ECD')
            (ecd if kind == 'ECD' else icd).append((a, b))
    return ecd, icd


def topology_sides(tm_segments, config):
    """Infer (nterm, cterm, first_loop) sides from the first inter-TM loop.

    For a single TM there is no inter-TM loop; explicit ``nterm`` / ``cterm``
    in the config are authoritative.  When absent, the common type-I
    orientation (N-out / C-in) is assumed and a warning is returned through
    ``resolve_topology``.
    """
    first_loop = str(config.get('first_loop', 'ECD')).upper()
    first_loop = 'ECD' if first_loop.startswith('E') else 'ICD'
    if len(tm_segments) <= 1:
        nterm = _normalize_side(config.get('nterm'), 'ECD', 'nterm')
        cterm = _normalize_side(config.get('cterm'), 'ICD', 'cterm')
        return nterm, cterm, first_loop
    other = 'ICD' if first_loop == 'ECD' else 'ECD'
    nterm = _normalize_side(config.get('nterm'), other, 'nterm')
    cterm = _normalize_side(
        config.get('cterm'), other if len(tm_segments) % 2 == 0 else first_loop,
        'cterm')
    return nterm, cterm, first_loop


def _validate_segments(tm, ecd, icd, ref_lengths=None):
    for name, segs in (('tm_segments', tm), ('ecd_segments', ecd),
                       ('icd_segments', icd)):
        for i, (a, b) in enumerate(segs):
            if a < 1 or b < 1:
                raise ValueError(f'{name}[{i}] must use 1-based positions')
            if ref_lengths:
                for c, n in ref_lengths.items():
                    if b > n:
                        raise ValueError(
                            f'{name}[{i}]={a}-{b} exceeds chain {c} length {n}')
    def overlaps(segs):
        for i in range(len(segs)):
            for j in range(i + 1, len(segs)):
                if not (segs[i][1] < segs[j][0] or segs[j][1] < segs[i][0]):
                    return segs[i], segs[j]
        return None
    for name, segs in (('tm_segments', tm), ('ecd_segments', ecd),
                       ('icd_segments', icd)):
        hit = overlaps(segs)
        if hit:
            raise ValueError(f'{name} contains overlapping segments {hit}')
    for name, segs in (('ecd_segments', ecd), ('icd_segments', icd)):
        for e in segs:
            for t in tm:
                if not (e[1] < t[0] or t[1] < e[0]):
                    raise ValueError(
                        f'{name} segment {e} overlaps TM segment {t}')


def resolve_topology(config, ref_lengths=None):
    """Return a normalized topology dict.

    The dict always contains ``class``, ``tm_segments``, ``ecd_segments``,
    ``icd_segments``, ``nterm``, ``cterm``, ``first_loop``, ``tm_type`` and
    ``signal_peptide``.  ``warnings`` is a list of human-readable issues that
    do not invalidate the config (e.g. an inferred single-pass orientation).
    """
    membrane_class = normalize_membrane_class(config)
    tm = _segments(config.get('tm_segments'), 'tm_segments')
    tm.sort()
    ecd = _segments(config.get('ecd_segments'), 'ecd_segments')
    icd = _segments(config.get('icd_segments'), 'icd_segments')
    first_loop = str(config.get('first_loop', 'ECD')).upper()
    first_loop = 'ECD' if first_loop.startswith('E') else 'ICD'
    warnings = []

    if not ecd and not icd and tm:
        ecd, icd = auto_loop_topology(tm, first_loop)

    nterm, cterm, first_loop = topology_sides(tm, config)

    if membrane_class in ('single_pass', 'tail_anchored'):
        if len(tm) != 1:
            raise ValueError(
                f'membrane_class={membrane_class} requires exactly one '
                f'TM segment, got {len(tm)}')
        if not (config.get('nterm') and config.get('cterm')):
            warnings.append(
                'single-TM target without explicit nterm/cterm; assuming '
                'type-I N-ECD / C-ICD orientation')
        if membrane_class == 'tail_anchored':
            anchor = tm[0]
            length = max(ref_lengths.values()) if ref_lengths else None
            distance_to_terminus = anchor[0] - 1
            if length is not None:
                distance_to_terminus = min(distance_to_terminus,
                                           length - anchor[1])
            if distance_to_terminus > 40:
                warnings.append(
                    'tail_anchored target has a TM segment far from both '
                    'termini; check the membrane_class')
    elif membrane_class == 'beta_barrel':
        if not tm:
            raise ValueError('beta_barrel requires tm_segments (beta-strands)')
        if not (config.get('nterm') and config.get('cterm')):
            warnings.append(
                'beta_barrel target without explicit nterm/cterm; using '
                'first_loop parity (verify against an OPM/PDBTM annotation)')
    elif membrane_class == 'peripheral':
        if tm:
            warnings.append('peripheral target has tm_segments; they will be '
                            'treated as ordinary regions')
    else:
        if not tm:
            raise ValueError('alpha_helical_multi requires tm_segments')

    if len(tm) > 1:
        expected_nterm = 'ICD' if first_loop == 'ECD' else 'ECD'
        if config.get('nterm') and nterm != expected_nterm:
            warnings.append(
                f'nterm={nterm} is inconsistent with first_loop={first_loop} '
                f'(expected {expected_nterm})')

    _validate_segments(tm, ecd, icd, ref_lengths)

    signal = config.get('signal_peptide')
    signal_segments = []
    if signal:
        if isinstance(signal, (list, tuple)) and len(signal) == 2 \
                and not isinstance(signal[0], (list, tuple)):
            signal_segments = [(int(signal[0]), int(signal[1]))]
        else:
            signal_segments = _segments(signal, 'signal_peptide')

    return {
        'class': membrane_class,
        'tm_type': 'strand' if membrane_class == 'beta_barrel' else 'helix',
        'tm_segments': tm,
        'ecd_segments': ecd,
        'icd_segments': icd,
        'nterm': nterm,
        'cterm': cterm,
        'first_loop': first_loop,
        'signal_peptide': signal_segments,
        'warnings': warnings,
    }


def resolve_topologies(config, ref_lengths=None, chain_ids=None):
    """Resolve one topology per chain.

    A config may provide a global topology (v1) or per-chain overrides::

        topology:
          A: {tm_segments: [...], first_loop: ECD, ...}
          B: {tm_segments: [...], first_loop: ICD, ...}

    Chains without an override inherit the global keys, so existing configs
    keep working unchanged.
    """
    chains = list(chain_ids or config.get('chains') or [])
    if not chains:
        return {}
    overrides = config.get('topology')
    out = {}
    for cid in chains:
        local = dict(config)
        if isinstance(overrides, dict):
            override = overrides.get(cid)
            if isinstance(override, dict):
                local.update(override)
        lengths = None
        if ref_lengths:
            lengths = ({cid: ref_lengths[cid]} if cid in ref_lengths
                       else dict(ref_lengths))
        out[cid] = resolve_topology(local, lengths)
    return out


def region_at(refn, topology):
    """Classify a 1-based reference position as SIGNAL/TMD/ECD/ICD/LOOP."""
    for a, b in topology.get('signal_peptide') or []:
        if a <= refn <= b:
            return 'SIGNAL'
    for a, b in topology['tm_segments']:
        if a <= refn <= b:
            return 'TMD'
    for a, b in topology['ecd_segments']:
        if a <= refn <= b:
            return 'ECD'
    for a, b in topology['icd_segments']:
        if a <= refn <= b:
            return 'ICD'
    tm = topology['tm_segments']
    if not tm:
        return 'SOLUBLE'
    if refn < tm[0][0]:
        return topology['nterm']
    if refn > tm[-1][1]:
        return topology['cterm']
    return 'LOOP'


def region_map(topology, length):
    return OrderedDict((i, region_at(i, topology)) for i in range(1, length + 1))


def topology_summary(topology):
    regions = {'SIGNAL': 0, 'TMD': 0, 'ECD': 0, 'ICD': 0, 'LOOP': 0,
               'SOLUBLE': 0}
    for a, b in topology.get('signal_peptide') or []:
        regions['SIGNAL'] += b - a + 1
    for a, b in topology['tm_segments']:
        regions['TMD'] += b - a + 1
    for a, b in topology['ecd_segments']:
        regions['ECD'] += b - a + 1
    for a, b in topology['icd_segments']:
        regions['ICD'] += b - a + 1
    return {
        'class': topology['class'],
        'tm_type': topology['tm_type'],
        'nterm': topology['nterm'],
        'cterm': topology['cterm'],
        'first_loop': topology['first_loop'],
        'n_tm': len(topology['tm_segments']),
        'regions': {k: v for k, v in regions.items() if v},
        'warnings': list(topology.get('warnings') or []),
    }


# Backwards-compatible alias used by older code/tests.
classify_region = region_at
