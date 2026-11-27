"""Generic topology model tests (no STEAP1-specific assumptions)."""
from tmkit.config import validate_config
from tmkit.topology import auto_loop_topology, resolve_topology, region_at, \
    topology_sides


def test_six_tm_first_loop_ecd():
    tm = [[71, 91], [119, 139], [164, 184],
          [218, 238], [258, 278], [291, 311]]
    assert topology_sides(tm, {'first_loop': 'ECD'}) == ('ICD', 'ICD', 'ECD')
    ecd, icd = auto_loop_topology(tm, 'ECD')
    assert (92, 118) in ecd and (140, 163) in icd


def test_type1_gpcr_first_loop_icd():
    tm = [[40, 62], [78, 100], [118, 140], [160, 182],
          [205, 227], [250, 272], [285, 307]]
    assert topology_sides(tm, {'first_loop': 'ICD'}) == \
        ('ECD', 'ICD', 'ICD')


def test_single_pass_explicit_sides_no_warning():
    cfg = {
        'membrane_class': 'single_pass',
        'tm_segments': [[100, 120]],
        'nterm': 'ECD', 'cterm': 'ICD',
        'reference_sequences': {'A': 'M' * 200},
    }
    topo = resolve_topology(cfg, {'A': 200})
    assert topo['nterm'] == 'ECD' and topo['cterm'] == 'ICD'
    assert not topo['warnings']


def test_single_pass_inferred_orientation_warns():
    cfg = {
        'membrane_class': 'single_pass',
        'tm_segments': [[100, 120]],
        'reference_sequences': {'A': 'M' * 200},
    }
    topo = resolve_topology(cfg, {'A': 200})
    assert topo['nterm'] == 'ECD' and topo['cterm'] == 'ICD'
    assert any('explicit nterm/cterm' in w for w in topo['warnings'])


def test_beta_barrel_without_pore_is_warned():
    cfg = {
        'name': 'barrel', 'chains': ['A'],
        'membrane_class': 'beta_barrel',
        'tm_segments': [[10, 18], [25, 33], [40, 48]],
        'first_loop': 'ECD',
        'reference_sequences': {'A': 'M' * 100},
    }
    topo = resolve_topology(cfg, {'A': 100})
    assert topo['tm_type'] == 'strand'
    warnings = validate_config(cfg, {'A': 100})
    assert any('pore' in w.lower() for w in warnings)


def test_region_at_uses_terminal_sides():
    topo = {
        'tm_segments': [(10, 20)],
        'ecd_segments': [],
        'icd_segments': [],
        'nterm': 'ECD', 'cterm': 'ICD',
        'signal_peptide': [],
    }
    assert region_at(5, topo) == 'ECD'
    assert region_at(15, topo) == 'TMD'
    assert region_at(30, topo) == 'ICD'
