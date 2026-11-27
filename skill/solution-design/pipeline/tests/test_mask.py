"""Mask + reference-numbering tests on the real 8UCD template."""
import numpy as np
import pytest

from tests.conftest import TM, ECD


def test_reference_mapping_identity(target, ref_seq):
    """PDB numbering equals UniProt numbering for 8UCD (verified 1:1)."""
    m = target.reference_map['A']
    assert m.get(67) == 67
    assert m.get(310) == 310
    assert m.get(268) == 268
    # unmodeled termini must not be mapped (chain only covers 67..310)
    assert 66 not in m.values()
    assert min(m.values()) == 67
    assert max(m.values()) == 310


def test_mask_tm_only(target, mask):
    """Designable residues are all inside TMD, none in ECD/ICD."""
    for cid in target.reference_map:
        for refn, info in mask[cid].items():
            if info['designable']:
                assert any(a <= refn <= b for a, b in TM), refn
            if info['region'] == 'ECD':
                assert not info['designable']


def test_heme_pocket_protected(target, mask):
    """H175/H268 (heme coordination) must never be designable."""
    for cid in target.reference_map:
        assert not mask[cid][175]['designable']
        assert not mask[cid][268]['designable']
        assert mask[cid][175]['reason'] == 'protected'
        assert mask[cid][268]['reason'] == 'protected'


def test_mask_symmetric(target, mask):
    """C3 symmetric mask: identical designable set across chains."""
    sets = [sorted(refn for refn, info in mask[c].items()
                   if info['designable'])
            for c in target.reference_map]
    assert sets[0] == sets[1] == sets[2]
    assert len(sets[0]) > 20  # meaningful surface, not degenerate


def test_interface_fixed(target, mask):
    """Oligomer interface residues are hard-fixed."""
    for cid in target.reference_map:
        for refn, info in mask[cid].items():
            if info['interface']:
                assert not info['designable']


def test_ecd_fully_covered(target, mask):
    """All modeled ECD residues are present and fixed."""
    ecd_refs = set()
    for a, b in ECD:
        ecd_refs |= set(range(a, b + 1))
    for cid in target.reference_map:
        fixed_ecd = [r for r, i in mask[cid].items() if i['region'] == 'ECD']
        assert fixed_ecd
        assert all(r in ecd_refs for r in fixed_ecd)


def test_topology_type1_gpcr():
    """7TM type-I GPCR: first loop intracellular, N-term ectodomain, C-term
    intracellular (e.g. a type-I GPCR)."""
    from tmkit.mask import auto_loop_topology, topology_sides
    tm = [[40, 62], [78, 100], [118, 140], [160, 182],
          [205, 227], [250, 272], [285, 307]]
    nterm, cterm, first = topology_sides(tm, {'first_loop': 'ICD'})
    assert (nterm, cterm, first) == ('ECD', 'ICD', 'ICD')
    ecd, icd = auto_loop_topology(tm, first)
    assert (63, 77) in icd and (101, 117) in ecd


def test_topology_type1_6tm():
    """6TM (e.g. STEAP1): first loop ECD => both termini intracellular."""
    from tmkit.mask import topology_sides
    tm = [[71, 91], [119, 139], [164, 184], [218, 238], [258, 278], [291, 311]]
    assert topology_sides(tm, {'first_loop': 'ECD'}) == ('ICD', 'ICD', 'ECD')


def test_tm_score_identity(target, ref_seq):
    """TM-score of a chain against itself is 1.0."""
    from tmkit.metrics import tm_score
    from tmkit.pdbio import ca_coords
    A, refs = ca_coords(target, 'A')
    assert np.isclose(tm_score(A, A), 1.0, atol=1e-6)


def test_per_chain_topology_override(cfg):
    from copy import deepcopy
    from tmkit.topology import resolve_topologies
    conf = deepcopy(cfg)
    conf['topology'] = {
        'B': {'tm_segments': [[70, 90]], 'ecd_segments': [],
              'icd_segments': []},
    }
    lengths = {c: len(conf['reference_sequences'][c]) for c in conf['chains']}
    tops = resolve_topologies(conf, lengths, conf['chains'])
    assert len(tops['A']['tm_segments']) == 6
    assert tops['B']['tm_segments'] == [(70, 90)]
    assert tops['C']['tm_segments'] == tops['A']['tm_segments']


def test_symmetry_groups_only_for_equivalence(cfg):
    from tmkit.constraints import symmetry_groups
    assert symmetry_groups(cfg, cfg['chains']) == [cfg['chains']]
    hetero = {'reference_sequences': {'A': 'AAA', 'B': 'CCC'}}
    assert symmetry_groups(hetero, ['A', 'B']) == []
    explicit = {'symmetry': {'groups': [['B', 'A']]}}
    assert symmetry_groups(explicit, ['A', 'B', 'C']) == [['B', 'A']]


def test_validate_config_rejects_bad_grouping(cfg):
    from copy import deepcopy
    from tmkit.config import normalize_config, validate_config
    conf = normalize_config(deepcopy(cfg))
    lengths = {c: len(conf['reference_sequences'][c]) for c in conf['chains']}
    assert isinstance(validate_config(conf, lengths), list)
    bad = deepcopy(conf)
    bad['symmetry']['groups'] = [['A', 'B'], ['B', 'C']]
    with pytest.raises(ValueError):
        validate_config(bad, lengths)
    bad = deepcopy(conf)
    bad['interface_pairs'] = [['A', 'A']]
    with pytest.raises(ValueError):
        validate_config(bad, lengths)


def test_cofactor_contacts_are_fixed(target, cfg):
    from copy import deepcopy
    from tmkit.mask import build_mask
    conf = deepcopy(cfg)
    conf['fixed_residues'] = []
    conf['cofactors'] = [
        {'name': 'heme', 'resnames': ['HEM'], 'shell': 4.5},
    ]
    result = build_mask(target, cfg['chains'], conf)
    contacts = [(c, r) for c in cfg['chains']
                for r, info in result[c].items()
                if info['reason'] == 'cofactor']
    assert contacts


def test_relative_sasa_map_and_mask(target, cfg, tmp_path):
    from copy import deepcopy
    from tmkit.pdbio import relative_sasa_map
    from tmkit.mask import build_mask
    rel = relative_sasa_map(target, cfg['chains'], workdir=str(tmp_path))
    assert rel and all(0.0 <= value for value in rel.values())
    conf = deepcopy(cfg)
    conf['relative_sasa'] = True
    conf['surface_rsa_threshold'] = 0.20
    conf['workdir'] = str(tmp_path)
    result = build_mask(target, cfg['chains'], conf)
    assert any(info['designable'] for c in cfg['chains']
               for info in result[c].values())


def test_topology_classes_and_peripheral_scope(cfg):
    from copy import deepcopy
    from tmkit.config import normalize_config
    from tmkit.topology import resolve_topology
    beta = resolve_topology({
        'membrane_class': 'beta_barrel',
        'tm_segments': [[1, 10], [21, 30]],
        'nterm': 'ECD', 'cterm': 'ECD',
    })
    assert beta['tm_type'] == 'strand'
    single = resolve_topology({
        'membrane_class': 'single_pass',
        'tm_segments': [[50, 70]], 'nterm': 'ECD', 'cterm': 'ICD',
    })
    assert single['tm_segments'] == [(50, 70)]
    peripheral = normalize_config(deepcopy({
        'chains': ['A'], 'membrane_class': 'peripheral',
        'reference_sequences': {'A': 'AA'},
    }))
    assert peripheral['design_scope'] == 'surface'


def test_membrane_geometries_are_per_chain(target, cfg):
    from copy import deepcopy
    from tmkit.membrane import membrane_geometries
    from tmkit.topology import resolve_topologies
    conf = deepcopy(cfg)
    conf['membrane'] = {
        'normal': 'auto', 'thickness': 30.0, 'headgroup': 15.0,
    }
    lengths = {c: len(conf['reference_sequences'][c]) for c in conf['chains']}
    tops = resolve_topologies(conf, lengths, conf['chains'])
    geoms = membrane_geometries(target, conf['chains'], tops, conf)
    assert set(geoms) == set(conf['chains'])
    assert all(g is not None for g in geoms.values())


def test_tagged_pdb_has_designable_atoms(target, cfg, mask, tmp_path):
    from tmkit.report import write_tagged_pdb
    p = tmp_path / 'tagged.pdb'
    write_tagged_pdb(target, cfg['chains'], mask, p)
    lines = [ln for ln in p.read_text().splitlines() if ln.startswith('ATOM')]
    assert len(lines) > 100
    assert any(float(ln[60:66]) == 1.0 for ln in lines)
    assert any(float(ln[60:66]) == 0.0 for ln in lines)
