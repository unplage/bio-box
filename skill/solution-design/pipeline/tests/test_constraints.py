"""Constraint builder tests: ECD, pockets, interface, fixed regions."""
import copy

from tmkit.constraints import build_constraints, detect_disulfides, \
    detect_ptm_sequons
from tmkit.mask import build_mask
from tmkit.topology import resolve_topology


def test_ptm_sequon_detection():
    assert detect_ptm_sequons('ANVST') == [2]
    assert detect_ptm_sequons('ANPST') == []


def test_disulfide_detection_returns_reference_pairs(target, cfg):
    pairs = detect_disulfides(target, cfg['chains'])
    assert isinstance(pairs, list)
    for ca, ra, cb, rb in pairs:
        assert ca in cfg['chains'] and cb in cfg['chains']
        assert isinstance(ra, int) and isinstance(rb, int)


def test_constraints_fix_ecd_and_pocket(target, cfg):
    topo = resolve_topology(cfg, {c: len(cfg['reference_sequences'][c])
                                  for c in cfg['chains']})
    con = build_constraints(target, cfg['chains'], cfg, topo)
    assert con['fix_ecd']
    assert 175 in con['fixed_reasons']['A']
    assert 268 in con['fixed_reasons']['A']
    assert 175 in con['pocket_shell']['A']


def test_interface_core_subset(target, cfg):
    conf = copy.deepcopy(cfg)
    conf['interface_mode'] = 'core'
    topo = resolve_topology(conf, {c: len(conf['reference_sequences'][c])
                                   for c in conf['chains']})
    con = build_constraints(target, conf['chains'], conf, topo)
    assert con['interface_mode'] == 'core'
    assert con['interface_core']
    assert con['interface_core'] <= set(con['interface_counts'])
    mask = build_mask(target, conf['chains'], conf)
    for cid in conf['chains']:
        for refn, info in mask[cid].items():
            if (cid, info['pdb_residue']) in con['interface_core']:
                assert not info['designable']


def test_fixed_region_is_not_designable(target, cfg):
    conf = copy.deepcopy(cfg)
    conf['fixed_regions'] = ['ICD']
    conf['design_scope'] = ['TMD', 'ICD']
    mask = build_mask(target, conf['chains'], conf)
    icd = [info for info in mask['A'].values() if info['region'] == 'ICD']
    assert icd
    assert all(not info['designable'] and info['reason'] == 'icd'
               for info in icd)


def test_pore_exclusion_without_fixing_seed(target, cfg):
    conf = copy.deepcopy(cfg)
    conf['design_scope'] = 'tm_core'
    conf['pore_residues'] = [80]
    conf['fix_pore'] = False
    conf['exclude_pore'] = True
    mask = build_mask(target, conf['chains'], conf)
    assert not mask['A'][80]['designable']
    assert mask['A'][80]['reason'] == 'pore'
