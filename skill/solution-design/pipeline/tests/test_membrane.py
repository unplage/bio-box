"""Membrane slab / lipid-facing mask tests."""
import numpy as np

from tmkit.mask import build_mask, summary
from tmkit.membrane import membrane_geometry, membrane_exclusion, \
    signed_distance


def test_membrane_geometry_explicit(target, cfg):
    conf = dict(cfg)
    conf['membrane'] = {'normal': [0, 0, 1], 'thickness': 30.0,
                        'headgroup': 15.0}
    geom = membrane_geometry(target, cfg['chains'], conf)
    assert geom is not None
    assert np.allclose(geom['normal'], [0, 0, 1])
    assert geom['core'] == 15.0 and geom['head'] == 15.0
    assert np.isfinite(geom['center'])


def test_membrane_geometry_auto_normal(target, cfg):
    conf = dict(cfg)
    conf['membrane_auto'] = True
    geom = membrane_geometry(target, cfg['chains'], conf)
    assert geom is not None
    assert np.isclose(np.linalg.norm(geom['normal']), 1.0)


def test_membrane_geometry_disabled(target, cfg):
    assert membrane_geometry(target, cfg['chains'], cfg) is None


def test_membrane_exclusion_reasons():
    geom = {'normal': np.array([0, 0, 1.0]), 'center': 0.0,
            'core': 15.0, 'head': 15.0,
            'exclude_aromatic_belt': True, 'restrict_to_core': False}
    assert membrane_exclusion(np.array([0, 0, 10.0]), geom, 'ALA') is None
    assert membrane_exclusion(np.array([0, 0, 40.0]), geom, 'ALA') == \
        'outside_membrane'
    assert membrane_exclusion(np.array([0, 0, 20.0]), geom, 'TRP') == \
        'aromatic_belt'
    assert membrane_exclusion(np.array([0, 0, 20.0]), geom, 'ALA') is None

    core_only = dict(geom, restrict_to_core=True)
    assert membrane_exclusion(np.array([0, 0, 20.0]), core_only, 'ALA') == \
        'headgroup'


def test_signed_distance():
    geom = {'normal': np.array([0, 0, 1.0]), 'center': 5.0}
    assert signed_distance(np.array([1, 2, 8.0]), geom) == 3.0


def test_membrane_slab_reduces_designable(target, cfg):
    base = summary(build_mask(target, cfg['chains'], cfg))
    conf = dict(cfg)
    conf['membrane'] = {'normal': [0, 0, 1], 'center': 1000.0,
                        'thickness': 20.0, 'headgroup': 10.0}
    masked = summary(build_mask(target, cfg['chains'], conf))
    assert base['designable'] > 0
    assert masked['designable'] == 0
    assert masked['reasons'].get('outside_membrane', 0) > 0
