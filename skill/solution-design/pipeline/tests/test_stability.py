"""Solubility / stability heuristic metrics."""
import math

from tmkit.stability import net_charge, isoelectric_point, \
    hydrophobic_moment, max_hydrophobic_window, max_charge_run, \
    max_polar_window, sequence_properties


def test_net_charge_ph_dependence():
    seq = 'K' * 10 + 'A' * 10
    assert net_charge(seq, ph=7.0) > 8.0
    assert net_charge(seq, ph=13.0) < net_charge(seq, ph=7.0)
    acidic = 'D' * 10 + 'A' * 10
    assert net_charge(acidic, ph=7.0) < -8.0


def test_isoelectric_point_range():
    pI = isoelectric_point('A' * 50)
    assert 0.0 < pI < 14.0
    assert isoelectric_point('K' * 20) > isoelectric_point('A' * 20)
    assert isoelectric_point('D' * 20) < isoelectric_point('A' * 20)


def test_hydrophobic_moment_helix():
    hydro = hydrophobic_moment('L' * 40, tm_segments=[(1, 40)])
    polar = hydrophobic_moment('S' * 40, tm_segments=[(1, 40)])
    assert hydro > polar
    assert math.isfinite(hydro)


def test_max_hydrophobic_window_patchiness():
    hydro = max_hydrophobic_window('ILFVWM' * 3, list(range(1, 16)), window=9)
    polar = max_hydrophobic_window('S' * 20, list(range(1, 16)), window=9)
    assert hydro == 1.0
    assert polar == 0.0


def test_sequence_properties_keys(cfg):
    ref = cfg['reference_sequences']['A']
    props = sequence_properties(ref, list(range(1, 40)),
                                cfg['tm_segments'])
    for key in ('net_charge', 'abs_net_charge', 'isoelectric_point',
                'hydrophobic_moment', 'max_hydro_window',
                'max_charge_run', 'max_polar_window'):
        assert key in props and math.isfinite(props[key])


def test_step_stability_adds_properties(target, cfg):
    from tmkit.pipeline import step_stability
    ref = {c: cfg['reference_sequences'][c] for c in cfg['chains']}
    row = {'name': 'd1', 'chains_full': ref}
    out = step_stability([row], target, cfg)[0]
    for key in ('net_charge', 'abs_net_charge', 'isoelectric_point',
                'hydrophobic_moment', 'max_hydro_window',
                'max_charge_run', 'max_polar_window'):
        assert key in out and math.isfinite(out[key])
