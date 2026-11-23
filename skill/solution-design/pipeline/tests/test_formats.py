"""External-model formats: mmCIF input, predictor naming, confidence JSON,
and single-chain (monomer) generality."""
import json
import shutil

import numpy as np
import pytest

from tmkit.pdbio import as_pdb, parse_pdb, harmonize_chains, reference_numbering
from tmkit.mask import build_mask, summary
from tmkit.predict import find_rank1_models, evaluate_models
from tmkit.confidence import load_confidence

from Bio.PDB.mmcifio import MMCIFIO


@pytest.fixture()
def template_cif(target, tmp_path):
    path = tmp_path / 'template.cif'
    io = MMCIFIO()
    io.set_structure(target.structure)
    io.save(str(path))
    return path


def test_mmcif_template_parses(template_cif, cfg):
    assert as_pdb(template_cif).exists()
    t = parse_pdb(template_cif, chain_ids=cfg['chains'])
    harmonize_chains(t, cfg['chains'])
    for c in cfg['chains']:
        reference_numbering(t, c, cfg['reference_sequences'][c])
    assert list(t.chains) == cfg['chains']
    assert len(t.chains['A'].residues) > 200


def test_find_models_boltz_and_af3(template_cif, tmp_path):
    boltz = tmp_path / 'boltz_results_job' / 'predictions' / 'd1'
    boltz.mkdir(parents=True)
    shutil.copy(template_cif, boltz / 'd1_model_0.cif')
    af3 = tmp_path / 'af3_out' / 'd2'
    af3.mkdir(parents=True)
    shutil.copy(template_cif, af3 / 'd2_model.cif')
    found = find_rank1_models(tmp_path, ['d1', 'd2'])
    assert found['d1'].name == 'd1_model_0.cif'
    assert found['d2'].name == 'd2_model.cif'


def test_confidence_per_design_not_crossing(tmp_path):
    for name, iptm in (('d1', 0.91), ('d2', 0.55)):
        (tmp_path / f'{name}_model_0.cif').write_text('')
        (tmp_path / f'confidence_{name}_model_0.json').write_text(
            json.dumps({'iptm': iptm, 'ptm': 0.8}))
    c1 = load_confidence(tmp_path / 'd1_model_0.cif')
    c2 = load_confidence(tmp_path / 'd2_model_0.cif')
    assert c1['iptm'] == pytest.approx(0.91)
    assert c2['iptm'] == pytest.approx(0.55)


def test_af3_summary_confidences(tmp_path):
    folder = tmp_path / 'job' / 'd1'
    folder.mkdir(parents=True)
    (folder / 'd1_model.cif').write_text('')
    (folder / 'summary_confidences.json').write_text(
        json.dumps({'iptm': 0.77, 'ptm': 0.81, 'ranking_score': 0.7}))
    conf = load_confidence(folder / 'd1_model.cif')
    assert conf['iptm'] == pytest.approx(0.77)
    assert conf['ptm'] == pytest.approx(0.81)
    assert conf['ranking_score'] == pytest.approx(0.7)


def test_evaluate_external_cif_reports_iptm_and_tm_hydro(
        target, cfg, template_cif, tmp_path):
    mask = build_mask(target, cfg['chains'], cfg)
    ref = {c: cfg['reference_sequences'][c] for c in cfg['chains']}
    model = tmp_path / 'd1_model_0.cif'
    shutil.copy(template_cif, model)
    (tmp_path / 'confidence_d1_model_0.json').write_text(
        json.dumps({'iptm': 0.85, 'ptm': 0.9}))
    rec = {'name': 'd1', 'chains_full': ref, 'designed_positions': {},
           'score': 0.9, 'seq_recovery': 1.0, 'temperature': '0.15'}
    rows = evaluate_models({'d1': model}, [rec], cfg['chains'], target, mask,
                           cfg['ecd_segments'], ref, tmp_path,
                           pocket_residues=cfg.get('pocket_residues'),
                           tm_segments=cfg.get('tm_segments'))
    row = rows[0]
    assert row['iptm'] == pytest.approx(0.85)
    assert row['ptm'] == pytest.approx(0.9)
    assert row['tm_score'] > 0.99
    assert row['ecd_rmsd'] < 0.01
    assert 0.0 <= row['surf_hydro_tm'] <= 1.0
    assert np.isfinite(row['chain_symmetry_rmsd'])
    assert row['interface_recovery'] == pytest.approx(1.0)


def test_interface_recovery_drops_when_chain_separated(
        target, cfg, tmp_path):
    """Shifting one chain away removes its interface contacts."""
    from Bio.PDB.PDBIO import PDBIO
    from tmkit.pdbio import parse_pdb
    shifted = parse_pdb(target.path, chain_ids=cfg['chains'])
    for residue in shifted.chains['C'].residues:
        for atom in residue:
            atom.coord = atom.coord + np.array([100.0, 0.0, 0.0])
    model = tmp_path / 'd1_model_0.pdb'
    io = PDBIO()
    io.set_structure(shifted.structure)
    io.save(str(model))

    mask = build_mask(target, cfg['chains'], cfg)
    ref = {c: cfg['reference_sequences'][c] for c in cfg['chains']}
    rec = {'name': 'd1', 'chains_full': ref, 'designed_positions': {},
           'score': 0.9, 'seq_recovery': 1.0, 'temperature': '0.15'}
    row = evaluate_models({'d1': model}, [rec], cfg['chains'], target, mask,
                          cfg['ecd_segments'], ref, tmp_path)[0]
    assert 0.0 <= row['interface_recovery'] < 0.5


def test_external_input_emission(cfg, tmp_path):
    from tmkit.predict import write_external_inputs
    ref = {c: cfg['reference_sequences'][c] for c in cfg['chains']}
    rec = {'name': 'd1', 'chains_full': ref}
    write_external_inputs([rec], cfg['chains'], tmp_path)
    boltz = (tmp_path / 'boltz' / 'd1.yaml').read_text()
    af3 = json.loads((tmp_path / 'af3' / 'd1.json').read_text())
    assert 'id: A' in boltz and 'sequence:' in boltz
    assert boltz.index('id: A') < boltz.index('id: B') < boltz.index('id: C')
    assert [s['protein']['id'] for s in af3['sequences']] == cfg['chains']
    assert af3['name'] == 'd1'


def test_monomer_mask(target, cfg):
    single = dict(cfg)
    single['chains'] = ['A']
    single['symmetric_mask'] = False
    t = parse_pdb(target.path, chain_ids=['A'])
    harmonize_chains(t, ['A'])
    reference_numbering(t, 'A', cfg['reference_sequences']['A'])
    mask = build_mask(t, ['A'], single)
    info = summary(mask)
    assert info['designable'] > 0
    # no cross-chain interface residues can exist for a single chain
    assert all(not m['interface'] for m in mask['A'].values())
