"""Evaluation robustness tests: alignment-driven pairing, P0 fix, gates."""
import numpy as np
import pytest

from tmkit.predict import find_rank1_models, evaluate_models, multimer_fasta
from tmkit.report import passes_filters
from tmkit.stability import sequence_solubility, designable_positions


def make_model_like(target, cfg, outdir, resnum_start=1):
    """Write a synthetic 'predicted' PDB using the template CA coords but
    with renumbered residues (resnum_start..) per chain, mimicking ColabFold
    output numbering that differs from the template."""
    from Bio.PDB.PDBIO import PDBIO
    struct = target.structure
    for model in struct:
        for chain in model:
            if chain.id not in cfg['chains']:
                continue
            new_id = resnum_start
            for residue in list(chain):
                chain.detach_child(residue.id)
                residue.id = (' ', new_id, ' ')
                chain.add(residue)
                new_id += 1
    path = outdir / 'renumbered.pdb'
    io = PDBIO()
    io.set_structure(struct)
    io.save(str(path))
    return path


def test_find_rank1_models(tmp_path):
    (tmp_path / 'd1_rank_001_alphafold2_multimer_v3_model_1_seed_000.pdb') \
        .write_text('')
    (tmp_path / 'd2_pred_0.pdb').write_text('')
    found = find_rank1_models(tmp_path, ['d1', 'd2', 'd3'])
    assert 'd1' in found and 'd2' in found and 'd3' not in found


def test_evaluate_models_alignment_driven(target, cfg, tmp_path):
    """ECD RMSD must be near 0 even when the model chain is renumbered
    (1..N) -- pairing is alignment-driven, not number-driven."""
    from tmkit.mask import build_mask
    mask = build_mask(target, cfg['chains'], cfg)
    ref = cfg['reference_sequences']['A']
    model = make_model_like(target, cfg, tmp_path, resnum_start=1)
    rec = {
        'name': 'd1', 'chains_full': {c: ref for c in cfg['chains']},
        'designed_positions': {}, 'score': 0.9, 'seq_recovery': 1.0,
        'temperature': '0.15',
    }
    rows = evaluate_models({'d1': model}, [rec], cfg['chains'], target,
                           mask, cfg['ecd_segments'], {c: ref for c in cfg['chains']},
                           tmp_path / 'pred',
                           pocket_residues=cfg.get('pocket_residues'))
    r = rows[0]
    assert r['ecd_rmsd'] < 0.01, r['ecd_rmsd']
    assert r['tm_score'] > 0.95, r['tm_score']
    assert r['pocket_rmsd'] < 0.01, r['pocket_rmsd']
    assert r['cys_des'] == 0


def test_evaluate_models_shifted_renumbering(target, cfg, tmp_path):
    """Robust even when model numbering is offset relative to template."""
    from tmkit.mask import build_mask
    mask = build_mask(target, cfg['chains'], cfg)
    ref = cfg['reference_sequences']['A']
    model = make_model_like(target, cfg, tmp_path, resnum_start=100)
    rec = {
        'name': 'd2', 'chains_full': {c: ref for c in cfg['chains']},
        'designed_positions': {}, 'score': 0.9, 'seq_recovery': 1.0,
        'temperature': '0.15',
    }
    rows = evaluate_models({'d2': model}, [rec], cfg['chains'], target,
                           mask, cfg['ecd_segments'], {c: ref for c in cfg['chains']},
                           tmp_path / 'pred',
                           pocket_residues=cfg.get('pocket_residues'))
    assert rows[0]['ecd_rmsd'] < 0.01


def test_evaluate_models_missing_model(target, cfg, tmp_path):
    from tmkit.mask import build_mask
    mask = build_mask(target, cfg['chains'], cfg)
    ref = cfg['reference_sequences']['A']
    rec = {'name': 'd1', 'chains_full': {c: ref for c in cfg['chains']},
           'designed_positions': {}}
    rows = evaluate_models({}, [rec], cfg['chains'], target, mask,
                           cfg['ecd_segments'],
                           {c: ref for c in cfg['chains']}, tmp_path / 'pred')
    r = rows[0]
    assert np.isnan(r['tm_score'])
    assert r['model'] == ''


def test_filters_record_failures(cfg):
    row = {'name': 'x', 'tm_score': 0.5, 'plddt': 90.0}
    ok = passes_filters(row, cfg['filters'])
    assert not ok
    assert 'tm_score' in row['fail_gates']
    assert 'plddt' not in row['fail_gates']


def test_sequence_solubility_mixed(target, cfg):
    from tmkit.mask import build_mask
    mask = build_mask(target, cfg['chains'], cfg)
    des = designable_positions(mask, 'A')
    ref = cfg['reference_sequences']['A']
    seq = list(ref)
    for p in des:
        seq[p - 1] = 'S'  # polarize everything
    stats = sequence_solubility(''.join(seq), ref, des)
    assert stats['soluble_score'] > 0.8
    assert stats['residual_hydro'] == 0.0


def test_surf_hydro_des_metric(target, cfg, tmp_path):
    """surf_hydro_des is restricted to designable positions only."""
    from tmkit.mask import build_mask
    from tmkit.metrics import surface_hydrophobic_fraction_at
    mask = build_mask(target, cfg['chains'], cfg)
    ref = cfg['reference_sequences']['A']
    model = make_model_like(target, cfg, tmp_path, resnum_start=1)
    from tmkit.pdbio import map_model_to_reference, residue_sasa_from_file
    m_map = map_model_to_reference(model, 'A', ref)
    refn_to_model = {v: k for k, v in m_map.items()}
    sasa = residue_sasa_from_file(model)
    aa_map = {}
    with open(model) as f:
        for line in f:
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                aa_map[(line[21], int(line[22:26]))] = 'A'
    des_pos = [('A', refn_to_model[r]) for r, i in mask['A'].items()
               if i['designable'] and r in refn_to_model]
    frac = surface_hydrophobic_fraction_at(sasa, aa_map, des_pos)
    assert 0.0 <= frac <= 1.0
