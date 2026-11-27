"""Predictor command construction: ColabFold, AlphaFold2 and Boltz-2."""
from pathlib import Path

from tmkit.predict import predict_command


def test_colabfold_uses_directory_for_many_fastas(tmp_path):
    paths = [tmp_path / 'a.fa', tmp_path / 'b.fa']
    cmd = predict_command({'predict_tool': 'colabfold'}, paths, tmp_path / 'out')
    assert isinstance(cmd, list)
    assert cmd[0] == 'colabfold_batch'
    assert cmd[1] == str(tmp_path)
    assert ',' not in cmd[1]


def test_colabfold_single_fasta_is_passed_as_file(tmp_path):
    fa = tmp_path / 'a.fa'
    cmd = predict_command({'predict_tool': 'colabfold'}, [fa], tmp_path / 'out')
    assert isinstance(cmd, list)
    assert cmd[1] == str(fa)


def test_alphafold_fasta_paths_are_space_separated(tmp_path):
    paths = [tmp_path / 'a.fa', tmp_path / 'b.fa']
    cmd = predict_command(
        {'predict_tool': 'alphafold', 'alphafold_script': '/opt/af/run.py',
         'data_dir': '/data'}, paths, tmp_path / 'out')
    assert isinstance(cmd, list)
    idx = cmd.index('--fasta_paths')
    assert cmd[idx + 1] == f'{paths[0]} {paths[1]}'
    assert '--max_template_date' in cmd


def test_boltz_command_uses_yaml_dir(tmp_path):
    paths = [tmp_path / 'a.fa']
    cmd = predict_command(
        {'predict_tool': 'boltz', 'boltz_use_msa_server': True,
         'boltz_no_kernels': True},
        paths, tmp_path / 'out')
    assert isinstance(cmd, list)
    assert cmd[:3] == ['boltz', 'predict', str(tmp_path / 'boltz')]
    assert '--use_msa_server' in cmd
    assert '--out_dir' in cmd
    assert '--no_kernels' in cmd


def test_external_predictor_has_no_command(tmp_path):
    assert predict_command({'predict_tool': 'af3'},
                           [tmp_path / 'a.fa'], tmp_path / 'out') is None


def test_find_rank1_models_af2_layout(tmp_path):
    from tmkit.predict import find_rank1_models
    d = tmp_path / 'design_0001'
    d.mkdir()
    (d / 'ranked_0.pdb').write_text('END\n')
    found = find_rank1_models(tmp_path, ['design_0001'])
    assert found['design_0001'].name == 'ranked_0.pdb'


def test_executable_overrides(tmp_path):
    cmd = predict_command(
        {'predict_tool': 'boltz', 'boltz_executable': '/fake/boltz'},
        [tmp_path / 'a.fa'], tmp_path / 'out')
    assert isinstance(cmd, list) and cmd[0] == '/fake/boltz'
    cmd = predict_command(
        {'predict_tool': 'colabfold',
         'colabfold_executable': '/fake/colabfold_batch'},
        [tmp_path / 'a.fa'], tmp_path / 'out')
    assert isinstance(cmd, list) and cmd[0] == '/fake/colabfold_batch'
