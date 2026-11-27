"""Config normalization, validation and scaffold-backbone discovery."""
from tmkit.config import normalize_config, resolve_template_path, validate_config
from tmkit.design import list_backbones, rfdiffusion_command


def _base():
    return {
        'name': 'x', 'chains': ['A'],
        'reference_sequences': {'A': 'M' * 120},
        'tm_segments': [[30, 50]],
        'membrane_class': 'single_pass',
        'nterm': 'ECD', 'cterm': 'ICD',
    }


def test_normalize_defaults_are_generic():
    cfg = normalize_config(_base())
    assert cfg['fix_ecd'] is True
    assert cfg['interface_mode'] == 'fixed'
    assert cfg['design_scope'] == 'tm_surface'
    assert cfg['strict_predict'] is True
    assert cfg['assembly'] == 'as_is'


def test_validate_warns_without_membrane_slab():
    cfg = _base()
    warnings = validate_config(cfg, {'A': 120})
    assert any('membrane slab' in w for w in warnings)


def test_assembly_path_overrides_template(tmp_path):
    p = tmp_path / 'assembly.pdb'
    p.write_text('END\n')
    cfg = _base()
    cfg['assembly'] = str(p)
    assert resolve_template_path(cfg, 'template.pdb') == str(p)


def test_list_backbones_and_rfdiffusion_command(tmp_path):
    (tmp_path / 'b.pdb').write_text('END\n')
    (tmp_path / 'a.pdb').write_text('END\n')
    (tmp_path / 'ignore.txt').write_text('')
    files = list_backbones({'backbone_dir': str(tmp_path)})
    assert [p.name for p in files] == ['a.pdb', 'b.pdb']

    cmd = rfdiffusion_command(
        {'template_pdb': 't.pdb', 'rfdiffusion_script': '/opt/rfd/run.py',
         'rfdiffusion_contigs': '[A1-120]', 'rfdiffusion_num_designs': 3},
        tmp_path / 'out', name='bb')
    assert cmd[0].endswith('python') or 'python' in cmd[0]
    assert 'inference.input_pdb=t.pdb' in cmd
    assert 'contigmap.contigs=[A1-120]' in cmd
    assert 'inference.num_designs=3' in cmd


def test_cli_overrides_are_ephemeral():
    from types import SimpleNamespace
    from tmkit.pipeline import apply_cli_overrides
    args = SimpleNamespace(proteinmpnn_dir='/fake/mpnn',
                           proteinmpnn_python=None, design_mode=None,
                           predict_tool='boltz',
                           rfdiffusion_script=None, rfdiffusion_dir=None,
                           predict_executable=None, colabfold_executable=None,
                           boltz_executable=None, boltz_devices=None,
                           boltz_use_msa_server=None, boltz_no_kernels=None,
                           require_predictions=None, strict_predict=None,
                           alphafold_script=None)
    original = _base()
    out = apply_cli_overrides(original, args)
    assert out['proteinmpnn_dir'] == '/fake/mpnn'
    assert out['predict_tool'] == 'boltz'
    assert 'proteinmpnn_dir' not in original
