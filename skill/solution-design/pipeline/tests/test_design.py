"""Design-input generation + design parsing tests."""
import json
from pathlib import Path

from tmkit.design import prepare_inputs, parse_design_fasta, mpnn_command


def test_prepare_inputs_writes_jsonls(target, cfg, tmp_path):
    from tmkit.mask import build_mask
    mask = build_mask(target, cfg['chains'], cfg)
    inputs = prepare_inputs(target, cfg['chains'], mask, cfg, tmp_path,
                            name=cfg['name'])
    assert inputs['pdb'].exists()
    fix = json.loads(inputs['fixed'].read_text())
    ties = json.loads(inputs['ties'].read_text())
    design = json.loads(inputs['design_chains'].read_text())
    key = cfg['name']
    # fixed + designable == modeled length per chain
    for c in cfg['chains']:
        n_des = sum(1 for info in mask[c].values() if info['designable'])
        assert len(fix[key][c]) == len(target.chains[c].residues) - n_des
        assert len(design[key][c]) == len(target.chains[c].residues)
    # C3 symmetry: single ties group across three chains
    assert ties[key]['0'][0]['A'] == ties[key]['0'][0]['B'] == \
        ties[key]['0'][0]['C']


def test_mpnn_command_flags(target, cfg, tmp_path):
    from tmkit.mask import build_mask
    from tmkit.design import prepare_inputs
    mask = build_mask(target, cfg['chains'], cfg)
    inputs = prepare_inputs(target, cfg['chains'], mask, cfg, tmp_path,
                            name=cfg['name'])
    cmd = mpnn_command(inputs, cfg, tmp_path, num_seqs=10,
                       temperature='0.1 0.2')
    assert '--use_soluble_model' in cmd
    assert '--omit_AAs' in cmd and 'C' in cmd
    assert '--fixed_positions_jsonl' in cmd
    assert '--tied_positions_jsonl' in cmd
    assert '--num_seq_per_target' in cmd and '10' in cmd


def test_parse_design_fasta_roundtrip(tmp_path):
    fa = tmp_path / 'seqs' / 'x.fa'
    fa.parent.mkdir(parents=True)
    fa.write_text(
        '>x,T=0.15,sample=1,score=0.95,seq_recovery=0.42\n'
        'AAAAAAAAAA/BBBBBBBBBB/CCCCCCCCCC\n')
    designs = parse_design_fasta(fa, ['A', 'B', 'C'])
    assert len(designs) == 1
    d = designs[0]
    assert d['chains']['A'] == 'AAAAAAAAAA'
    assert d['score'] == 0.95
    assert d['seq_recovery'] == 0.42
    assert d['temperature'] == '0.15'


def test_assemble_designs_unique_names(target, cfg, tmp_path):
    """Multiple MPNN samples share the PDB-name header prefix; assembled
    records must get unique names or prediction outputs would overwrite
    each other."""
    from tmkit.design import assemble_designs
    from tmkit.mask import build_mask
    mask = build_mask(target, cfg['chains'], cfg)
    from tmkit.pdbio import map_to_reference, AMINO3
    designs = []
    for i in range(3):
        chains_seq = {}
        for c in cfg['chains']:
            chain_seq = ''.join(
                AMINO3.get(r.resname, 'X')
                for r in target.chains[c].residues)
            chains_seq[c] = chain_seq
        designs.append({'name': '8UCD,T=0.15,sample=%d' % (i + 1),
                        'chains': chains_seq, 'score': 0.9 - i * 0.1,
                        'seq_recovery': 0.5, 'temperature': '0.15'})
    assembled = assemble_designs(designs, cfg, target, cfg['chains'],
                                 tmp_path, name=cfg['name'])
    names = [r['name'] for r in assembled]
    assert len(set(names)) == 3, names
    assert names[0].startswith('8UCD_0000')


def test_full_length_sequence_preserves_unmodeled(target, cfg):
    """Unmodeled regions (N/C termini, gap 143-149) stay at reference."""
    from tmkit.mask import build_mask
    from tmkit.design import assemble_designs
    from tmkit.pdbio import map_to_reference
    mask = build_mask(target, cfg['chains'], cfg)
    ref = cfg['reference_sequences']['A']
    # simulate one design: redesign all mask positions to 'S'
    design = {'name': 'd1', 'chains': {}, 'score': 0.9,
              'seq_recovery': 0.5, 'temperature': '0.15'}
    chains_full = {}
    for c in cfg['chains']:
        ordinal = {r.id[1]: i + 1 for i, r in
                   enumerate(target.chains[c].residues)}
        mapping = map_to_reference(target, c)
        des_refs = sorted(r for r, i in mask[c].items() if i['designable'])
        # chain sequence: keep modeled residues, flip designable -> S
        chain_seq = []
        for r in target.chains[c].residues:
            refn = mapping.get(r.id[1], r.id[1])
            aa = r.resname
            if refn in des_refs:
                chain_seq.append('S')
            else:
                from tmkit.pdbio import AMINO3
                chain_seq.append(AMINO3.get(aa, 'X'))
        design['chains'][c] = ''.join(chain_seq)
    rec = assemble_designs([design], cfg, target, cfg['chains'],
                           Path('/tmp'), name=cfg['name'])[0]
    full = rec['chains_full']['A']
    assert len(full) == len(ref)
    # unmodeled positions keep reference identity
    for pos in (1, 50, 66, 143, 144, 311, 339):
        assert full[pos - 1] == ref[pos - 1], pos
