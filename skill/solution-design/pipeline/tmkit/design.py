import json
import subprocess
import sys
from pathlib import Path

from .pdbio import write_structure, write_fasta, full_length_sequence, \
    mpnn_positions


def prepare_inputs(target, chain_ids, mask, config, outdir, name='target'):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    clean_pdb = outdir / f'{name}.pdb'
    write_structure(target, clean_pdb, chains=chain_ids,
                    omit_resnames=config.get('omit_resnames', []),
                    omit_resseq=config.get('omit_resseq', []))
    fix = fixed_positions_dict(target, chain_ids, mask)
    ties = symmetry_ties_dict(target, chain_ids, mask)
    design_chains = {}
    for c in chain_ids:
        positions = mpnn_positions(target, c)
        n = max(positions.values()) if positions else 0
        design_chains[c] = list(range(1, n + 1))
    fix_name = outdir / f'{name}_fixed.jsonl'
    ties_name = outdir / f'{name}_ties.jsonl'
    des_name = outdir / f'{name}_designchains.jsonl'
    fix_name.write_text(json.dumps({name: fix}) + '\n')
    ties_name.write_text(json.dumps({name: ties}) + '\n')
    des_name.write_text(json.dumps({name: design_chains}) + '\n')
    return {
        'pdb': clean_pdb, 'fixed': fix_name, 'ties': ties_name,
        'design_chains': des_name, 'name': name,
    }


def fixed_positions_dict(target, chain_ids, mask):
    from .mask import resolve
    _, fixed = resolve(target, chain_ids, mask)
    return fixed


def symmetry_ties_dict(target, chain_ids, mask):
    from .mask import symmetry_ties_dict as _t
    return _t(target, chain_ids, mask)


def mpnn_command(inputs, config, outdir, num_seqs=8, temperature='0.15 0.2'):
    prog = Path(config.get('proteinmpnn_dir', 'proteinmpnn'))
    cmd = [
        sys.executable, str(prog / 'protein_mpnn_run.py'),
        '--pdb_path', str(inputs['pdb']),
        '--fixed_positions_jsonl', str(inputs['fixed']),
        '--tied_positions_jsonl', str(inputs['ties']),
        '--chain_id_jsonl', str(inputs['design_chains']),
        '--num_seq_per_target', str(num_seqs),
        '--sampling_temp', temperature,
        '--out_folder', str(outdir),
        '--omit_AAs', 'C',
        '--suppress_print', '0',
    ]
    if config.get('soluble_model'):
        cmd.append('--use_soluble_model')
    model = config.get('model_name')
    if model:
        cmd += ['--model_name', model]
    backbone_noise = config.get('backbone_noise')
    if backbone_noise is not None:
        cmd += ['--backbone_noise', str(backbone_noise)]
    bias = config.get('bias_aa')
    if bias:
        bias_json = Path(outdir) / 'bias.jsonl'
        bias_json.write_text(json.dumps({inputs['name']: bias}) + '\n')
        cmd += ['--bias_AA_jsonl', str(bias_json)]
    return cmd


def parse_design_fasta(fa_path, chain_ids):
    records = []
    current = None
    for line in Path(fa_path).read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            current = line[1:]
            records.append([current, []])
        else:
            records[-1][1].append(line)
    designs = []
    for header, chunks in records:
        seq = ''.join(chunks)
        if '/' not in seq:
            continue
        chains = seq.split('/')
        if len(chains) != len(chain_ids):
            continue
        attr = {}
        for token in header.split(','):
            pair = token.strip().split('=')
            if len(pair) == 2:
                attr[pair[0].strip()] = pair[1].strip()
        # ProteinMPNN writes a leading non-sampled record (input sequence,
        # no T=/sample= attrs); skip it so it is not treated as a design.
        if 'T' not in attr and 'sample' not in attr:
            continue
        designs.append({
            'name': header.split(',')[0].strip(),
            'chains': {c: s for c, s in zip(chain_ids, chains)},
            'score': float(attr.get('score', 'nan')),
            'seq_recovery': float(attr.get('seq_recovery', 'nan')),
            'temperature': attr.get('T'),
        })
    return designs


def assemble_designs(designs, config, target, chain_ids, outdir, name='target'):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    assembled = []
    seen = {}
    for d in designs:
        base = d['name'].split(',')[0].strip() or name
        n = seen.get(base, 0)
        seen[base] = n + 1
        uname = f'{base}_{n:04d}'
        rec = {
            'name': uname, 'score': d['score'],
            'seq_recovery': d['seq_recovery'], 'temperature': d['temperature'],
            'chains_full': {},
            'designed_positions': {},
        }
        for cid in chain_ids:
            ref = config['reference_sequences'][cid]
            refnums = []
            seq_by_ref = {}
            from .pdbio import map_to_reference
            mapping = map_to_reference(target, cid)
            positions = mpnn_positions(target, cid)
            chain_seq = d['chains'][cid]
            if len(chain_seq) == len(target.chains[cid].residues):
                positions = {r.id[1]: i + 1 for i, r in
                             enumerate(target.chains[cid].residues)}
            for r in target.chains[cid].residues:
                pdb_res = r.id[1]
                refn = mapping.get(pdb_res, pdb_res)
                pos = positions.get(pdb_res)
                if pos is None or pos > len(chain_seq):
                    continue
                seq_by_ref[refn] = chain_seq[pos - 1]
            rec['designed_positions'][cid] = sorted(seq_by_ref)
            full = full_length_sequence(ref, sorted(seq_by_ref),
                                        [seq_by_ref[r] for r in sorted(seq_by_ref)])
            rec['chains_full'][cid] = full
        assembled.append(rec)
    return assembled