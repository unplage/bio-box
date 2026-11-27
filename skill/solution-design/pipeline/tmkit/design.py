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
    ties = symmetry_ties_dict(target, chain_ids, mask, config)
    # ProteinMPNN's --chain_id_jsonl format is {pdb: [masked, visible]}, not
    # {chain: [positions]}.  The old format was silently ignored because
    # --pdb_path overrides chain_id_dict; keep it correct for --jsonl_path.
    designed = [c for c in (config.get('mpnn_designed_chains') or chain_ids)
                if c in chain_ids]
    visible = [c for c in chain_ids if c not in designed]
    design_chains = [designed, visible]
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


def symmetry_ties_dict(target, chain_ids, mask, config=None):
    from .mask import symmetry_ties_dict as _t
    return _t(target, chain_ids, mask, config)


def list_backbones(config):
    """List user-supplied RFdiffusion/scaffold backbones, sorted by name."""
    bdir = config.get('backbone_dir')
    if not bdir:
        return []
    base = Path(bdir)
    if not base.exists():
        raise FileNotFoundError(f'backbone_dir does not exist: {base}')
    globs = config.get('backbone_globs') or ['*.pdb', '*.cif', '*.mmcif']
    files = sorted({p for pattern in globs for p in base.rglob(pattern)})
    limit = config.get('max_backbones')
    return files[:int(limit)] if limit else files


def rfdiffusion_script(config):
    script = config.get('rfdiffusion_script')
    if script:
        return Path(script)
    root = Path(config.get('rfdiffusion_dir', 'RFdiffusion'))
    return root / 'scripts' / 'run_inference.py'


def rfdiffusion_command(config, outdir, name='backbone'):
    """Build an RFdiffusion ``run_inference.py`` command.

    The adapter is intentionally thin: RFdiffusion options vary by version,
    so ``rfdiffusion_extra`` can carry any additional Hydra overrides.  The
    produced backbones are then sequenced by ProteinMPNN in ``step_design``.
    """
    script = rfdiffusion_script(config)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    template = config.get('template_pdb') or config.get('pdb')
    if not template:
        raise ValueError('RFdiffusion needs template_pdb/pdb in the config')
    cmd = [
        sys.executable, str(script),
        f'inference.input_pdb={template}',
        f'inference.output_prefix={outdir / name}',
        f'inference.num_designs={int(config.get("rfdiffusion_num_designs", 8))}',
    ]
    if config.get('rfdiffusion_contigs'):
        cmd.append(f'contigmap.contigs={config["rfdiffusion_contigs"]}')
    for key, hydra in [('rfdiffusion_noise_ca', 'denoiser.noise_scale_ca'),
                       ('rfdiffusion_noise_frame',
                        'denoiser.noise_scale_frame')]:
        if config.get(key) is not None:
            cmd.append(f'{hydra}={config[key]}')
    cmd += [str(x) for x in (config.get('rfdiffusion_extra') or [])]
    return cmd


def prepare_backbones(config, outdir, dry_run=False):
    """Return scaffold backbones, optionally running RFdiffusion first."""
    backbones = list_backbones(config)
    if backbones:
        return backbones
    mode = str(config.get('design_mode', 'sequence')).lower()
    if mode not in ('rf_diffusion', 'rfdiffusion', 'scaffold',
                    'partial_diffusion'):
        return []
    script = rfdiffusion_script(config)
    if not dry_run and not script.exists():
        raise FileNotFoundError(
            f'RFdiffusion script not found: {script}; set rfdiffusion_script '
            f'or rfdiffusion_dir')
    cmd = rfdiffusion_command(config, Path(outdir) / 'backbones',
                              name=config.get('name', 'target'))
    print('RFdiffusion command:')
    print('  ' + ' \\\n  '.join(cmd))
    if dry_run:
        return []
    proc = subprocess.run(cmd)
    if proc.returncode != 0:
        raise RuntimeError(
            f'RFdiffusion failed with exit code {proc.returncode}')
    return list_backbones({**config, 'backbone_dir': Path(outdir) / 'backbones'})


def _as_list(value):
    if value is None:
        return [None]
    if isinstance(value, (list, tuple)):
        return list(value) if value else [None]
    return [value]


def mpnn_runs(config):
    """Expand scalar or list model/noise/seed settings into individual runs.

    Each run gets its own output folder so FASTA files never overwrite one
    another. ``seeds`` falls back to ``seed`` for a single explicit seed.
    """
    models = _as_list(config.get('model_name'))
    noises = _as_list(config.get('backbone_noise'))
    seeds = _as_list(config.get('seeds', config.get('seed')))
    runs = []
    for model in models:
        for noise in noises:
            for seed in seeds:
                runs.append({
                    'tag': f'run{len(runs):02d}',
                    'model_name': model,
                    'backbone_noise': noise,
                    'seed': seed,
                })
    return runs


def mpnn_command(inputs, config, outdir, num_seqs=8, temperature='0.15 0.2',
                 model_name=None, backbone_noise=None, seed=None):
    prog = Path(config.get('proteinmpnn_dir', 'proteinmpnn'))
    python = config.get('proteinmpnn_python') or sys.executable
    cmd = [
        python, str(prog / 'protein_mpnn_run.py'),
        '--pdb_path', str(inputs['pdb']),
        '--fixed_positions_jsonl', str(inputs['fixed']),
        '--tied_positions_jsonl', str(inputs['ties']),
        '--chain_id_jsonl', str(inputs['design_chains']),
        '--num_seq_per_target', str(num_seqs),
        '--sampling_temp', temperature,
        '--out_folder', str(outdir),
        '--omit_AAs', str(config.get('omit_AAs', 'C')),
        '--suppress_print', '0',
    ]
    if config.get('soluble_model'):
        cmd.append('--use_soluble_model')
    designed = config.get('mpnn_designed_chains')
    if designed:
        cmd += ['--pdb_path_chains', ' '.join(str(c) for c in designed)]
    model = model_name if model_name is not None else config.get('model_name')
    if model:
        cmd += ['--model_name', model]
    noise = backbone_noise if backbone_noise is not None \
        else config.get('backbone_noise')
    if noise is not None:
        cmd += ['--backbone_noise', str(noise)]
    if seed is not None:
        cmd += ['--seed', str(seed)]
    bias = config.get('bias_aa')
    if bias:
        bias_json = Path(outdir) / 'bias.jsonl'
        bias_json.write_text(json.dumps({inputs['name']: bias}) + '\n')
        cmd += ['--bias_AA_jsonl', str(bias_json)]
    return cmd


def parse_design_fasta(fa_path, chain_ids):
    records = []
    for line in Path(fa_path).read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            records.append([line[1:], []])
        else:
            records[-1][1].append(line)
    # The leading non-sampled record carries the job/PDB name; sampled records
    # begin with "T=...,sample=...". Use the job name for monomer designs whose
    # sampled headers have no name prefix.
    job_name = ''
    if records:
        head = records[0][0].split(',')[0].strip()
        if not head.startswith('T='):
            job_name = head
    designs = []
    for header, chunks in records:
        seq = ''.join(chunks)
        attr = {}
        for token in header.split(','):
            pair = token.strip().split('=')
            if len(pair) == 2:
                attr[pair[0].strip()] = pair[1].strip()
        # ProteinMPNN writes a leading non-sampled record (input sequence,
        # no T=/sample= attrs); skip it so it is not treated as a design.
        if 'T' not in attr and 'sample' not in attr:
            continue
        if '/' in seq:
            chains = seq.split('/')
        elif len(chain_ids) == 1:
            chains = [seq]
        else:
            continue
        if len(chains) != len(chain_ids):
            continue
        name = header.split(',')[0].strip()
        if not name or name.startswith('T='):
            name = job_name
        designs.append({
            'name': name,
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