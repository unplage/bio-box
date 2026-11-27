import argparse
import importlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

from .pdbio import load_target, harmonize_chains, \
    reference_numbering, ca_coords, residue_sasa_from_file, \
    full_length_sequence
from .config import normalize_config, resolve_template_path, validate_config
from .topology import resolve_topology, resolve_topologies
from .mask import build_mask, summary
from .design import prepare_inputs, prepare_backbones, mpnn_command, \
    mpnn_runs, parse_design_fasta, assemble_designs
from .predict import multimer_fasta, write_external_inputs, predict_command, \
    find_rank1_models, evaluate_models, _aa1
from .report import passes_filters, write_tagged_pdb, write_report
from .constraints import symmetry_groups
from .stability import sequence_solubility, sequence_properties, \
    designable_positions


_CLI_CONFIG_KEYS = (
    'proteinmpnn_dir', 'proteinmpnn_python', 'design_mode',
    'rfdiffusion_script',
    'rfdiffusion_dir',
    'predict_tool', 'predict_executable', 'colabfold_executable',
    'boltz_executable', 'boltz_devices', 'boltz_use_msa_server',
    'boltz_no_kernels', 'require_predictions', 'strict_predict',
    'alphafold_script',
)


def apply_cli_overrides(cfg, args):
    """Overlay tool locations supplied on the CLI (never persisted to YAML)."""
    out = dict(cfg or {})
    for key in _CLI_CONFIG_KEYS:
        value = getattr(args, key, None)
        if value is not None:
            out[key] = value
    return out


def resolve_config(cfg):
    chains = cfg['chains']
    ref = {c: cfg['reference_sequences'][c] for c in chains}
    tm = cfg.get('tm_segments')
    ecd = cfg.get('ecd_segments', [])
    ecd_lookup = [tuple(int(x) for x in s) for s in ecd]
    return chains, ref, tm, ecd, ecd_lookup


def step_annotate(target, cfg):
    chains, ref, tm, ecd, ecd_lookup = resolve_config(cfg)
    harmonize_chains(target, chains)
    topologies = resolve_topologies(cfg, {c: len(ref[c]) for c in chains},
                                    chains)
    n_digits = {c: len(str(len(ref[c]))) for c in chains}
    payload = {
        'name': cfg['name'], 'pdb_chain': chains,
        'membrane_class': next(iter(topologies.values()))['class']
        if topologies else cfg.get('membrane_class'),
        'topology_by_chain': {},
        'ref_seqlen': {c: len(ref[c]) for c in chains},
        'gaps': {},
        'alignment': '{' + ','.join(
            f'&CGO[{n_digits[c]}]{c}' for c in chains) + '}',
        'note': 'CGO positional numberings; n_digits:chain in chain_match; '
                'tm order 1=transmembrane, 2=combine, 3=exclude all.',
    }
    for c in chains:
        topology = topologies[c]
        n_mapped = len(reference_numbering(target, c, ref[c]))
        payload['gaps'][c] = len(ref[c]) - n_mapped
        payload['topology_by_chain'][c] = {
            'nterm': topology['nterm'], 'cterm': topology['cterm'],
            'tm_segments': [list(s) for s in topology['tm_segments']],
            'ecd_segments': [list(s) for s in topology['ecd_segments']],
            'icd_segments': [list(s) for s in topology['icd_segments']],
            'tm_type': topology['tm_type'],
            'warnings': list(topology.get('warnings') or []),
        }
        print(f'chain {c}: {len(list(target.chains[c].residues))} residues modeled, '
              f'{payload["gaps"][c]} UCD gaps, fragment rule: '
              f'{cfg.get("assign_map", "")}')
        print(f'  topology: class={topology["class"]} '
              f'nterm={topology["nterm"]} cterm={topology["cterm"]} '
              f'TM={len(topology["tm_segments"])} '
              f'ECD={len(topology["ecd_segments"])} '
              f'ICD={len(topology["icd_segments"])}')
        for warning in topology.get('warnings') or []:
            print(f'  [topology warning] {warning}', file=sys.stderr)
    print(json.dumps(payload, indent=2))


def _resolve_seg_bounds(seg, ref, chains):
    if isinstance(seg, list) and len(seg) == 2:
        return int(seg[0]), int(seg[1])
    if isinstance(seg, dict):
        for c in chains:
            s = seg.get(c)
            if s is not None:
                return int(s[0]), int(s[1])
    raise ValueError(f'cannot resolve segment {seg!r}')


def step_mask(target, cfg, outdir):
    chains, ref, tm, ecd, ecd_lookup = resolve_config(cfg)
    topologies = resolve_topologies(cfg, {c: len(ref[c]) for c in chains},
                                    chains)
    topology = topologies.get(chains[0]) if chains else None
    for warning in (topology or {}).get('warnings') or []:
        print(f'[topology warning] {warning}', file=sys.stderr)
    mask = build_mask(target, chains, cfg)
    print(summary(mask))
    tagged = outdir / 'outputs' / f'{cfg["name"]}_designable.pdb'
    tagged.parent.mkdir(parents=True, exist_ok=True)
    write_tagged_pdb(target, chains, mask, tagged)
    print(f'tagged PDB (B=1 designable) -> {tagged}')
    return mask


def _design_on_backbone(target, cfg, outdir, num_seqs, temperature, dry_run,
                        tag):
    chains, ref, tm, ecd, ecd_lookup = resolve_config(cfg)
    mask = build_mask(target, chains, cfg)
    des_dir = outdir / 'outputs' / 'design_inputs' / tag
    inputs = prepare_inputs(target, chains, mask, cfg, des_dir,
                            name=f'{cfg["name"]}_{tag}')
    out = outdir / 'outputs' / 'designs' / tag
    out.mkdir(parents=True, exist_ok=True)
    runs = mpnn_runs(cfg)
    execute = not dry_run
    if execute:
        in_dir = Path(cfg.get('proteinmpnn_dir', 'proteinmpnn'))
        execute = (in_dir / 'protein_mpnn_run.py').exists()
        if not execute:
            print('WARNING: protein_mpnn_run.py not found; dry-run only.',
                  file=sys.stderr)
    designs = []
    for run in runs:
        run_out = out / run['tag']
        run_out.mkdir(parents=True, exist_ok=True)
        cmd = mpnn_command(inputs, cfg, run_out, num_seqs=num_seqs,
                           temperature=temperature,
                           model_name=run['model_name'],
                           backbone_noise=run['backbone_noise'],
                           seed=run['seed'])
        print(f'ProteinMPNN command ({tag}/{run["tag"]} '
              f'model={run["model_name"]} noise={run["backbone_noise"]} '
              f'seed={run["seed"]}):')
        print('  ' + ' \\\n  '.join(cmd))
        if not execute:
            continue
        proc = subprocess.run(cmd)
        if proc.returncode != 0:
            raise RuntimeError(
                f'ProteinMPNN failed with exit code {proc.returncode} '
                f'({tag}/{run["tag"]})')
        fa = run_out / 'seqs' / f'{inputs["name"]}.fa'
        if not fa.exists():
            print(f'WARNING: no design FASTA at {fa}', file=sys.stderr)
            continue
        parsed = parse_design_fasta(fa, chains)
        print(f'parsed {len(parsed)} designs from {fa}')
        for d in parsed:
            d['name'] = f'{tag}_{d["name"]}'
        designs.extend(parsed)
    return assemble_designs(designs, cfg, target, chains, out,
                            name=cfg['name'])


def step_design(target, cfg, outdir, num_seqs, temperature, dry_run):
    chains, ref, tm, ecd, ecd_lookup = resolve_config(cfg)
    backbones = prepare_backbones(cfg, outdir, dry_run=dry_run)
    if not backbones:
        backbones = [None]
    all_designs = []
    topologies = resolve_topologies(cfg, {c: len(ref[c]) for c in chains},
                                    chains)
    topology = topologies.get(chains[0]) if chains else None
    for i, backbone in enumerate(backbones):
        if backbone is None:
            t, tag = target, 'template'
        else:
            t, _ = load_target(str(backbone), topology['ecd_segments'])
            harmonize_chains(t, chains)
            for c in chains:
                reference_numbering(t, c, ref[c])
            tag = f'b{i:02d}_{Path(backbone).stem}'
        all_designs.extend(_design_on_backbone(
            t, cfg, outdir, num_seqs, temperature, dry_run, tag))
    return all_designs


def step_native(target, cfg, outdir):
    chains, ref, tm, ecd, ecd_lookup = resolve_config(cfg)
    from .metrics import tm_score, surface_hydrophobic_fraction, \
        surface_hydrophobic_fraction_at
    from .pdbio import write_structure
    mask = build_mask(target, chains, cfg)
    topologies = resolve_topologies(cfg, {c: len(ref[c]) for c in chains},
                                    chains)
    tmp = Path(outdir) / 'outputs' / '_native_target.pdb'
    tmp.parent.mkdir(parents=True, exist_ok=True)
    write_structure(target, tmp, chains=chains)
    all_sasa = residue_sasa_from_file(tmp)
    aa_map = {}
    with open(tmp) as f:
        for line in f:
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                aa_map[(line[21], int(line[22:26]))] = _aa1(line[17:20].strip())
    des_pos, tm_pos = [], []
    for cid in chains:
        tm_set = set(r for a, b in topologies[cid]['tm_segments']
                     for r in range(int(a), int(b) + 1))
        mapping = target.reference_map.get(cid, {})
        for residue in target.chains[cid].residues:
            resnum = residue.id[1]
            refn = mapping.get(resnum, resnum)
            if mask.get(cid, {}).get(refn, {}).get('designable'):
                des_pos.append((cid, resnum))
            if refn in tm_set:
                tm_pos.append((cid, resnum))
    tm_total = []
    for cid in chains:
        A, _ = ca_coords(target, cid)
        if len(A):
            tm_total.append(tm_score(A, A))
    full = {}
    for cid in chains:
        ref_seq = ref[cid]
        des_positions = [i + 1 for i, aa in enumerate(ref_seq)]
        full[cid] = full_length_sequence(ref_seq, des_positions,
                                         [aa for aa in ref_seq])
    row = {
        'name': 'native', 'pass': None,
        'tm_score': float(np.mean(tm_total)) if tm_total else float('nan'),
        'complex_tm_score': float(np.mean(tm_total)) if tm_total else float('nan'),
        'plddt': float('nan'), 'plddt_ecd': float('nan'),
        'plddt_tm': float('nan'), 'plddt_icd': float('nan'),
        'plddt_loop': float('nan'),
        'ecd_rmsd': 0.0, 'ecd_internal_rmsd': 0.0,
        'ecd_centroid_shift': 0.0, 'ecd_tm_angle': 0.0,
        'ligand_rmsd': float('nan'),
        'surf_hydro': surface_hydrophobic_fraction(all_sasa, aa_map),
        'surf_hydro_des': surface_hydrophobic_fraction_at(
            all_sasa, aa_map, des_pos),
        'surf_hydro_tm': surface_hydrophobic_fraction_at(
            all_sasa, aa_map, tm_pos),
        'iptm': float('nan'), 'ptm': float('nan'),
        'chain_symmetry_rmsd': float('nan'), 'interface_recovery': float('nan'),
        'seq_recovery': 1.0, 'score': float('nan'), 'mutations': 0,
        'chains_full': full, 'model': str(target.path),
        'predicted': True, 'status': 'native_control',
    }
    return [row]


def step_stability(rows, target, cfg):
    chains = cfg['chains']
    mask = build_mask(target, chains, cfg)
    ref = {c: cfg['reference_sequences'][c] for c in chains}
    topologies = resolve_topologies(cfg, {c: len(ref[c]) for c in chains},
                                    chains)
    for r in rows:
        scores, resid, overp = [], [], []
        props = {k: [] for k in ('net_charge', 'abs_net_charge',
                                 'isoelectric_point', 'hydrophobic_moment',
                                 'max_hydro_window', 'max_charge_run',
                                 'max_polar_window')}
        interf = 0
        for c in chains:
            des_pos = designable_positions(mask, c)
            seq = r['chains_full'][c]
            stats = sequence_solubility(seq, ref[c], des_pos)
            scores.append(stats['soluble_score'])
            resid.append(stats['residual_hydro'])
            overp.append(stats['over_polar'])
            topo = topologies.get(c) or {}
            extra = sequence_properties(seq, des_pos,
                                        topo.get('tm_segments'),
                                        tm_type=topo.get('tm_type', 'helix'))
            for k in props:
                props[k].append(extra[k])
            for refn, info in mask[c].items():
                if info['interface'] and refn - 1 < len(r['chains_full'][c]) \
                        and refn - 1 < len(ref[c]):
                    if r['chains_full'][c][refn - 1] != ref[c][refn - 1]:
                        interf += 1
        r['soluble_score'] = float(np.mean([s for s in scores if s == s])) \
            if any(s == s for s in scores) else float('nan')
        r['residual_hydro'] = float(np.mean([s for s in resid if s == s])) \
            if any(s == s for s in resid) else float('nan')
        r['over_polar'] = float(np.mean([s for s in overp if s == s])) \
            if any(s == s for s in overp) else float('nan')
        for k, vals in props.items():
            vals = [v for v in vals if v == v]
            if k in ('max_hydro_window', 'max_charge_run',
                     'max_polar_window'):
                r[k] = float(np.max(vals)) if vals else float('nan')
            else:
                r[k] = float(np.mean(vals)) if vals else float('nan')
        r['interface_mut'] = interf
    return rows


def check_dependencies(cfg=None, pdb_path=None):
    """Report availability of python libs, external tools and configured paths."""
    def fail(msg):
        print(f'  [FAIL] {msg}')
        return 1
    def ok(msg):
        print(f'  [OK]   {msg}')
        return 0

    n_fail = 0
    print('[python libraries]')
    for module, note in [('numpy', 'metrics'), ('Bio', 'PDB parsing'),
                         ('yaml', 'config'), ('freesasa', 'SASA')]:
        try:
            importlib.import_module(module)
            n_fail += ok(f'{module} ({note})')
        except ImportError:
            n_fail += fail(f'{module} ({note}) -- pip install -r requirements.txt')
    for module, note in [('scipy', 'optional KD-tree')]:
        try:
            importlib.import_module(module)
            n_fail += ok(f'{module} ({note})')
        except ImportError:
            print('  [WARN] scipy missing (optional; falls back to Biopython KD-tree)')

    if pdb_path:
        print('[input]')
        n_fail += ok(str(pdb_path)) if Path(pdb_path).exists() \
            else fail(f'input PDB missing: {pdb_path}')

    if cfg is None:
        print()
        print(f'result: {"ALL OK" if n_fail == 0 else f"{n_fail} FAILURE(S)"}')
        return n_fail

    print('[config paths]')
    if cfg.get('proteinmpnn_python'):
        py = Path(cfg['proteinmpnn_python'])
        n_fail += ok(f'ProteinMPNN python: {py}') if py.exists() \
            else fail(f'proteinmpnn_python missing: {py}')
    mpnn = cfg.get('proteinmpnn_dir')
    if mpnn:
        mpnn_dir = Path(mpnn)
        if (mpnn_dir / 'protein_mpnn_run.py').exists():
            n_fail += ok(f'ProteinMPNN: {mpnn_dir / "protein_mpnn_run.py"}')
        else:
            n_fail += fail(f'protein_mpnn_run.py not found under proteinmpnn_dir={mpnn}')
        if cfg.get('soluble_model'):
            if (mpnn_dir / 'soluble_model_weights').exists():
                n_fail += ok('soluble_model_weights present')
            else:
                n_fail += fail('soluble_model requires soluble_model_weights under proteinmpnn_dir')
    if cfg.get('results_dir'):
        rd = Path(cfg['results_dir'])
        n_fail += ok(f'results_dir: {rd}') if rd.exists() \
            else fail(f'results_dir set but missing: {rd}')
    if cfg.get('backbone_dir'):
        bd = Path(cfg['backbone_dir'])
        n_fail += ok(f'backbone_dir: {bd}') if bd.exists() \
            else fail(f'backbone_dir set but missing: {bd}')
    mode = str(cfg.get('design_mode', 'sequence')).lower()
    if mode in ('scaffold', 'rf_diffusion', 'rfdiffusion',
                'partial_diffusion'):
        script = cfg.get('rfdiffusion_script')
        if not script:
            script = Path(cfg.get('rfdiffusion_dir', 'RFdiffusion')) / \
                'scripts' / 'run_inference.py'
        n_fail += ok(f'RFdiffusion: {script}') if Path(script).exists() \
            else fail(f'RFdiffusion run_inference.py not found: {script}')
    if cfg.get('predict_executable'):
        n_fail += ok(f'predict_executable: {cfg["predict_executable"]}') \
            if Path(cfg['predict_executable']).exists() \
            else fail(f'predict_executable missing: {cfg["predict_executable"]}')
    if cfg.get('chains') and cfg.get('reference_sequences'):
        try:
            ref_lengths = {c: len(cfg['reference_sequences'][c])
                           for c in cfg['chains']}
            for warning in validate_config(cfg, ref_lengths):
                print(f'  [WARN] {warning}')
        except Exception as exc:
            n_fail += fail(f'config validation: {exc}')

    tool = cfg.get('predict_tool', 'colabfold')
    print('[structure prediction]')
    if tool == 'colabfold':
        exe = cfg.get('colabfold_executable') or \
            shutil.which('colabfold_batch')
        n_fail += ok(f'colabfold executable: {exe}') if exe \
            and Path(exe).exists() \
            else fail('colabfold_batch not found; pass '
                      '--colabfold-executable /path/to/colabfold_batch')
        model = cfg.get('colabfold_model', 'alpha2_multimer_v3')
        if model not in ('alpha2_multimer_v3', 'alphafold2_multimer_v3', 'alpha2_multimer_v2', 'alphafold2_multimer_v2', 'alphafold2_multimer_v1'):
            print(f'  [WARN] unusual colabfold_model: {model}')
    elif tool == 'alphafold':
        script = cfg.get('alphafold_script')
        n_fail += ok(f'alphafold_script: {script}') \
            if script and Path(script).exists() \
            else fail('alphafold_script missing (set absolute path to run_alphafold.py)')
        for key, flag in [('data_dir', '--data_dir'), ('uniref90', '--uniref90_database_path'),
                          ('mgnify', '--mgnify_database_path'), ('bfd', '--bfd_database_path'),
                          ('uniref30', '--uniref30_database_path'),
                          ('pdb_seqres', '--pdb_seqres_database_path'),
                          ('template_mmcif', '--template_mmcif_dir'),
                          ('obsolete', '--obsolete_pdbs_path')]:
            if cfg.get(key):
                p = Path(cfg[key])
                n_fail += ok(f'{flag}: {p}') if p.exists() \
                    else fail(f'{flag} set but missing: {p}')
    elif tool in ('boltz', 'boltz2'):
        exe = cfg.get('boltz_executable') or shutil.which('boltz')
        n_fail += ok(f'boltz executable: {exe}') if exe \
            and Path(exe).exists() \
            else fail('boltz not found; pass --boltz-executable /path/to/boltz')
    elif tool in ('af3', 'chai', 'external', 'results', 'none'):
        print(f'  [INFO] predict_tool={tool}; run the predictor externally '
              f'and pass --results-dir (PDB or mmCIF models)')

    print()
    print(f'result: {"ALL OK" if n_fail == 0 else f"{n_fail} FAILURE(S)"}')
    return n_fail


def run_target(cfg, pdb_path, args, outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    cfg = apply_cli_overrides(cfg, args)
    if args.chain_ids:
        cfg['chains'] = [c.strip() for c in args.chain_ids.split(',')]
    if args.target_name:
        cfg['name'] = args.target_name
    cfg.setdefault('name', Path(pdb_path).stem)
    if args.results_dir:
        cfg['results_dir'] = args.results_dir
    cfg = normalize_config(cfg)
    ref_lengths = {c: len(cfg['reference_sequences'][c])
                   for c in cfg['chains']}
    warnings = validate_config(cfg, ref_lengths)
    for warning in warnings:
        print(f'[config warning] {warning}', file=sys.stderr)
    topologies = resolve_topologies(cfg, ref_lengths, cfg['chains'])
    pdb_path = resolve_template_path(cfg, pdb_path)
    filters = cfg.get('filters') or {}
    if len(cfg['chains']) <= 1:
        for key in ('iptm', 'chain_symmetry_rmsd', 'interface_mut',
                   'interface_recovery'):
            filters.pop(key, None)
    else:
        groups = symmetry_groups(cfg, cfg['chains'])
        if not groups:
            filters.pop('chain_symmetry_rmsd', None)
    if not any(topologies[c].get('ecd_segments') for c in cfg['chains']):
        filters.pop('ecd_rmsd', None)
        filters.pop('plddt_ecd', None)
    num_seqs = args.num_seqs if args.num_seqs is not None \
        else int(cfg.get('num_seqs', 8))
    temperature = args.temperature or str(cfg.get('temperature', '0.15 0.2'))

    ecd_segments = []
    for topology in topologies.values():
        for segment in topology.get('ecd_segments') or []:
            if tuple(segment) not in {tuple(s) for s in ecd_segments}:
                ecd_segments.append(segment)
    target, ecd_lookup = load_target(pdb_path, ecd_segments)
    chains, ref, tm, ecd, _ec = resolve_config(cfg)
    harmonize_chains(target, chains)
    for c in chains:
        reference_numbering(target, c, ref[c])
    mask = None
    rows = []
    for step in [s.strip() for s in args.steps.split(',') if s.strip()]:
        if step == 'annotate':
            step_annotate(target, cfg)
        elif step == 'mask':
            mask = step_mask(target, cfg, outdir)
        elif step == 'design':
            rows = step_design(target, cfg, outdir, num_seqs,
                               temperature, args.dry_run)
        elif step == 'native':
            rows += step_native(target, cfg, outdir)
        elif step == 'predict':
            if mask is None:
                mask = build_mask(target, chains, cfg)
            if not rows:
                print('no designs to predict; run design/native first',
                      file=sys.stderr)
                continue
            predict_rows = rows
            if not cfg.get('predict_native', False):
                predict_rows = [r for r in rows if r.get('name') != 'native']
            if not predict_rows:
                print('nothing to predict (native excluded by '
                      'predict_native=false)', file=sys.stderr)
                continue
            predict_in = outdir / 'outputs' / 'predict_in'
            fasta = multimer_fasta(predict_rows, chains, predict_in)
            write_external_inputs(predict_rows, chains, predict_in)
            cmd = predict_command(cfg, fasta, outdir / 'outputs' / 'predict')
            print('Prediction command:')
            print('  ' + ' \\\n  '.join(cmd or []))
            names = [r['name'].replace(':', '_') for r in predict_rows]
            patterns = cfg.get('results_globs')
            if args.dry_run:
                model_paths = {}
            elif cfg.get('results_dir'):
                model_paths = find_rank1_models(cfg['results_dir'], names,
                                                patterns=patterns)
            elif cmd:
                proc = subprocess.run(cmd)
                if proc.returncode != 0:
                    msg = (f'prediction command failed with exit code '
                           f'{proc.returncode}: {cmd[0]}')
                    if cfg.get('strict_predict', True):
                        raise RuntimeError(msg)
                    print(f'WARNING: {msg}', file=sys.stderr)
                model_paths = find_rank1_models(
                    outdir / 'outputs' / 'predict', names, patterns=patterns)
            else:
                model_paths = {}
            missing = [n for n in names if n not in model_paths]
            if missing:
                print(f'WARNING: no rank-1 model found for '
                      f'{len(missing)}/{len(names)} designs '
                      f'(e.g. {missing[0]})', file=sys.stderr)
                if cfg.get('require_predictions', False):
                    raise RuntimeError(
                        f'missing predictions for {len(missing)} design(s); '
                        f'see {outdir / "outputs" / "predict"}')
            evaluated = evaluate_models(model_paths, predict_rows, chains,
                                        target, mask, ecd_lookup, ref,
                                        outdir / 'outputs' / 'predict',
                                        pocket_residues=cfg.get(
                                            'pocket_residues'),
                                         tm_segments=cfg.get('tm_segments'),
                                         pocket_shell_dist=float(cfg.get(
                                             'pocket_shell_dist', 8.0)),
                                         interface_recovery_dist=float(cfg.get(
                                             'interface_recovery_dist', 8.0)),
                                         topologies=topologies,
                                         symmetry_groups=symmetry_groups(
                                             cfg, chains),
                                         interface_pairs=cfg.get(
                                             'interface_pairs'),
                                         sasa_threshold=float(cfg.get(
                                             'surface_rsa_threshold', 0.25)
                                             if cfg.get('relative_sasa')
                                             else cfg.get(
                                                 'surface_sasa_threshold',
                                                 40.0)),
                                         relative_sasa=bool(cfg.get(
                                             'relative_sasa', False)),
                                         cofactors=cfg.get('cofactors') or [])
            by_name = {r['name']: r for r in evaluated}
            rows = [by_name.get(r['name'], r) for r in rows]
        elif step == 'stability':
            rows = step_stability(rows, target, cfg)
        elif step == 'report':
            for r in rows:
                passes_filters(r, cfg.get('filters', {}))
            write_report(rows, cfg, outdir / 'outputs')
            if mask is None:
                mask = build_mask(target, chains, cfg)
            write_tagged_pdb(target, chains, mask,
                             outdir / 'outputs' /
                             f'{cfg["name"]}_designable.pdb')
        else:
            raise ValueError(f'unknown step {step}')
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description='Soluble TM redesign pipeline')
    ap.add_argument('--config', default=None,
                    help='yaml config file, or a directory for batch runs')
    ap.add_argument('--config-dir', default=None,
                    help='directory of *.yaml targets (batch)')
    ap.add_argument('--pdb', default=None,
                    help='template structure (PDB or mmCIF); defaults to '
                         'config template_pdb')
    ap.add_argument('--chain-ids', default=None, help='comma list, override config')
    ap.add_argument('--steps', default='annotate,mask,native,predict,stability,report',
                    help='comma list: annotate,mask,design,native,predict,stability,report')
    ap.add_argument('--dry-run', action='store_true',
                    help='print external commands, do not run')
    ap.add_argument('--out', default='.')
    ap.add_argument('--target-name', default=None)
    ap.add_argument('--results-dir', default=None,
                    help='directory of externally predicted models '
                         '(ColabFold/AlphaFold/Boltz/AF3; PDB or mmCIF)')
    ap.add_argument('--num-seqs', type=int, default=None,
                    help='designs per target (default: config num_seqs)')
    ap.add_argument('--temperature', default=None,
                    help='MPNN sampling temperatures (default: config)')
    ap.add_argument('--check-deps', action='store_true',
                    help='check libraries/tools/config paths and exit')
    ap.add_argument('--proteinmpnn-dir', dest='proteinmpnn_dir', default=None,
                    help='override config proteinmpnn_dir (not persisted)')
    ap.add_argument('--proteinmpnn-python', dest='proteinmpnn_python',
                    default=None,
                    help='override Python used for protein_mpnn_run.py')
    ap.add_argument('--design-mode', dest='design_mode', default=None,
                    help='override config design_mode (sequence/scaffold/rf_diffusion)')
    ap.add_argument('--rfdiffusion-script', dest='rfdiffusion_script',
                    default=None,
                    help='override config RFdiffusion run_inference.py')
    ap.add_argument('--rfdiffusion-dir', dest='rfdiffusion_dir', default=None,
                    help='override config RFdiffusion root directory')
    ap.add_argument('--predict-tool', dest='predict_tool', default=None,
                    help='override config predict_tool')
    ap.add_argument('--predict-executable', dest='predict_executable',
                    default=None,
                    help='override predictor executable path')
    ap.add_argument('--colabfold-executable', dest='colabfold_executable',
                    default=None,
                    help='override colabfold_batch path')
    ap.add_argument('--boltz-executable', dest='boltz_executable',
                    default=None, help='override boltz path')
    ap.add_argument('--boltz-devices', dest='boltz_devices', default=None,
                    help='number of devices for Boltz (e.g. 2)')
    ap.add_argument('--boltz-use-msa-server', dest='boltz_use_msa_server',
                    action='store_true', default=None,
                    help='generate MSAs via the MMSeqs2 server')
    ap.add_argument('--boltz-no-kernels', dest='boltz_no_kernels',
                    action='store_true', default=None,
                    help='disable optional Boltz CUDA kernels')
    ap.add_argument('--require-predictions', dest='require_predictions',
                    action='store_true', default=None,
                    help='fail if any rank-1 model is missing')
    ap.add_argument('--alphafold-script', dest='alphafold_script',
                    default=None, help='override run_alphafold.py path')
    args = ap.parse_args(argv)

    import yaml
    batch_dir = args.config_dir or (args.config
                                    if args.config and Path(args.config).is_dir()
                                    else None)

    if args.check_deps:
        if batch_dir:
            print(f'[INFO] batch config dir: {batch_dir}')
            return 0
        cfg = None
        if args.config:
            cfg = yaml.safe_load(Path(args.config).read_text())
            cfg = apply_cli_overrides(cfg, args)
        return check_dependencies(cfg, args.pdb)

    if not args.config and not batch_dir:
        ap.error('--config (file or dir) is required '
                 '(use --check-deps for a lightweight check)')

    if batch_dir:
        base = Path(batch_dir)
        files = sorted(p for p in base.iterdir()
                       if p.suffix.lower() in ('.yaml', '.yml'))
        rc = 0
        for conf in files:
            cfg = yaml.safe_load(conf.read_text())
            pdb = args.pdb or cfg.get('template_pdb') or cfg.get('pdb')
            if not pdb:
                print(f'[skip] {conf.name}: no template_pdb/pdb configured',
                      file=sys.stderr)
                continue
            out = Path(args.out) / str(cfg.get('name', conf.stem))
            print(f'=== target {cfg.get("name", conf.stem)} ===')
            rc |= run_target(cfg, pdb, args, out)
        return rc

    cfg = yaml.safe_load(Path(args.config).read_text())
    pdb = args.pdb or cfg.get('template_pdb') or cfg.get('pdb')
    if not pdb:
        ap.error('--pdb is required (or set template_pdb in the config)')
    return run_target(cfg, pdb, args, args.out)


if __name__ == '__main__':
    sys.exit(main())