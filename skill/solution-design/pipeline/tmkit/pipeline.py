import argparse
import importlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

from .pdbio import load_target, harmonize_chains, reference_numbering, ca_coords, \
    residue_sasa_from_file, full_length_sequence
from .mask import build_mask, summary
from .design import prepare_inputs, mpnn_command, mpnn_runs, \
    parse_design_fasta, assemble_designs
from .predict import multimer_fasta, write_external_inputs, predict_command, \
    find_rank1_models, evaluate_models, _aa1
from .report import passes_filters, write_tagged_pdb, write_report
from .stability import sequence_solubility, sequence_properties, \
    designable_positions


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
    n_digits = {c: len(str(len(ref[c]))) for c in chains}
    report_tm = []
    for seg in tm or []:
        start, end = _resolve_seg_bounds(seg, ref, chains)
        report_tm.append((start, end))
    payload = {
        'name': cfg['name'], 'pdb_chain': chains,
        'ref_seqlen': {c: len(ref[c]) for c in chains},
        'gaps': {}, 'tm_segments': [list(s) for s in report_tm],
        'alignment': '{' + ','.join(
            f'&CGO[{n_digits[c]}]{c}' for c in chains) + '}',
        'note': 'CGO positional numberings; n_digits:chain in chain_match; '
                'tm order 1=transmembrane, 2=combine, 3=exclude all.',
    }
    for c in chains:
        n_mapped = len(reference_numbering(target, c, ref[c]))
        payload['gaps'][c] = len(ref[c]) - n_mapped
        print(f'chain {c}: {len(list(target.chains[c].residues))} residues modeled, '
              f'{payload["gaps"][c]} UCD gaps, fragment rule: '
              f'{cfg.get("assign_map", "")}')
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
    mask = build_mask(target, chains, cfg)
    print(summary(mask))
    tagged = outdir / 'outputs' / f'{cfg["name"]}_designable.pdb'
    tagged.parent.mkdir(parents=True, exist_ok=True)
    write_tagged_pdb(target, chains, mask, tagged)
    print(f'tagged PDB (B=1 designable) -> {tagged}')
    return mask


def step_design(target, cfg, outdir, num_seqs, temperature, dry_run):
    chains, ref, tm, ecd, ecd_lookup = resolve_config(cfg)
    mask = build_mask(target, chains, cfg)
    des_dir = outdir / 'outputs' / 'design_inputs'
    inputs = prepare_inputs(target, chains, mask, cfg, des_dir,
                            name=cfg['name'])
    out = outdir / 'outputs' / 'designs'
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
        print(f'ProteinMPNN command ({run["tag"]} '
              f'model={run["model_name"]} noise={run["backbone_noise"]} '
              f'seed={run["seed"]}):')
        print('  ' + ' \\\n  '.join(cmd))
        if not execute:
            continue
        subprocess.run(cmd, check=False)
        fa = run_out / 'seqs' / f'{inputs["name"]}.fa'
        if not fa.exists():
            print(f'WARNING: no design FASTA at {fa}', file=sys.stderr)
            continue
        parsed = parse_design_fasta(fa, chains)
        print(f'parsed {len(parsed)} designs from {fa}')
        designs.extend(parsed)
    return assemble_designs(designs, cfg, target, chains, out, name=cfg['name'])


def step_native(target, cfg, outdir):
    chains, ref, tm, ecd, ecd_lookup = resolve_config(cfg)
    from .metrics import tm_score, surface_hydrophobic_fraction, \
        surface_hydrophobic_fraction_at
    from .pdbio import write_structure
    mask = build_mask(target, chains, cfg)
    tmp = Path(outdir) / 'outputs' / '_native_target.pdb'
    tmp.parent.mkdir(parents=True, exist_ok=True)
    write_structure(target, tmp, chains=chains)
    all_sasa = residue_sasa_from_file(tmp)
    aa_map = {}
    with open(tmp) as f:
        for line in f:
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                aa_map[(line[21], int(line[22:26]))] = _aa1(line[17:20].strip())
    tm_set = set(r for a, b in (tm or []) for r in range(int(a), int(b) + 1))
    des_pos, tm_pos = [], []
    for cid in chains:
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
        'plddt': float('nan'), 'ecd_rmsd': 0.0,
        'surf_hydro': surface_hydrophobic_fraction(all_sasa, aa_map),
        'surf_hydro_des': surface_hydrophobic_fraction_at(
            all_sasa, aa_map, des_pos),
        'surf_hydro_tm': surface_hydrophobic_fraction_at(
            all_sasa, aa_map, tm_pos),
        'iptm': float('nan'), 'ptm': float('nan'),
        'chain_symmetry_rmsd': float('nan'), 'interface_recovery': float('nan'),
        'seq_recovery': 1.0, 'score': float('nan'), 'mutations': 0,
        'chains_full': full, 'model': str(target.path),
    }
    return [row]


def step_stability(rows, target, cfg):
    chains = cfg['chains']
    mask = build_mask(target, chains, cfg)
    ref = {c: cfg['reference_sequences'][c] for c in chains}
    tm = cfg.get('tm_segments')
    for r in rows:
        scores, resid, overp = [], [], []
        props = {k: [] for k in ('net_charge', 'abs_net_charge',
                                 'isoelectric_point', 'hydrophobic_moment',
                                 'max_hydro_window')}
        interf = 0
        for c in chains:
            des_pos = designable_positions(mask, c)
            seq = r['chains_full'][c]
            stats = sequence_solubility(seq, ref[c], des_pos)
            scores.append(stats['soluble_score'])
            resid.append(stats['residual_hydro'])
            overp.append(stats['over_polar'])
            extra = sequence_properties(seq, des_pos, tm)
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
            if k == 'max_hydro_window':
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

    tool = cfg.get('predict_tool', 'colabfold')
    print('[structure prediction]')
    if tool == 'colabfold':
        n_fail += ok('colabfold_batch on PATH') \
            if shutil.which('colabfold_batch') \
            else fail('colabfold_batch not found on PATH (conda: conda install -c conda-forge colabfold)')
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
    elif tool in ('boltz', 'af3', 'chai', 'external', 'results', 'none'):
        print(f'  [INFO] predict_tool={tool}; run the predictor externally '
              f'and pass --results-dir (PDB or mmCIF models)')

    print()
    print(f'result: {"ALL OK" if n_fail == 0 else f"{n_fail} FAILURE(S)"}')
    return n_fail


def run_target(cfg, pdb_path, args, outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    cfg = dict(cfg)
    if args.chain_ids:
        cfg['chains'] = [c.strip() for c in args.chain_ids.split(',')]
    if args.target_name:
        cfg['name'] = args.target_name
    cfg.setdefault('name', Path(pdb_path).stem)
    if args.results_dir:
        cfg['results_dir'] = args.results_dir
    filters = cfg.get('filters') or {}
    if len(cfg['chains']) <= 1:
        for key in ('iptm', 'chain_symmetry_rmsd', 'interface_mut',
                   'interface_recovery'):
            filters.pop(key, None)
    num_seqs = args.num_seqs if args.num_seqs is not None \
        else int(cfg.get('num_seqs', 8))
    temperature = args.temperature or str(cfg.get('temperature', '0.15 0.2'))

    target, ecd_lookup = load_target(pdb_path, cfg.get('ecd_segments') or [])
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
                subprocess.run(cmd, check=False)
                model_paths = find_rank1_models(
                    outdir / 'outputs' / 'predict', names, patterns=patterns)
            else:
                model_paths = {}
            missing = [n for n in names if n not in model_paths]
            if missing:
                print(f'WARNING: no rank-1 model found for '
                      f'{len(missing)}/{len(names)} designs '
                      f'(e.g. {missing[0]})', file=sys.stderr)
            evaluated = evaluate_models(model_paths, predict_rows, chains,
                                        target, mask, ecd_lookup, ref,
                                        outdir / 'outputs' / 'predict',
                                        pocket_residues=cfg.get(
                                            'pocket_residues'),
                                         tm_segments=cfg.get('tm_segments'),
                                         pocket_shell_dist=float(cfg.get(
                                             'pocket_shell_dist', 8.0)),
                                         interface_recovery_dist=float(cfg.get(
                                             'interface_recovery_dist', 8.0)))
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