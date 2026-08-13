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
from .design import prepare_inputs, mpnn_command, parse_design_fasta, assemble_designs
from .predict import multimer_fasta, predict_command, find_rank1_models, \
    evaluate_models, _aa1
from .report import passes_filters, write_tagged_pdb, write_report
from .stability import sequence_solubility, designable_positions


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
    cmd = mpnn_command(inputs, cfg, out, num_seqs=num_seqs, temperature=temperature)
    print('ProteinMPNN command:')
    print('  ' + ' \\\n  '.join(cmd))
    run = not dry_run
    if run:
        in_dir = Path(cfg.get('proteinmpnn_dir', 'proteinmpnn'))
        run = (in_dir / 'protein_mpnn_run.py').exists()
        if not run:
            print('WARNING: protein_mpnn_run.py not found; dry-run only.',
                  file=sys.stderr)
    if run:
        subprocess.run(cmd, check=False)
        fa = Path(out) / 'seqs' / f'{inputs["name"]}.fa'
        designs = parse_design_fasta(fa, chains)
        print(f'parsed {len(designs)} designs from {fa}')
    else:
        designs = []
    return assemble_designs(designs, cfg, target, chains, out, name=cfg['name'])


def step_native(target, cfg, outdir):
    chains, ref, tm, ecd, ecd_lookup = resolve_config(cfg)
    from .metrics import tm_score, surface_hydrophobic_fraction
    all_sasa = residue_sasa_from_file(target.path)
    aa_map = {}
    with open(target.path) as f:
        for line in f:
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                aa_map[(line[21], int(line[22:26]))] = _aa1(line[17:20].strip())
    tm_total = []
    for cid in chains:
        A, _ = ca_coords(target, cid)
        if len(A):
            tm_total.append(tm_score(A, A))
    full = {}
    for cid in chains:
        ref_seq = ref[cid]
        des_pos = [i + 1 for i, aa in enumerate(ref_seq)]
        full[cid] = full_length_sequence(ref_seq, des_pos,
                                         [aa for aa in ref_seq])
    row = {
        'name': 'native', 'pass': None,
        'tm_score': float(np.mean(tm_total)) if tm_total else float('nan'),
        'plddt': float('nan'), 'ecd_rmsd': 0.0,
        'surf_hydro': surface_hydrophobic_fraction(all_sasa, aa_map),
        'seq_recovery': 1.0, 'score': float('nan'), 'mutations': 0,
        'chains_full': full, 'model': str(target.path),
    }
    return [row]


def step_stability(rows, target, cfg):
    chains = cfg['chains']
    mask = build_mask(target, chains, cfg)
    ref = {c: cfg['reference_sequences'][c] for c in chains}
    for r in rows:
        scores, resid, overp = [], [], []
        interf = 0
        for c in chains:
            des_pos = designable_positions(mask, c)
            stats = sequence_solubility(r['chains_full'][c], ref[c], des_pos)
            scores.append(stats['soluble_score'])
            resid.append(stats['residual_hydro'])
            overp.append(stats['over_polar'])
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
    elif tool == 'none':
        print('  [INFO] predict_tool=none; structure prediction disabled')

    print()
    print(f'result: {"ALL OK" if n_fail == 0 else f"{n_fail} FAILURE(S)"}')
    return n_fail


def main(argv=None):
    ap = argparse.ArgumentParser(description='Soluble TM redesign pipeline')
    ap.add_argument('--config', default=None, help='yaml config (optional for --check-deps)')
    ap.add_argument('--pdb', default=None)
    ap.add_argument('--chain-ids', default=None, help='comma list, override config')
    ap.add_argument('--steps', default='annotate,mask,native,predict,stability,report',
                    help='comma list: annotate,mask,design,native,predict,stability,report')
    ap.add_argument('--dry-run', action='store_true',
                    help='print external commands, do not run')
    ap.add_argument('--out', default='.')
    ap.add_argument('--target-name', default=None)
    ap.add_argument('--num-seqs', type=int, default=None,
                    help='designs per target (default: config num_seqs)')
    ap.add_argument('--temperature', default=None,
                    help='MPNN sampling temperatures (default: config)')
    ap.add_argument('--check-deps', action='store_true',
                    help='check libraries/tools/config paths and exit')
    args = ap.parse_args(argv)

    if args.check_deps:
        cfg = None
        if args.config:
            import yaml
            cfg = yaml.safe_load(Path(args.config).read_text())
        return check_dependencies(cfg, args.pdb)

    if not args.config or not args.pdb:
        ap.error('--config and --pdb are required (use --check-deps for a lightweight check)')

    import yaml
    cfg = yaml.safe_load(Path(args.config).read_text())
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    if args.chain_ids:
        cfg['chains'] = [c.strip() for c in args.chain_ids.split(',')]
    if args.target_name:
        cfg['name'] = args.target_name
    cfg.setdefault('name', Path(args.pdb).stem)
    num_seqs = args.num_seqs if args.num_seqs is not None \
        else int(cfg.get('num_seqs', 8))
    temperature = args.temperature or str(cfg.get('temperature', '0.15 0.2'))

    target, ecd_lookup = load_target(args.pdb, cfg.get('ecd_segments') or [])
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
            fasta = multimer_fasta(predict_rows, chains,
                                   outdir / 'outputs' / 'predict_in')
            cmd = predict_command(cfg, fasta, outdir / 'outputs' / 'predict')
            print('Prediction command:')
            print('  ' + ' \\\n  '.join(cmd or []))
            names = [r['name'].replace(':', '_') for r in predict_rows]
            if args.dry_run:
                model_paths = {}
            elif cfg.get('results_dir'):
                model_paths = find_rank1_models(cfg['results_dir'], names)
            elif cmd:
                subprocess.run(cmd, check=False)
                model_paths = find_rank1_models(
                    outdir / 'outputs' / 'predict', names)
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
                                            'pocket_residues'))
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


if __name__ == '__main__':
    sys.exit(main())