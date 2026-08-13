import json
import csv
from pathlib import Path

import numpy as np


def passes_filters(row, filters):
    checks = {}
    for key, cond in filters.items():
        val = row.get(key)
        if val is None or (isinstance(val, float) and np.isnan(val)):
            checks[key] = False
            continue
        op, thr = cond[0], float(cond[1])
        if op == '>':
            checks[key] = val > thr
        elif op == '<':
            checks[key] = val < thr
        elif op == '>=':
            checks[key] = val >= thr
        elif op == '<=':
            checks[key] = val <= thr
        elif op == '==':
            checks[key] = val == thr
        elif op == '!=':
            checks[key] = val != thr
        else:
            checks[key] = True
    row['pass'] = all(checks.values()) if checks else False
    row['fail_gates'] = [k for k, v in checks.items() if not v]
    return row['pass']


def write_tagged_pdb(target, chain_ids, mask, outpath):
    m = target.structure
    keys = [('designable', 1.0), ('fixed', 0.0)]
    lines = []
    serial = 1
    for model in m:
        for chain in model:
            if chain.id not in chain_ids:
                continue
            mapping = target.reference_map.get(chain.id, {})
            for residue in chain:
                if residue.id[0] != ' ':
                    continue
                refn = mapping.get(residue.id[1], residue.id[1])
                b = 0.0
                info = mask.get(chain.id, {}).get(refn)
                if info and info['designable']:
                    b = 1.0
                resname = residue.resname
                for atom in residue:
                    coord = atom.get_coord()
                    lines.append(
                        f'ATOM  {serial:5d} {atom.name:>4} {resname:>3} '
                        f'{chain.id:>1}{residue.id[1]:4d}    '
                        f'{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}'
                        f'{1.0:6.2f}{b:6.2f}          {atom.element:>2}')
                    serial += 1
    Path(outpath).write_text('\n'.join(lines) + 'END\n')


def write_report(rows, config, outdir, top_n=10):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    core = ['name', 'pass', 'fail_gates', 'tm_score', 'plddt', 'ecd_rmsd',
            'surf_hydro', 'surf_hydro_des', 'plddt_ecd', 'pocket_rmsd',
            'exposed_arom', 'cys_des', 'seq_recovery', 'score', 'mutations',
            'soluble_score', 'residual_hydro', 'over_polar', 'interface_mut',
            'model']
    extra = []
    for r in rows:
        for k, v in r.items():
            if k in core or k in extra or k in ('chains_full', 'designed_positions'):
                continue
            if isinstance(v, (dict, list, tuple)):
                continue
            extra.append(k)
    fieldnames = [k for k in core if any(k in r for r in rows)] + extra
    with open(outdir / 'report.csv', 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        for r in rows:
            writer.writerow(r)
    passing = [r for r in rows if r.get('pass') and
               r.get('name') != 'native']
    passing.sort(key=lambda r: (-(r.get('tm_score') or -1),
                                r.get('plddt') or -1))
    top = passing[:top_n]
    mds = []
    mds.append('# Soluble TM redesign report\n')
    mds.append(f'- target: {config.get("name", "")}')
    mds.append(f'- designs: {len(rows)}, passing: {len(passing)}')
    native_rows = [r for r in rows if r.get('name') == 'native']
    for r in native_rows:
        mds.append(
            f'- native control: ecd_rmsd={r.get("ecd_rmsd"):.2f}A '
            f'tm={r.get("tm_score"):.3f} '
            f'surf_hydro_des={r.get("surf_hydro_des"):.3f} '
            f'plddt={r.get("plddt"):.1f} '
            f'(baseline to calibrate the gates above)')
    if not passing and len(rows):
        from collections import Counter
        counts = Counter()
        for r in rows:
            for g in r.get('fail_gates', []):
                counts[g] += 1
        if counts:
            mds.append('- no passing designs; most frequent failed gates: '
                       + ', '.join(f'{g} ({n})' for g, n in counts.most_common()))
    if passing:
        mds.append('\n## Passing designs\n')
    for r in passing:
        mds.append(
            f'- **{r["name"]}** tm={r.get("tm_score"):.3f} '
            f'plddt={r.get("plddt"):.1f} ecd_rmsd={r.get("ecd_rmsd"):.2f}A '
            f'surf_hydro={r.get("surf_hydro"):.3f} '
            f'mutations={r.get("mutations")}')
    Path(outdir / 'report.md').write_text('\n'.join(mds) + '\n')
    records = []
    for r in top:
        name = r['name']
        for c, seq in r['chains_full'].items():
            records.append((f'{name}_{c}', seq))
    from .pdbio import write_fasta
    write_fasta(outdir / 'variants_top.fa', records)
    for r in passing:
        for c, seq in r['chains_full'].items():
            records.append((f'{r["name"]}_{c}_pass', seq))
    write_fasta(outdir / 'variants_all_passing.fa',
                [(f'{r["name"]}_{c}_pass', r['chains_full'][c])
                 for r in passing for c in r['chains_full']])
    summary = {
        'config': config.get('name'),
        'n_total': len(rows),
        'n_passing': len(passing),
        'filters': config.get('filters', {}),
        'variants_top': [{'name': r['name'], 'seq': r['chains_full'][list(r['chains_full'])[0]]}
                         for r in top],
    }
    Path(outdir / 'summary.json').write_text(
        json.dumps(summary, indent=2, default=str))