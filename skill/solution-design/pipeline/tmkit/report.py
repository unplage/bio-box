import json
import csv
from pathlib import Path

import numpy as np

REPORT_SCHEMA_VERSION = '2.0'


def _is_missing(value):
    return value is None or (isinstance(value, float) and np.isnan(value))


def _as_number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _cmp(value, op, threshold):
    if _is_missing(value):
        return False
    a = _as_number(value)
    b = _as_number(threshold)
    if a is None or b is None:
        if op == '==':
            return value == threshold
        if op == '!=':
            return value != threshold
        return False
    if op == '>':
        return a > b
    if op == '<':
        return a < b
    if op == '>=':
        return a >= b
    if op == '<=':
        return a <= b
    if op == '==':
        return a == b
    if op == '!=':
        return a != b
    return False


def _eval_value(value, cond):
    """Evaluate one value against a condition (no cross-key access)."""
    if isinstance(cond, dict):
        if 'any' in cond:
            return any(_eval_value(value, c) for c in cond['any'])
        if 'all' in cond:
            return all(_eval_value(value, c) for c in cond['all'])
        if 'between' in cond:
            lo, hi = cond['between']
            return _cmp(value, '>=', lo) and _cmp(value, '<=', hi)
        if 'in' in cond:
            return value in cond['in']
        if cond.get('not_nan'):
            return not _is_missing(value)
        if 'op' in cond:
            return _cmp(value, cond['op'], cond.get('value'))
        return False
    if isinstance(cond, (list, tuple)):
        if cond and isinstance(cond[0], (list, tuple, dict)):
            return all(_eval_value(value, c) for c in cond)
        if cond and isinstance(cond[0], str):
            op = cond[0]
            if op in ('between', 'range') and len(cond) >= 3:
                return (_cmp(value, '>=', cond[1])
                        and _cmp(value, '<=', cond[2]))
            if op in ('in', 'any_of') and len(cond) == 2:
                return value in cond[1]
            if op == 'not_nan':
                return not _is_missing(value)
            if len(cond) == 2:
                return _cmp(value, op, cond[1])
            if len(cond) == 4:
                return (_cmp(value, cond[0], cond[1])
                        and _cmp(value, cond[2], cond[3]))
        return value in cond
    if isinstance(cond, (int, float)):
        return _cmp(value, '>=', cond)
    if isinstance(cond, str):
        return value == cond
    return False


def _eval_condition(row, cond, default_value):
    """Evaluate a condition that may reference other row keys."""
    if isinstance(cond, dict):
        if 'any' in cond:
            return any(_eval_condition(row, c, default_value)
                       for c in cond['any'])
        if 'all' in cond:
            return all(_eval_condition(row, c, default_value)
                       for c in cond['all'])
        if 'key' in cond:
            value = row.get(cond['key'])
            sub = {k: v for k, v in cond.items() if k != 'key'}
            return _eval_value(value, sub)
    return _eval_value(default_value, cond)


def _row_chains(row):
    import re
    chains = set()
    for key in row:
        match = re.match(r'mutations_(.+)$', key)
        if match:
            chains.add(match.group(1))
        match = re.match(r'interface_recovery_(.+)$', key)
        if match:
            chains.update(match.group(1).split('_'))
    return sorted(chains)


def passes_filters(row, filters):
    if row.get('predicted') is False:
        row['pass'] = None
        row['fail_gates'] = []
        return None
    checks = {}
    for key, cond in (filters or {}).items():
        if isinstance(key, (list, tuple)):
            keys = list(key)
            mode = cond.get('mode', 'all') if isinstance(cond, dict) else 'all'
            results = [_eval_condition(row, cond, row.get(k)) for k in keys]
            checks[str(key)] = any(results) if mode == 'any' else all(results)
            continue
        if '{chain}' in key:
            chains = _row_chains(row)
            if not chains:
                checks[key] = False
                continue
            mode = cond.get('mode', 'all') if isinstance(cond, dict) else 'all'
            results = [_eval_condition(
                row, cond, row.get(key.format(chain=c))) for c in chains]
            checks[key] = any(results) if mode == 'any' else all(results)
            continue
        checks[key] = _eval_condition(row, cond, row.get(key))
    row['pass'] = all(checks.values()) if checks else False
    row['fail_gates'] = [k for k, v in checks.items() if not v]
    return row['pass']


def write_tagged_pdb(target, chain_ids, mask, outpath):
    """Write the harmonized target with B=1 on designable reference positions."""
    lines = []
    serial = 1
    for cid in chain_ids:
        data = target.chains.get(cid)
        if data is None:
            continue
        mapping = target.reference_map.get(cid, {})
        for residue in data.residues:
            if residue.id[0] != ' ':
                continue
            refn = mapping.get(residue.id[1], residue.id[1])
            info = mask.get(cid, {}).get(refn)
            b = 1.0 if info and info.get('designable') else 0.0
            for atom in residue:
                coord = atom.get_coord()
                lines.append(
                    f'ATOM  {serial:5d} {atom.name:>4} {residue.resname:>3} '
                    f'{cid:>1}{residue.id[1]:4d}    '
                    f'{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}'
                    f'{1.0:6.2f}{b:6.2f}          {atom.element:>2}')
                serial += 1
    Path(outpath).write_text('\n'.join(lines) + 'END\n')


def _fmt(value, spec='.3f'):
    if isinstance(value, (int, float)) and not (
            isinstance(value, float) and np.isnan(value)):
        return f'{value:{spec}}'
    return 'NA'


def write_report(rows, config, outdir, top_n=10):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    core = ['name', 'pass', 'predicted', 'status', 'fail_gates',
            'tm_score', 'complex_tm_score', 'plddt',
            'plddt_ecd', 'plddt_tm', 'plddt_icd', 'plddt_loop',
            'ecd_rmsd', 'ecd_internal_rmsd', 'ecd_centroid_shift',
            'ecd_tm_angle', 'ligand_rmsd',
            'surf_hydro', 'surf_hydro_des', 'surf_hydro_tm', 'plddt_ecd',
            'pocket_rmsd', 'iptm', 'ptm', 'chain_symmetry_rmsd',
            'interface_recovery', 'exposed_arom', 'cys_des', 'seq_recovery',
            'score', 'mutations', 'soluble_score', 'residual_hydro',
            'over_polar', 'interface_mut', 'net_charge', 'abs_net_charge',
            'isoelectric_point', 'hydrophobic_moment', 'max_hydro_window',
            'max_charge_run', 'max_polar_window',
            'model']
    extra = []
    for r in rows:
        for k, v in r.items():
            if k in core or k in extra or k in ('chains_full', 'designed_positions'):
                continue
            if isinstance(v, (dict, list, tuple)):
                continue
            extra.append(k)
    fieldnames = list(dict.fromkeys(
        [k for k in core if any(k in r for r in rows)] + extra))
    with open(outdir / 'report.csv', 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        for r in rows:
            writer.writerow(r)
    passing = [r for r in rows if r.get('pass') and
               r.get('name') != 'native']

    def _sort_num(value):
        return value if isinstance(value, (int, float)) and value == value \
            else -1.0

    passing.sort(key=lambda r: (-_sort_num(r.get('tm_score')),
                                -_sort_num(r.get('plddt'))))
    top = passing[:top_n]
    mds = []
    mds.append('# Soluble TM redesign report\n')
    mds.append(f'- schema_version: {REPORT_SCHEMA_VERSION}')
    mds.append(f'- target: {config.get("name", "")}')
    missing = [r for r in rows if r.get('predicted') is False]
    mds.append(f'- designs: {len(rows)}, evaluated: {len(rows) - len(missing)}, '
               f'passing: {len(passing)}')
    if missing:
        mds.append(f'- **missing predictions: {len(missing)} design(s)** '
                   f'(not counted as failures; see report.csv status=no_model)')
    native_rows = [r for r in rows if r.get('name') == 'native']
    for r in native_rows:
        mds.append(
            f'- native control: ecd_rmsd={_fmt(r.get("ecd_rmsd"), ".2f")} '
            f'tm={_fmt(r.get("tm_score"))} '
            f'surf_hydro_des={_fmt(r.get("surf_hydro_des"))} '
            f'plddt={_fmt(r.get("plddt"), ".1f")} '
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
            f'- **{r["name"]}** tm={_fmt(r.get("tm_score"))} '
            f'plddt={_fmt(r.get("plddt"), ".1f")} '
            f'ecd_rmsd={_fmt(r.get("ecd_rmsd"), ".2f")} '
            f'surf_hydro={_fmt(r.get("surf_hydro"))} '
            f'mutations={_fmt(r.get("mutations"), "d") if isinstance(r.get("mutations"), int) else r.get("mutations")}')
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
        'schema_version': REPORT_SCHEMA_VERSION,
        'config': config.get('name'),
        'n_total': len(rows),
        'n_evaluated': len(rows) - len(missing),
        'n_missing_models': len(missing),
        'n_passing': len(passing),
        'filters': config.get('filters', {}),
        'variants_top': [{'name': r['name'], 'seq': r['chains_full'][list(r['chains_full'])[0]]}
                         for r in top],
    }
    Path(outdir / 'summary.json').write_text(
        json.dumps(summary, indent=2, default=str))