import re
import subprocess
from pathlib import Path

import numpy as np
import yaml
from .pdbio import write_fasta, residue_sasa_from_file, ca_coords, \
    map_model_to_reference, as_pdb
from .metrics import tm_score, parse_plddt, rmsd, \
    surface_hydrophobic_fraction, surface_hydrophobic_fraction_at, \
    exposed_aromatics_at, cys_at, mean_plddt_at, mutation_calls
from .confidence import load_confidence

NAN = float('nan')

# Predictor output naming conventions. `{name}` is the design name.
DEFAULT_RESULT_PATTERNS = (
    '{name}*rank_001*.pdb', '{name}*rank_1*.pdb', '{name}*_pred_0*.pdb',
    '{name}*_model_0.cif', '{name}*_model_0.pdb',
    '{name}*_model.cif', '{name}*_model.pdb',
)


def _model_rank(path):
    name = path.name
    match = re.search(r'rank[_-]?(\d+)', name)
    if match:
        return (0, int(match.group(1)))
    match = re.search(r'_model[_-]?(\d+)', name)
    if match:
        return (1, int(match.group(1)))
    return (2, 0)


def multimer_fasta(assemble_records, chain_ids, outdir):
    """Write one FASTA file per design; each file has one ``>`` record per
    chain (AF2 multimer format: all records in a single file = one complex)."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    paths = []
    for rec in assemble_records:
        name = rec['name'].replace(':', '_')
        fasta_records = [(f'{name}_{c}', rec['chains_full'][c])
                         for c in chain_ids]
        fa_path = outdir / f'{name}.fa'
        write_fasta(fa_path, fasta_records)
        paths.append(fa_path)
    return paths


def write_external_inputs(assemble_records, chain_ids, outdir):
    """Emit predictor-native inputs next to the AF2 FASTA.

    Boltz-2 prefers YAML with explicit chain ids (``boltz/<name>.yaml``);
    AlphaFold3 uses a JSON spec (``af3/<name>.json``). Both keep the design
    name as the job name so predicted outputs are discoverable by
    ``find_rank1_models``.
    """
    import json
    outdir = Path(outdir)
    boltz_dir = outdir / 'boltz'
    af3_dir = outdir / 'af3'
    boltz_dir.mkdir(parents=True, exist_ok=True)
    af3_dir.mkdir(parents=True, exist_ok=True)
    for rec in assemble_records:
        name = rec['name'].replace(':', '_')
        sequences = [{'protein': {'id': c,
                                  'sequence': rec['chains_full'][c]}}
                     for c in chain_ids]
        payload = {'version': 1, 'sequences': sequences}
        (boltz_dir / f'{name}.yaml').write_text(
            yaml.safe_dump(payload, sort_keys=False))
        af3 = {'name': name, 'modelSeeds': [1], 'dialect': 'alphafold3',
               'version': 1, 'sequences': sequences}
        (af3_dir / f'{name}.json').write_text(json.dumps(af3, indent=2))


def predict_command(config, fasta, outdir, dry_run=True):
    tool = config.get('predict_tool', 'colabfold')
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    if isinstance(fasta, (list, tuple)):
        fasta_str = ','.join(str(f) for f in fasta)
    else:
        fasta_str = str(fasta)
    if tool == 'colabfold':
        cmd = [
            'colabfold_batch', fasta_str,
            '--model-type', config.get('colabfold_model', 'alpha2_multimer_v3'),
            '--num-recycle', str(config.get('num_recycles', 3)),
            '--outdir', str(outdir),
        ]
        if config.get('colabfold_amber', False):
            cmd.append('--amber')
        return cmd
    if tool == 'alphafold':
        cmd = [
            'python', str(config['alphafold_script']),
            '--fasta_paths', fasta_str,
            '--output_dir', str(outdir),
            '--model_preset', 'multimer',
            '--db_preset', config.get('db_preset', 'full_dbs'),
        ]
        for key, flag in [('data_dir', '--data_dir'), ('uniref90', '--uniref90_database_path'),
                          ('mgnify', '--mgnify_database_path'), ('bfd', '--bfd_database_path'),
                          ('uniref30', '--uniref30_database_path'), ('pdb_seqres', '--pdb_seqres_database_path'),
                          ('template_mmcif', '--template_mmcif_dir'), ('obsolete', '--obsolete_pdbs_path')]:
            if config.get(key):
                cmd += [flag, str(config[key])]
        return cmd
    # External predictors (boltz / af3 / chai / ...) are run by the user;
    # results are supplied via `results_dir` and discovered by filename.
    return None


def find_rank1_models(pred_dir, names, patterns=None):
    """Locate the best (rank-1 / model-0) structure per design.

    Recursive so nested layouts (e.g. ``boltz_results_*/predictions/<name>/``
    or AF3 ``<name>/<name>_model.cif``) work. Supports PDB and mmCIF.
    """
    pred_dir = Path(pred_dir)
    patterns = patterns or DEFAULT_RESULT_PATTERNS
    found = {}
    for name in names:
        candidates = set()
        for pattern in patterns:
            candidates.update(pred_dir.rglob(pattern.format(name=name)))
        keep = []
        for cand in candidates:
            base = cand.name
            if not base.startswith(name):
                continue
            if len(base) > len(name) and base[len(name)].isalnum():
                continue
            keep.append(cand)
        if keep:
            found[name] = sorted(keep, key=_model_rank)[0]
    return found


def _interface_contacts(by_chain, chain_ids, dist):
    """Cross-chain C-alpha contact set keyed by (chain, refnum) pairs."""
    pairs = set()
    cutoff2 = float(dist) ** 2
    for i in range(len(chain_ids)):
        for j in range(i + 1, len(chain_ids)):
            ci, cj = chain_ids[i], chain_ids[j]
            a, b = by_chain.get(ci, {}), by_chain.get(cj, {})
            if not a or not b:
                continue
            ra, rb = list(a), list(b)
            ca = np.asarray([a[r] for r in ra], dtype=float)
            cb = np.asarray([b[r] for r in rb], dtype=float)
            d2 = ((ca[:, None, :] - cb[None, :, :]) ** 2).sum(-1)
            ii, jj = np.where(d2 <= cutoff2)
            for x, y in zip(ii, jj):
                pairs.add((ci, ra[x], cj, rb[y]))
    return pairs


def _expand_pocket_shell(template_by_ref, pocket_seed, dist):
    """Grow the pocket residue set with template C-alpha neighbours."""
    shell = set(pocket_seed)
    if not pocket_seed:
        return shell
    seeds = np.asarray([c for refs in template_by_ref.values()
                        for r, c in refs.items() if r in pocket_seed])
    if len(seeds) == 0:
        return shell
    cutoff2 = float(dist) ** 2
    for refs in template_by_ref.values():
        for refn, coord in refs.items():
            if refn in shell:
                continue
            if ((seeds - coord) ** 2).sum(1).min() <= cutoff2:
                shell.add(refn)
    return shell


def evaluate_models(model_paths, design_records, chain_ids, target, mask,
                    ecd_lookup, ref_seqs, pred_dir, pocket_residues=None,
                    tm_segments=None, pocket_shell_dist=8.0,
                    interface_recovery_dist=8.0):
    """Evaluate predicted models against the template.

    All residue pairing is alignment-driven: model residue numbers are
    mapped to reference positions via `map_model_to_reference` (robust to
    renumbering/gaps), then compared to template C-alpha at common
    reference positions. Accepts PDB or mmCIF (Boltz-2 / AlphaFold3).
    """
    rows = []
    mask = mask or {}
    pocket_seed = set(int(x) for x in (pocket_residues or []))
    _ecds = [r for a, b in ecd_lookup for r in range(a, b + 1)]
    _tm = set(r for a, b in (tm_segments or [])
              for r in range(int(a), int(b) + 1))
    for rec in design_records:
        name = rec['name'].replace(':', '_')
        model_file = model_paths.get(name)
        if model_file is None:
            row = rec.copy()
            row.update({'tm_score': NAN, 'plddt': NAN, 'ecd_rmsd': NAN,
                        'surf_hydro': NAN, 'surf_hydro_des': NAN,
                        'surf_hydro_tm': NAN, 'plddt_ecd': NAN,
                        'pocket_rmsd': NAN, 'exposed_arom': NAN,
                         'cys_des': NAN, 'mutations': NAN, 'iptm': NAN,
                         'ptm': NAN, 'chain_symmetry_rmsd': NAN,
                         'interface_recovery': NAN, 'model': ''})
            rows.append(row)
            continue
        model_pdb = as_pdb(model_file)
        model_res_maps = {c: map_model_to_reference(model_pdb, c, ref_seqs[c])
                          for c in chain_ids}
        plddt = parse_plddt(model_pdb)
        all_sasa = residue_sasa_from_file(model_pdb)
        aa_map = {}
        with open(model_pdb) as f:
            for line in f:
                if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                    ch = line[21]
                    res = int(line[22:26])
                    resname = line[17:20].strip()
                    aa_map[(ch, res)] = _aa1(resname)
        tm_total = []
        ecd_vals = []
        des_pos = []
        ecd_pos = []
        model_by_ref = {}
        template_by_ref = {}
        for cid in chain_ids:
            coords_t, refs_t_raw = ca_coords(target, cid)
            refs_t = [target.reference_map.get(cid, {}).get(r, r)
                      for r in refs_t_raw]
            if len(coords_t) == 0:
                continue
            t_by_ref = {r: c for r, c in zip(refs_t, coords_t)}
            template_by_ref[cid] = t_by_ref
            m_map = model_res_maps.get(cid, {})
            m_coords, m_refs_raw = ca_coords_from_pdb(model_pdb, cid)
            m_by_ref = {}
            for r, c in zip(m_refs_raw, m_coords):
                refn = m_map.get(r)
                if refn is not None:
                    m_by_ref[refn] = c
            model_by_ref[cid] = m_by_ref
            common = sorted(set(t_by_ref) & set(m_by_ref))
            if len(common) >= 5:
                tm_total.append(tm_score(
                    np.asarray([m_by_ref[r] for r in common]),
                    np.asarray([t_by_ref[r] for r in common])))
            ecd_common = [r for r in common if r in _ecds]
            if len(ecd_common) >= 3:
                ecd_vals.append(rmsd(
                    np.asarray([m_by_ref[r] for r in ecd_common]),
                    np.asarray([t_by_ref[r] for r in ecd_common])))
            refn_to_model = {v: k for k, v in m_map.items()}
            for refn, info in (mask.get(cid) or {}).items():
                if not info.get('designable') or refn not in refn_to_model:
                    continue
                des_pos.append((cid, refn_to_model[refn]))
            for refn in _ecds:
                if refn in refn_to_model:
                    ecd_pos.append((cid, refn_to_model[refn]))

        pocket = _expand_pocket_shell(template_by_ref, pocket_seed,
                                      pocket_shell_dist)
        pocket_rmsd_vals = []
        for cid in chain_ids:
            common = sorted(set(template_by_ref.get(cid, {}))
                            & set(model_by_ref.get(cid, {})) & pocket)
            if len(common) >= 2:
                pocket_rmsd_vals.append(rmsd(
                    np.asarray([model_by_ref[cid][r] for r in common]),
                    np.asarray([template_by_ref[cid][r] for r in common])))

        tm_pos = []
        for cid in chain_ids:
            refn_to_model = {v: k for k, v in model_res_maps.get(cid, {}).items()}
            for refn in _tm:
                if refn in refn_to_model:
                    tm_pos.append((cid, refn_to_model[refn]))

        sym_rmsd = NAN
        if len(chain_ids) > 1 and all(model_by_ref.get(c) for c in chain_ids):
            shared = sorted(set.intersection(
                *[set(model_by_ref[c]) for c in chain_ids]))
            if len(shared) >= 5:
                vals = []
                for i in range(len(chain_ids)):
                    for j in range(i + 1, len(chain_ids)):
                        vals.append(rmsd(
                            np.asarray([model_by_ref[chain_ids[i]][r]
                                        for r in shared]),
                            np.asarray([model_by_ref[chain_ids[j]][r]
                                        for r in shared])))
                sym_rmsd = float(np.mean(vals))

        interface_recovery = NAN
        if len(chain_ids) > 1 and len(template_by_ref) == len(chain_ids) \
                and len(model_by_ref) == len(chain_ids):
            t_contacts = _interface_contacts(template_by_ref, chain_ids,
                                              interface_recovery_dist)
            m_contacts = _interface_contacts(model_by_ref, chain_ids,
                                              interface_recovery_dist)
            if t_contacts:
                interface_recovery = len(t_contacts & m_contacts) / len(t_contacts)

        conf = load_confidence(model_file)
        plddt_vals = list(plddt.values()) if plddt else []
        ref = ref_seqs[chain_ids[0]]
        des = rec['chains_full'][chain_ids[0]]
        row = rec.copy()
        row.update({
            'tm_score': float(np.mean(tm_total)) if tm_total else NAN,
            'ecd_rmsd': float(np.mean(ecd_vals)) if ecd_vals else NAN,
            'plddt': float(np.mean(plddt_vals)) if plddt_vals else NAN,
            'surf_hydro': surface_hydrophobic_fraction(all_sasa, aa_map),
            'surf_hydro_des': surface_hydrophobic_fraction_at(
                all_sasa, aa_map, des_pos),
            'surf_hydro_tm': surface_hydrophobic_fraction_at(
                all_sasa, aa_map, tm_pos),
            'plddt_ecd': mean_plddt_at(plddt, ecd_pos),
            'pocket_rmsd': float(np.mean(pocket_rmsd_vals))
                if pocket_rmsd_vals else NAN,
            'exposed_arom': exposed_aromatics_at(all_sasa, aa_map, des_pos),
            'cys_des': cys_at(aa_map, des_pos),
            'mutations': len(mutation_calls(ref, des)),
            'iptm': conf['iptm'] if conf['iptm'] is not None else NAN,
            'ptm': conf['ptm'] if conf['ptm'] is not None else NAN,
            'chain_symmetry_rmsd': sym_rmsd,
            'interface_recovery': interface_recovery,
            'model': str(model_file),
        })
        rows.append(row)
    return rows


def _aa1(resname):
    from .pdbio import AMINO3
    return AMINO3.get(resname, 'X')


def ca_coords_from_pdb(pdb_path, chain_id):
    coords, refs = [], []
    with open(as_pdb(pdb_path)) as f:
        for line in f:
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                if line[21] == chain_id:
                    coords.append([float(line[30:38]), float(line[38:46]),
                                   float(line[46:54])])
                    refs.append(int(line[22:26]))
    return np.asarray(coords), refs