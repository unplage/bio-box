import re
import subprocess
from pathlib import Path

import numpy as np
import yaml
from Bio.PDB.PDBParser import PDBParser
from .pdbio import write_fasta, residue_sasa_from_file, ca_coords, \
    map_model_to_reference, as_pdb, hetero_residues, MAX_ASA
from .metrics import tm_score, parse_plddt, rmsd, kabsch, transform, \
    centroid, angle_between, normalize_plddt_scale, \
    surface_hydrophobic_fraction, surface_hydrophobic_fraction_at, \
    exposed_aromatics_at, cys_at, mean_plddt_at, mutation_calls
from .confidence import load_confidence

NAN = float('nan')

# Predictor output naming conventions. `{name}` is the design name.
DEFAULT_RESULT_PATTERNS = (
    '{name}*rank_001*.pdb', '{name}*rank_1*.pdb', '{name}*_pred_0*.pdb',
    '{name}*_model_0.cif', '{name}*_model_0.pdb',
    '{name}*_model.cif', '{name}*_model.pdb',
    '{name}*/ranked_0.pdb', '{name}*ranked_0.pdb',
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
    """Build a predictor command for one or more design FASTA files.

    ColabFold takes a directory (or a single FASTA); AlphaFold2 needs a
    space-separated ``--fasta_paths`` argument.  Boltz-2 is driven from the
    per-design YAML written by ``write_external_inputs``.  External tools can
    still be run by the user and supplied through ``results_dir``.
    """
    tool = config.get('predict_tool', 'colabfold')
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    paths = [Path(f) for f in (fasta if isinstance(fasta, (list, tuple))
                               else [fasta])]
    if not paths:
        return None
    if tool == 'colabfold':
        # A directory is the only reliable way to give ColabFold many jobs;
        # a single file is also accepted.
        input_arg = str(paths[0]) if len(paths) == 1 else str(paths[0].parent)
        exe = config.get('colabfold_executable') or \
            config.get('predict_executable') or 'colabfold_batch'
        cmd = [
            exe, input_arg,
            '--model-type', config.get('colabfold_model', 'alpha2_multimer_v3'),
            '--num-recycle', str(config.get('num_recycles', 3)),
            '--outdir', str(outdir),
        ]
        if config.get('colabfold_amber', False):
            cmd.append('--amber')
        if config.get('colabfold_templates', False):
            cmd.append('--templates')
        return cmd
    if tool == 'alphafold':
        cmd = [
            'python', str(config['alphafold_script']),
            '--fasta_paths', ' '.join(str(p) for p in paths),
            '--output_dir', str(outdir),
            '--model_preset', 'multimer',
            '--db_preset', config.get('db_preset', 'full_dbs'),
            '--max_template_date',
            str(config.get('max_template_date', '2020-05-14')),
        ]
        for key, flag in [('data_dir', '--data_dir'),
                          ('uniref90', '--uniref90_database_path'),
                          ('mgnify', '--mgnify_database_path'),
                          ('bfd', '--bfd_database_path'),
                          ('uniref30', '--uniref30_database_path'),
                          ('pdb_seqres', '--pdb_seqres_database_path'),
                          ('template_mmcif', '--template_mmcif_dir'),
                          ('obsolete', '--obsolete_pdbs_path')]:
            if config.get(key):
                cmd += [flag, str(config[key])]
        return cmd
    if tool in ('boltz', 'boltz2'):
        boltz_in = config.get('boltz_input_dir') or \
            (paths[0].parent / 'boltz')
        exe = config.get('boltz_executable') or \
            config.get('predict_executable') or 'boltz'
        cmd = [exe, 'predict', str(boltz_in),
               '--out_dir', str(outdir)]
        if config.get('boltz_use_msa_server', False):
            cmd.append('--use_msa_server')
        if config.get('boltz_no_kernels', False):
            cmd.append('--no_kernels')
        if config.get('boltz_cache'):
            cmd += ['--cache', str(config['boltz_cache'])]
        if config.get('boltz_devices'):
            cmd += ['--devices', str(config['boltz_devices'])]
        return cmd
    # External predictors (af3 / chai / ...) are run by the user; results
    # are supplied via `results_dir` and discovered by filename.
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
            # Accept name-prefixed files (ColabFold/Boltz) and AF2's
            # <name>/ranked_0.pdb layout where the design name is the parent.
            in_parent = any(part == name or part.startswith(f'{name}_')
                            or part.startswith(f'{name}-')
                            for part in cand.parts[:-1])
            if not (base.startswith(name) or in_parent):
                continue
            if base.startswith(name) and len(base) > len(name) \
                    and base[len(name)].isalnum():
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


def _ligand_rmsd(model_pdb, target, cofactors, R, t):
    """RMSD of matching hetero residues after the protein superposition."""
    if not cofactors:
        return NAN
    try:
        model = PDBParser(QUIET=True).get_structure(
            'model', str(as_pdb(model_pdb)))[0]
    except Exception:
        return NAN

    def keyed(items):
        out = {}
        counts = {}
        for item in items:
            key = item['resname'].upper()
            counts[key] = counts.get(key, 0) + 1
            out[(key, counts[key])] = item['coords']
        return out

    values = []
    for spec in cofactors:
        resnames = spec.get('resnames') if isinstance(spec, dict) else [spec]
        if not resnames:
            continue
        template = keyed(hetero_residues(target.structure, resnames))
        predicted = keyed(hetero_residues(model, resnames))
        for key in set(template) & set(predicted):
            a = template[key]
            b = transform(predicted[key], R, t)
            values.append(rmsd(b, a) if a.shape == b.shape
                          else float(np.linalg.norm(
                              centroid(b) - centroid(a))))
    return float(np.mean(values)) if values else NAN


def evaluate_models(model_paths, design_records, chain_ids, target, mask,
                    ecd_lookup, ref_seqs, pred_dir, pocket_residues=None,
                    tm_segments=None, pocket_shell_dist=8.0,
                    interface_recovery_dist=8.0, topologies=None,
                    symmetry_groups=None, interface_pairs=None,
                    sasa_threshold=40.0, normalize_plddt=True,
                    relative_sasa=False, cofactors=None):
    """Evaluate predicted models against the template.

    All residue pairing is alignment-driven: model residue numbers are
    mapped to reference positions via ``map_model_to_reference`` (robust to
    renumbering/gaps), then compared to template C-alpha at common reference
    positions.  PDB and mmCIF are accepted.

    ``ecd_rmsd`` is measured **after** a TMD/whole-complex superposition, so
    it detects rigid-body displacement of the ECD relative to the TMD.
    ``ecd_internal_rmsd`` is the old ECD-only Kabsch value and is reported
    separately for fold quality, not as the ECD-fixed gate.
    """
    rows = []
    mask = mask or {}
    topologies = topologies or {}
    pocket_seed = set(int(x) for x in (pocket_residues or []))
    groups = [list(g) for g in (symmetry_groups or []) if len(g) > 1]
    if not groups:
        seqs = [ref_seqs.get(c) for c in chain_ids]
        if len(chain_ids) > 1 and all(s and s == seqs[0] for s in seqs):
            groups = [list(chain_ids)]
    pairs = [tuple(p) for p in (interface_pairs or []) if len(p) == 2]
    if not pairs and len(chain_ids) > 1:
        pairs = [(chain_ids[i], chain_ids[j])
                 for i in range(len(chain_ids))
                 for j in range(i + 1, len(chain_ids))]

    def _segments(cid):
        topo = topologies.get(cid) or {}
        ecd_segs = topo.get('ecd_segments')
        if ecd_segs is None:
            ecd_segs = ecd_lookup or []
        tm_segs = topo.get('tm_segments')
        if tm_segs is None:
            tm_segs = tm_segments or []
        ecd = {r for a, b in ecd_segs for r in range(int(a), int(b) + 1)}
        tm = {r for a, b in tm_segs for r in range(int(a), int(b) + 1)}
        icd = {r for a, b in (topo.get('icd_segments') or [])
               for r in range(int(a), int(b) + 1)}
        return ecd, tm, icd

    for rec in design_records:
        name = rec['name'].replace(':', '_')
        model_file = model_paths.get(name)
        if model_file is None:
            row = rec.copy()
            row.update({'tm_score': NAN, 'complex_tm_score': NAN,
                        'plddt': NAN, 'plddt_ecd': NAN, 'plddt_tm': NAN,
                        'plddt_icd': NAN, 'plddt_loop': NAN,
                        'ecd_rmsd': NAN, 'ecd_internal_rmsd': NAN,
                        'ecd_centroid_shift': NAN, 'ecd_tm_angle': NAN,
                        'surf_hydro': NAN, 'surf_hydro_des': NAN,
                        'surf_hydro_tm': NAN, 'pocket_rmsd': NAN,
                        'ligand_rmsd': NAN, 'exposed_arom': NAN,
                        'cys_des': NAN, 'mutations': NAN, 'iptm': NAN,
                        'ptm': NAN, 'chain_symmetry_rmsd': NAN,
                        'interface_recovery': NAN, 'model': '',
                        'predicted': False, 'status': 'no_model'})
            for cid in chain_ids:
                row[f'mutations_{cid}'] = NAN
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
        if relative_sasa:
            all_sasa = {
                key: area / MAX_ASA.get(aa, 1.0)
                for key, area in all_sasa.items()
                for aa in [aa_map.get(key, 'X')]
            }

        conf = load_confidence(model_file)
        conf_plddt = conf.get('plddt')
        bfactors = [v for v in plddt.values()
                    if isinstance(v, (int, float)) and v == v]
        if conf_plddt is not None and (not bfactors or max(bfactors) < 1e-6):
            value = float(conf_plddt)
            if value <= 1.5:
                value *= 100.0
            plddt = {key: value for key in aa_map}
        elif normalize_plddt:
            plddt = normalize_plddt_scale(plddt)

        model_by_ref = {}
        template_by_ref = {}
        common_by_chain = {}
        region_refs = {}
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
            common_by_chain[cid] = common
            region_refs[cid] = _segments(cid)

        # Global superposition: TMD if available, otherwise the whole complex.
        align_m, align_t = [], []
        for cid in chain_ids:
            ecd, tm, icd = region_refs[cid]
            use = [r for r in common_by_chain[cid] if r in tm]
            if not use:
                use = list(common_by_chain[cid])
            align_m.extend(model_by_ref[cid][r] for r in use)
            align_t.extend(template_by_ref[cid][r] for r in use)
        if len(align_m) >= 3:
            R, t = kabsch(np.asarray(align_m), np.asarray(align_t))
        else:
            R = np.eye(3)
            t = np.zeros(3)

        ligand_rmsd = _ligand_rmsd(model_pdb, target, cofactors, R, t)

        tm_total = []
        complex_m, complex_t = [], []
        for cid in chain_ids:
            common = common_by_chain[cid]
            if len(common) >= 5:
                m = np.asarray([model_by_ref[cid][r] for r in common])
                tt = np.asarray([template_by_ref[cid][r] for r in common])
                tm_total.append(tm_score(m, tt))
                complex_m.extend(m)
                complex_t.extend(tt)

        ecd_vals, ecd_internal, ecd_model_all, ecd_template_all = [], [], [], []
        pocket_rmsd_vals = []
        pocket = _expand_pocket_shell(template_by_ref, pocket_seed,
                                      pocket_shell_dist)
        for cid in chain_ids:
            ecd, tm, icd = region_refs[cid]
            common = common_by_chain[cid]
            ecd_common = [r for r in common if r in ecd]
            if len(ecd_common) >= 3:
                m = np.asarray([model_by_ref[cid][r] for r in ecd_common])
                tt = np.asarray([template_by_ref[cid][r] for r in ecd_common])
                ecd_internal.append(rmsd(m, tt))
                moved = transform(m, R, t)
                ecd_vals.append(float(np.sqrt(
                    ((moved - tt) ** 2).sum(1).mean())))
                ecd_model_all.extend(moved)
                ecd_template_all.extend(tt)
            p_common = [r for r in common if r in pocket]
            if len(p_common) >= 2:
                m = np.asarray([model_by_ref[cid][r] for r in p_common])
                tt = np.asarray([template_by_ref[cid][r] for r in p_common])
                moved = transform(m, R, t)
                pocket_rmsd_vals.append(float(np.sqrt(
                    ((moved - tt) ** 2).sum(1).mean())))

        ecd_centroid_shift = NAN
        ecd_tm_angle = NAN
        if ecd_model_all and ecd_template_all:
            ecd_centroid_shift = float(np.linalg.norm(
                centroid(ecd_model_all) - centroid(ecd_template_all)))
            if len(align_m) >= 3:
                ecd_tm_angle = angle_between(
                    centroid(ecd_model_all) - transform(centroid(align_m), R, t),
                    centroid(ecd_template_all) - centroid(align_t))

        des_pos = []
        tm_pos = []
        plddt_pos = {'ecd': [], 'tm': [], 'icd': [], 'loop': []}
        for cid in chain_ids:
            ecd, tm, icd = region_refs[cid]
            refn_to_model = {v: k for k, v in
                             model_res_maps.get(cid, {}).items()}
            for refn, info in (mask.get(cid) or {}).items():
                if not info.get('designable') or refn not in refn_to_model:
                    continue
                des_pos.append((cid, refn_to_model[refn]))
            for refn, model_res in refn_to_model.items():
                if refn in ecd:
                    plddt_pos['ecd'].append((cid, model_res))
                elif refn in tm:
                    tm_pos.append((cid, model_res))
                    plddt_pos['tm'].append((cid, model_res))
                elif refn in icd:
                    plddt_pos['icd'].append((cid, model_res))
                else:
                    plddt_pos['loop'].append((cid, model_res))

        sym_vals = []
        for group in groups:
            if not all(model_by_ref.get(c) for c in group):
                continue
            shared = sorted(set.intersection(
                *[set(model_by_ref[c]) for c in group]))
            if len(shared) < 5:
                continue
            for i in range(len(group)):
                for j in range(i + 1, len(group)):
                    sym_vals.append(rmsd(
                        np.asarray([model_by_ref[group[i]][r]
                                    for r in shared]),
                        np.asarray([model_by_ref[group[j]][r]
                                    for r in shared])))
        sym_rmsd = float(np.mean(sym_vals)) if sym_vals else NAN

        interface_recovery = NAN
        pair_recovery = {}
        if len(chain_ids) > 1 and len(template_by_ref) == len(chain_ids)                 and len(model_by_ref) == len(chain_ids):
            t_contacts = _interface_contacts(template_by_ref, chain_ids,
                                             interface_recovery_dist)
            m_contacts = _interface_contacts(model_by_ref, chain_ids,
                                             interface_recovery_dist)
            if t_contacts:
                interface_recovery = len(t_contacts & m_contacts) / len(t_contacts)
            for a, b in pairs:
                ta = {x for x in t_contacts
                      if {x[0], x[2]} == {a, b}}
                ma = {x for x in m_contacts
                      if {x[0], x[2]} == {a, b}}
                pair_recovery[f'{a}_{b}'] = (
                    len(ta & ma) / len(ta)) if ta else NAN

        mut_total = 0
        mut_by_chain = {}
        for cid in chain_ids:
            ref = ref_seqs.get(cid, '')
            des = rec.get('chains_full', {}).get(cid, '')
            n = len(mutation_calls(ref, des)) if ref and des else 0
            mut_by_chain[cid] = n
            mut_total += n

        plddt_vals = list(plddt.values()) if plddt else []
        row = rec.copy()
        row.update({
            'tm_score': float(np.mean(tm_total)) if tm_total else NAN,
            'complex_tm_score': tm_score(np.asarray(complex_m),
                                         np.asarray(complex_t))
            if len(complex_m) >= 5 else NAN,
            'ecd_rmsd': float(np.mean(ecd_vals)) if ecd_vals else NAN,
            'ecd_internal_rmsd': float(np.mean(ecd_internal))
            if ecd_internal else NAN,
            'ecd_centroid_shift': ecd_centroid_shift,
            'ecd_tm_angle': ecd_tm_angle,
            'plddt': float(np.mean(plddt_vals)) if plddt_vals else NAN,
            'plddt_ecd': mean_plddt_at(plddt, plddt_pos['ecd']),
            'plddt_tm': mean_plddt_at(plddt, plddt_pos['tm']),
            'plddt_icd': mean_plddt_at(plddt, plddt_pos['icd']),
            'plddt_loop': mean_plddt_at(plddt, plddt_pos['loop']),
            'surf_hydro': surface_hydrophobic_fraction(
                all_sasa, aa_map, sasa_threshold=sasa_threshold),
            'surf_hydro_des': surface_hydrophobic_fraction_at(
                all_sasa, aa_map, des_pos, sasa_threshold=sasa_threshold),
            'surf_hydro_tm': surface_hydrophobic_fraction_at(
                all_sasa, aa_map, tm_pos, sasa_threshold=sasa_threshold),
            'pocket_rmsd': float(np.mean(pocket_rmsd_vals))
            if pocket_rmsd_vals else NAN,
            'ligand_rmsd': ligand_rmsd,
            'exposed_arom': exposed_aromatics_at(all_sasa, aa_map, des_pos,
                                                 sasa_threshold=sasa_threshold),
            'cys_des': cys_at(aa_map, des_pos),
            'mutations': mut_total,
            'iptm': conf['iptm'] if conf['iptm'] is not None else NAN,
            'ptm': conf['ptm'] if conf['ptm'] is not None else NAN,
            'chain_symmetry_rmsd': sym_rmsd,
            'interface_recovery': interface_recovery,
            'model': str(model_file),
            'predicted': True, 'status': 'evaluated',
        })
        for cid, n in mut_by_chain.items():
            row[f'mutations_{cid}'] = n
        for key, value in pair_recovery.items():
            row[f'interface_recovery_{key}'] = value
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