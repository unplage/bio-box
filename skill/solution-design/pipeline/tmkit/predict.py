import subprocess
from pathlib import Path

import numpy as np
from .pdbio import write_fasta, residue_sasa_from_file, ca_coords, \
    map_model_to_reference
from .metrics import tm_score, parse_plddt, rmsd, \
    surface_hydrophobic_fraction, surface_hydrophobic_fraction_at, \
    exposed_aromatics_at, cys_at, mean_plddt_at, mutation_calls


def multimer_fasta(assemble_records, chain_ids, outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    fasta_records = []
    for rec in assemble_records:
        seqs = [rec['chains_full'][c] for c in chain_ids]
        name = rec['name'].replace(':', '_')
        fasta_records.append((name, ':'.join(seqs)))
    write_fasta(outdir / 'predict_input.fasta', fasta_records)
    return outdir / 'predict_input.fasta'


def predict_command(config, fasta, outdir, dry_run=True):
    tool = config.get('predict_tool', 'colabfold')
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    if tool == 'colabfold':
        cmd = [
            'colabfold_batch', str(fasta),
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
            '--fasta_paths', str(fasta),
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
    return None


def find_rank1_models(pred_dir, names):
    pred_dir = Path(pred_dir)
    found = {}
    for name in names:
        hits = sorted(pred_dir.glob(name + '*rank_001*.pdb'))
        if hits:
            found[name] = hits[0]
        else:
            hits = sorted(pred_dir.glob(name + '*_pred_0.pdb'))
            if hits:
                found[name] = hits[0]
    return found


def evaluate_models(model_paths, design_records, chain_ids, target, mask,
                    ecd_lookup, ref_seqs, pred_dir, pocket_residues=None):
    """Evaluate predicted models against the template.

    All residue pairing is alignment-driven: model residue numbers are
    mapped to reference positions via `map_model_to_reference` (robust to
    renumbering/gaps), then compared to template C-alpha at common
    reference positions.
    """
    rows = []
    pocket = set(int(x) for x in (pocket_residues or []))
    _ecds = [r for a, b in ecd_lookup for r in range(a, b + 1)]
    for rec in design_records:
        name = rec['name'].replace(':', '_')
        model_file = model_paths.get(name)
        if model_file is None:
            row = rec.copy()
            row.update({'tm_score': float('nan'), 'plddt': float('nan'),
                        'ecd_rmsd': float('nan'), 'surf_hydro': float('nan'),
                        'surf_hydro_des': float('nan'),
                        'plddt_ecd': float('nan'), 'pocket_rmsd': float('nan'),
                        'exposed_arom': float('nan'), 'cys_des': float('nan'),
                        'mutations': float('nan'), 'model': ''})
            rows.append(row)
            continue
        model_res_maps = {c: map_model_to_reference(model_file, c, ref_seqs[c])
                          for c in chain_ids}
        plddt = parse_plddt(model_file)
        all_sasa = residue_sasa_from_file(model_file)
        aa_map = {}
        with open(model_file) as f:
            for line in f:
                if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                    ch = line[21]
                    res = int(line[22:26])
                    resname = line[17:20].strip()
                    aa_map[(ch, res)] = _aa1(resname)
        tm_total = []
        ecd_vals = []
        pocket_rmsd_vals = []
        des_pos = []
        ecd_pos = []
        pocket_pos = []
        for cid in chain_ids:
            coords_t, refs_t_raw = ca_coords(target, cid)
            refs_t = [target.reference_map.get(cid, {}).get(r, r)
                      for r in refs_t_raw]
            if len(coords_t) == 0:
                continue
            t_by_ref = {r: c for r, c in zip(refs_t, coords_t)}
            m_map = model_res_maps.get(cid, {})
            m_coords, m_refs_raw = ca_coords_from_pdb(model_file, cid)
            m_by_ref = {}
            for r, c in zip(m_refs_raw, m_coords):
                refn = m_map.get(r)
                if refn is not None:
                    m_by_ref[refn] = c
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
            pocket_common = [r for r in common if r in pocket]
            if len(pocket_common) >= 2:
                pocket_rmsd_vals.append(rmsd(
                    np.asarray([m_by_ref[r] for r in pocket_common]),
                    np.asarray([t_by_ref[r] for r in pocket_common])))
            refn_to_model = {v: k for k, v in m_map.items()}
            for refn, info in (mask.get(cid) or {}).items():
                if not info.get('designable') or refn not in refn_to_model:
                    continue
                des_pos.append((cid, refn_to_model[refn]))
            for refn in _ecds:
                if refn in refn_to_model:
                    ecd_pos.append((cid, refn_to_model[refn]))
            for refn in pocket:
                if refn in refn_to_model:
                    pocket_pos.append((cid, refn_to_model[refn]))
        plddt_vals = list(plddt.values()) if plddt else []
        ref = ref_seqs[chain_ids[0]]
        des = rec['chains_full'][chain_ids[0]]
        row = rec.copy()
        row.update({
            'tm_score': float(np.mean(tm_total)) if tm_total else float('nan'),
            'ecd_rmsd': float(np.mean(ecd_vals)) if ecd_vals else float('nan'),
            'plddt': float(np.mean(plddt_vals)) if plddt_vals else float('nan'),
            'surf_hydro': surface_hydrophobic_fraction(all_sasa, aa_map),
            'surf_hydro_des': surface_hydrophobic_fraction_at(
                all_sasa, aa_map, des_pos),
            'plddt_ecd': mean_plddt_at(plddt, ecd_pos),
            'pocket_rmsd': float(np.mean(pocket_rmsd_vals))
                if pocket_rmsd_vals else float('nan'),
            'exposed_arom': exposed_aromatics_at(all_sasa, aa_map, des_pos),
            'cys_des': cys_at(aa_map, des_pos),
            'mutations': len(mutation_calls(ref, des)),
            'model': str(model_file),
        })
        rows.append(row)
    return rows


def _aa1(resname):
    from .pdbio import AMINO3
    return AMINO3.get(resname, 'X')


def ca_coords_from_pdb(pdb_path, chain_id):
    coords, refs = [], []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                if line[21] == chain_id:
                    coords.append([float(line[30:38]), float(line[38:46]),
                                   float(line[46:54])])
                    refs.append(int(line[22:26]))
    return np.asarray(coords), refs