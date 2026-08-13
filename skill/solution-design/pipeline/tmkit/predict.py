import subprocess
from pathlib import Path

import numpy as np
from .pdbio import write_fasta, residue_sasa_from_file, ca_coords
from .metrics import tm_score, parse_plddt, \
    surface_hydrophobic_fraction, mutation_calls


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
                    ecd_lookup, ref_seqs, pred_dir):
    rows = []
    _ecds = [r for a, b in ecd_lookup for r in range(a, b + 1)]
    for rec in design_records:
        name = rec['name'].replace(':', '_')
        model_file = model_paths.get(name)
        if model_file is None:
            row = rec.copy()
            row.update({'tm_score': float('nan'), 'plddt': float('nan'),
                        'ecd_rmsd': float('nan'), 'surf_hydro': float('nan'),
                        'mutations': float('nan'), 'model': ''})
            rows.append(row)
            continue
        mchains = chain_ids
        tm_total = []
        ecd_vals = []
        plddt_vals = []
        all_sasa = residue_sasa_from_file(model_file)
        aa_map = {}
        with open(model_file) as f:
            for line in f:
                if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                    ch = line[21]
                    res = int(line[22:26])
                    resname = line[17:20].strip()
                    aa_map[(ch, res)] = _aa1(resname)
        for idx, cid in enumerate(mchains):
            coords_m, _ = ca_coords_from_pdb(model_file, cid)
            coords_t, refs_t_raw = ca_coords(target, cid)
            refs_t = [target.reference_map.get(cid, {}).get(r, r)
                      for r in refs_t_raw]
            if len(coords_m) and len(coords_t):
                tm_total.append(tm_score(coords_m, coords_t))
            ecd_vals.append(_ecd_rmsd_for(model_file, cid, coords_t,
                                          refs_t, _ecds))
        plddt = parse_plddt(model_file)
        if plddt:
            plddt_vals = list(plddt.values())
        ref = ref_seqs[chain_ids[0]]
        des = rec['chains_full'][chain_ids[0]]
        row = rec.copy()
        row.update({
            'tm_score': float(np.mean(tm_total)) if tm_total else float('nan'),
            'ecd_rmsd': float(np.mean(ecd_vals)) if ecd_vals else float('nan'),
            'plddt': float(np.mean(plddt_vals)) if plddt_vals else float('nan'),
            'surf_hydro': surface_hydrophobic_fraction(all_sasa, aa_map),
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


def _ecd_rmsd_for(model_file, cid, t_coords, t_refs, ecd_refs):
    m_coords, m_refs = ca_coords_from_pdb(model_file, cid)
    if len(m_coords) == 0 or len(t_coords) == 0:
        return float('nan')
    m_map = {r: c for r, c in zip(m_refs, m_coords) if r in ecd_refs}
    t_map = {r: c for r, c in zip(t_refs, t_coords) if r in ecd_refs}
    common = sorted(set(m_map) & set(t_map))
    if not common:
        return float('nan')
    A = np.asarray([m_map[r] for r in common])
    B = np.asarray([t_map[r] for r in common])
    from .metrics import rmsd
    return rmsd(A, B)