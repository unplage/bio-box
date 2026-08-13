# Soluble TM Redesign pipeline

Soluble redesign of transmembrane targets for expression/detergent approaches:
mutate only the solvent-exposed, transferable TM-facing residues to polar/charged
amino acids while preserving TM structure, ECD geometry and interaction surfaces,
then validate designs in silico by structure prediction.

STEAP1 (8UCD) is the bundled reference target: a complete production config is
in `config.steap1.yaml` (6 TM helices, ECD EC1-3 fixed, heme pocket H175/H268
protected, C3-symmetric mask, AMG509/trimer-aware gates).

Pipeline (single config drives all steps):

```
python -m tmkit.pipeline --check-deps --config config.example.yaml   # environment check
python -m tmkit.pipeline --config config.example.yaml --pdb pdb/8UCD.pdb \
    --steps annotate,mask,design,native,predict,stability,report --out output
```

Installation prerequisites and tool setup: see `INSTALL.md`.

## Recommended workflow

```bash
# 1. environment sanity check (exit code != 0 if anything is missing)
python -m tmkit.pipeline --check-deps --config config.example.yaml

# 2. dry-run: inspect the exact external commands before running
python -m tmkit.pipeline --config config.example.yaml --pdb pdb/8UCD.pdb \
    --steps mask,design,native,predict,stability,report --dry-run --out output

# 3. full run (requires ProteinMPNN + ColabFold/AlphaFold, see INSTALL.md)
python -m tmkit.pipeline --config config.example.yaml --pdb pdb/8UCD.pdb \
    --steps mask,design,native,predict,stability,report --out output
```

No GPU available? Run the pipeline with `predict_tool: none` and point
`results_dir` at models computed elsewhere (`*_rank_001_*.pdb` per design) —
all evaluation/reporting still runs locally.

## Steps

| step | description |
|------|-------------|
| `annotate` | chain harmonization + UCD positional numbering; dumps reference/topology JSON |
| `mask` | SASA-based designability mask of the TMD; hard-fixes oligomer interface residues; writes `<name>_designable.pdb` (B-factor 1 = designable) |
| `design` | generates ProteinMPNN inputs (fixed/ties/design-chains JSONL) and runs `protein_mpnn_run.py` |
| `predict` | assembles multimer FASTA, runs ColabFold/AlphaFold (or reads `results_dir` models); when the pipeline runs the predictor itself, produced models are located and evaluated automatically |
| `stability` | stability gates: clamped solubility propensity, residual hydrophobicity, over-polarization, interface-mutation check |
| `report` | CES evaluation (TM-score, ECD RMSD, pLDDT, ECD-only pLDDT, redesigned-surface hydrophobicity, pocket RMSD, exposed aromatics, Cys audit, mutation count), CSV/Markdown/JSON + top variants FASTA + failed-gate summary |
| `native` | baseline row for the unmutated target (control metrics; predicted too when `predict_native: true`) |

`--dry-run` prints the exact external commands instead of executing them.

## Config (see `config.steap1.yaml` / `config.example.yaml`)

- `reference_sequences`: full-length sequence per chain; 1-based UCD positional numbering.
- `tm_segments` / `ecd_segments` / `icd_segments`: reference-numbered region bounds.
- `include`: `designable` (TMD only, default), `all`, `tm_plus_icd`.
- `fixed_residues` + `protect_pad`, `tm_junction_keep`: hard constraints
  (e.g. heme-coordinating H175/H268 for STEAP1).
- `fix_interface_contacts` (default on) hard-fixes residues contacting another
  chain (`interface_contact_dist`, 6 Å) to protect oligomer interfaces;
  `symmetric_mask` forces an identical designable set across chains for
  homooligomers.
- `surface_sasa_threshold`, `design_tm_surface_only`: only TM residues with
  solvent-accessible surface area above threshold are redesigned.
- `pocket_residues`: reference positions whose geometry is tracked by the
  `pocket_rmsd` gate (heme pocket preservation).
- `filters`: per-design gates. Structure/geometry (TM-score, pLDDT, ECD RMSD,
  surface hydrophobic fraction, pocket RMSD, exposed aromatics, Cys, mutation
  budget) plus stability gates (`soluble_score`, `residual_hydro`,
  `interface_mut`). A design passes only if every gate holds; failed gates are
  recorded per design in `fail_gates` and summarized in `report.md`.
- `proteinmpnn_dir`: folder containing `protein_mpnn_run.py`
  (e.g. /opt/proteinmpnn). `soluble_model: true` uses the soluble weights.
- `num_seqs` / `temperature`: MPNN sampling (config values are used unless
  overridden by CLI flags).
- `predict_tool`: `colabfold` (default), `alphafold`, or `none`.
  `results_dir` may point at precomputed `*_rank_001_*.pdb` models.
- `colabfold_amber: true` appends `--amber` (relax predicted models).
- `predict_native: true` also predicts the unmutated sequence: its metrics are
  the achievable baseline for calibrating the gates (native ECD RMSD /
  surface hydrophobicity).

## Outputs (`<out>/outputs/`)

| file | produced by | content |
|------|-------------|---------|
| `<name>_designable.pdb` | `mask` / `report` | template with B-factor 1 on designable positions (visualize mask in PyMOL) |
| `design_inputs/<name>*.jsonl` | `design` | ProteinMPNN fixed/ties/design-chains inputs |
| `designs/seqs/<name>.fa` | ProteinMPNN | sampled sequences |
| `predict_in/predict_input.fasta` | `predict` | multimer FASTA (chains joined by `:`; native included when `predict_native: true`) |
| `predict/` | ColabFold/AlphaFold | predicted models (`*_rank_001_*.pdb`) |
| `report.csv` | `report` | per-design metrics + pass/fail + `fail_gates` |
| `report.md` | `report` | human-readable pass list + failed-gate summary when nothing passes |
| `summary.json` | `report` | counts, filters, top variants |
| `variants_top.fa`, `variants_all_passing.fa` | `report` | passing design sequences for downstream synthesis |

## Requirements

`pip install -r requirements.txt` (biopython, numpy, freesasa, PyYAML).
ProteinMPNN and ColabFold/AlphaFold are executed as external tools.
`python -m tmkit.pipeline --check-deps` reports exactly which libraries, paths
and executables are missing and returns a non-zero exit code when any required
item fails.

## Tests

`python -m pytest tests` (masks, reference numbering, design inputs, MPNN
command flags, alignment-driven model evaluation, gates).