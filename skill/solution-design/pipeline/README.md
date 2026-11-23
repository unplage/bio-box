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

No GPU available? Run the pipeline with `predict_tool: none` (or `boltz` / `af3`
/ `external`) and point `--results-dir` at models computed elsewhere — all
evaluation/reporting still runs locally.

### External structure predictors (Boltz-2 / AlphaFold3 / ColabFold)

Templates and predicted models may be **PDB or mmCIF**. Model discovery is
recursive and understands the common layouts:

| tool | example model file |
|------|--------------------|
| ColabFold / AlphaFold | `<name>_unrelaxed_rank_001_*.pdb` |
| Boltz-2 | `boltz_results_*/predictions/<name>/<name>_model_0.cif` |
| AlphaFold3 | `<name>/<name>_model.cif` |

Confidence sidecars (`confidence_*.json` for Boltz, `summary_confidences.json`
for AF3, `*scores*.json` for ColabFold) are read for `ptm` / `iptm`; pLDDT is
read from the model B-factor. Custom filename patterns can be supplied with
`results_globs` in the config (or the defaults above are used).

Use the pipeline to emit the exact inputs to predict. The `predict` step writes,
for every design (+ native):

```
outputs/predict_in/<name>.fa        # ColabFold / AF2 multimer FASTA
outputs/predict_in/boltz/<name>.yaml  # Boltz-2 input (explicit chain ids)
outputs/predict_in/af3/<name>.json    # AlphaFold3 input spec
```

```bash
# 1) emit inputs (no predictor is run in --dry-run)
python -m tmkit.pipeline --config config.yaml --pdb template.pdb \
  --steps design,native,predict --dry-run --out output
# 2) run externally, e.g.
#    boltz predict output/outputs/predict_in/boltz --use_msa_server --out_dir <dir>
# 3) evaluate returned models (PDB or mmCIF) — flat or nested layouts both work
python -m tmkit.pipeline --config config.yaml --pdb template.pdb \
  --steps native,predict,report --results-dir /path/to/predictions --out output
```

## Steps

| step | description |
|------|-------------|
| `annotate` | chain harmonization + UCD positional numbering; dumps reference/topology JSON |
| `mask` | SASA-based designability mask of the TMD; optional membrane slab (lipid-facing only, aromatic belt kept native); hard-fixes oligomer interface residues; writes `<name>_designable.pdb` (B-factor 1 = designable) |
| `design` | generates ProteinMPNN inputs (fixed/ties/design-chains JSONL); expands `model_name`/`backbone_noise`/`seeds` sweeps into independent runs and runs `protein_mpnn_run.py` |
| `predict` | assembles multimer FASTA, runs ColabFold/AlphaFold (or reads `results_dir` models); when the pipeline runs the predictor itself, produced models are located and evaluated automatically |
| `stability` | stability gates: clamped solubility propensity, residual hydrophobicity, over-polarization, interface-mutation check; net charge, pI, helical hydrophobic moment and largest designed hydrophobic patch |
| `report` | CES evaluation (TM-score, ECD RMSD, pLDDT, ECD-only pLDDT, whole-TMD + redesigned-surface hydrophobicity, pocket RMSD, ipTM/pTM, chain-symmetry RMSD, interface-contact recovery, exposed aromatics, Cys audit, mutation count), CSV/Markdown/JSON + top variants FASTA + failed-gate summary |
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
- `membrane` (optional): restrict designability to the lipid-facing slab.
  `normal: [0,0,1]` (or `normal: auto` for PCA of TM C-alphas), `thickness`
  (hydrophobic core, default 30 Å), `headgroup` (extra interface zone, default
  15 Å). `exclude_aromatic_belt: true` keeps interfacial Trp/Phe/Tyr native;
  `restrict_to_core: true` forbids headgroup design entirely. Alternatively set
  `membrane_opm_pdb` to an OPM/PPM transformed structure (membrane plane z=0).
- `interface_recovery_dist`: C-alpha cutoff (default 8 Å) for the
  `interface_recovery` gate — the fraction of template cross-chain contacts
  reproduced by a predicted model. Dropped automatically for monomers.
- `pocket_residues`: reference positions whose geometry is tracked by the
  `pocket_rmsd` gate (heme pocket preservation).
- `filters`: per-design gates. Structure/geometry (TM-score, pLDDT, ECD RMSD,
  whole-TMD and redesigned-surface hydrophobic fraction, pocket RMSD, exposed
  aromatics, Cys, mutation budget), assembly gates (`iptm`, `chain_symmetry_rmsd`,
  `interface_mut`) plus stability gates (`soluble_score`, `residual_hydro`).
  A design passes only if every gate holds; failed gates are recorded per design
  in `fail_gates` and summarized in `report.md`.
  **`iptm`, `chain_symmetry_rmsd` and `interface_mut` are dropped automatically
  for single-chain (monomer) targets**, so the same config template works for
  monomers (e.g. GPCRs) and oligomers.
- `template_pdb`: template structure path (PDB **or mmCIF**); lets each config be
  self-contained for `--config-dir` batch runs (no `--pdb` needed).
- `first_loop`: side of the first inter-TM loop when loops are inferred from
  `tm_segments`. Default `ECD` (e.g. STEAP1). **Type-I GPCRs (large extracellular N-terminus) have
  an extracellular N-terminal ectodomain and an *intracellular* TM1–TM2 loop**,
  so use `first_loop: ICD` together with `nterm: ECD`, `cterm: ICD`.
  `nterm`/`cterm` default to `ICD`.
- `results_globs`: optional list of custom `{name}` filename patterns for
  external model discovery.
- `proteinmpnn_dir`: folder containing `protein_mpnn_run.py`
  (e.g. /opt/proteinmpnn). `soluble_model: true` uses the soluble weights.
- `num_seqs` / `temperature`: MPNN sampling (config values are used unless
  overridden by CLI flags).
- `model_name` / `backbone_noise` / `seeds`: each accepts a scalar or a list;
  every combination becomes an independent MPNN run (separate
  `designs/runNN/` folder and FASTA) so noise/seed/weight sweeps do not
  overwrite each other.
- `predict_tool`: `colabfold` (default), `alphafold`, or an external marker
  (`boltz`, `af3`, `chai`, `external`, `none`). For external predictors, point
  `--results-dir` (or `results_dir`) at the returned PDB/mmCIF models.
- `colabfold_amber: true` appends `--amber` (relax predicted models).
- `predict_native: true` also predicts the unmutated sequence: its metrics are
  the achievable baseline for calibrating the gates (native ECD RMSD /
  surface hydrophobicity).

## Outputs (`<out>/outputs/`)

| file | produced by | content |
|------|-------------|---------|
| `<name>_designable.pdb` | `mask` / `report` | template with B-factor 1 on designable positions (visualize mask in PyMOL) |
| `design_inputs/<name>*.jsonl` | `design` | ProteinMPNN fixed/ties/design-chains inputs |
| `designs/runNN/seqs/<name>.fa` | ProteinMPNN | sampled sequences (one `runNN` per sweep combination) |
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

`python -m pytest tests` (masks, membrane slab, reference numbering, design
inputs, MPNN command flags and sweeps, alignment-driven model evaluation,
interface recovery, stability heuristics, gates, external formats).