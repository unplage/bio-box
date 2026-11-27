# Soluble membrane-protein redesign pipeline

Generic soluble redesign of membrane proteins: mutate only transferable,
solvent-exposed residues to polar/charged amino acids while preserving the
extracellular region (ECD), TM structure, pockets and oligomer interfaces,
then validate designs in silico.

The tool is **not STEAP1-specific**. A target is described by six design axes:
membrane architecture (`alpha_helical_multi`, `single_pass`, `beta_barrel`,
`tail_anchored`, `peripheral`), topology, oligomer/symmetry, functional
constraints, product format and design mode. STEAP1 (8UCD) remains a fully
annotated reference config in `config.steap1.yaml`; start new targets from
`config.template.yaml`.

Pipeline (single config drives all steps):

```
python -m tmkit.pipeline --check-deps --config config.steap1.yaml   # environment check
python -m tmkit.pipeline --config config.steap1.yaml --pdb pdb/8UCD.pdb \
    --steps annotate,mask,design,native,predict,stability,report --out output
```

Installation prerequisites and tool setup: see `INSTALL.md`.

## Recommended workflow

```bash
# 1. environment sanity check (exit code != 0 if anything is missing)
python -m tmkit.pipeline --check-deps --config config.steap1.yaml

# 2. dry-run: inspect the exact external commands before running
python -m tmkit.pipeline --config config.steap1.yaml --pdb pdb/8UCD.pdb \
    --steps mask,design,native,predict,stability,report --dry-run --out output

# 3. full run (requires ProteinMPNN + ColabFold/AlphaFold, see INSTALL.md)
python -m tmkit.pipeline --config config.steap1.yaml --pdb pdb/8UCD.pdb \
    --steps mask,design,native,predict,stability,report --out output
```

For a new target: copy `config.template.yaml`, fill the six axes, then use the
same commands with `--pdb <template.pdb>`.

No GPU available? Run the pipeline with `predict_tool: none` (or `af3` /
`external`) and point `--results-dir` at models computed elsewhere — all
evaluation/reporting still runs locally.

`strict_predict: true` (default) makes a non-zero predictor exit fatal instead
of silently producing NaN metrics and a 0/N report. `require_predictions: true`
additionally fails when any rank-1 model is missing. ColabFold receives a
directory for multi-design jobs; AlphaFold2 receives space-separated
`--fasta_paths`; Boltz-2 is driven from `predict_in/boltz/*.yaml`.

```bash
python -m tmkit.pipeline --config config.steap1.yaml --pdb pdb/8UCD.pdb \
  --steps design,native,predict,stability,report --out output \
  --proteinmpnn-dir <MPNN_DIR> --proteinmpnn-python <MPNN_PYTHON> \
  --predict-tool boltz --boltz-executable <BOLTZ_BIN> --boltz-devices 2 \
  --boltz-use-msa-server --boltz-no-kernels --require-predictions
```

`--boltz-use-msa-server` is required when the YAML inputs contain no embedded
MSA; `--boltz-no-kernels` avoids the optional CUDA kernel dependency.

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
| `annotate` | chain harmonization + reference numbering; resolves membrane class, TM/ECD/ICD/LOOP topology and warns on ambiguous single-pass/beta-barrel/assembly cases |
| `mask` | designability mask from topology + constraints: hard-fixed ECD, fixed residues/pads, pockets, disulfides, PTMs, epitopes, pore, interface mode, optional membrane slab; writes `<name>_designable.pdb` (B-factor 1 = designable) |
| `design` | ProteinMPNN on the template or on `backbone_dir`; `design_mode: scaffold/rf_diffusion` can run RFdiffusion first, then sequence its backbones; expands `model_name`/`backbone_noise`/`seeds` sweeps into independent runs |
| `predict` | assembles multimer FASTA, runs ColabFold/AlphaFold/Boltz-2 (or reads `results_dir` models); predictor failures are fatal by default and produced models are located/evaluated automatically |
| `stability` | stability gates: clamped solubility propensity, residual hydrophobicity, over-polarization, interface-mutation check; net charge, pI, helical hydrophobic moment and largest designed hydrophobic patch |
| `report` | CES evaluation (TM-score, ECD RMSD, pLDDT, ECD-only pLDDT, whole-TMD + redesigned-surface hydrophobicity, pocket RMSD, ipTM/pTM, chain-symmetry RMSD, interface-contact recovery, exposed aromatics, Cys audit, mutation count), CSV/Markdown/JSON + top variants FASTA + failed-gate summary; designs with no rank-1 model are `status=no_model` and are not counted as failed gates |
| `native` | baseline row for the unmutated target (control metrics; predicted too when `predict_native: true`) |

`--dry-run` prints the exact external commands instead of executing them.

Tool locations can be supplied per invocation without editing the YAML: `--proteinmpnn-dir`, `--proteinmpnn-python`, `--design-mode`, `--rfdiffusion-script`, `--rfdiffusion-dir`, `--predict-tool`, `--predict-executable`, `--colabfold-executable`, `--boltz-executable`, `--boltz-devices`, `--boltz-use-msa-server`, `--boltz-no-kernels`, `--require-predictions`, `--alphafold-script`. These overrides are never written back to config files.

## Config (see `config.template.yaml`, `config.steap1.yaml`)

- `membrane_class`: `alpha_helical_multi` (default), `single_pass`,
  `beta_barrel`, `tail_anchored`, `peripheral`. Single-pass/tail-anchored
  targets require exactly one TM segment and should set `nterm`/`cterm`;
  beta-barrels should provide `pore_residues`.
- `reference_sequences`: full-length sequence per chain; all segment and
  constraint positions use this 1-based reference numbering.
- `tm_segments` / `ecd_segments` / `icd_segments`: reference-numbered region
  bounds. If loops are omitted, they alternate from `first_loop`.
- `first_loop`: side of the first inter-TM loop (`ECD` default). Type-I GPCRs
  (large extracellular N-terminus) use `first_loop: ICD` plus `nterm: ECD`,
  `cterm: ICD`; `nterm`/`cterm` are also authoritative for single-pass targets.
- `signal_peptide`: optional range; SIGNAL residues are always fixed.
- `design_scope`: `tm_surface` (default), `tm_core`, `tm_plus_icd`,
  `all_non_ecd`, `all`, or an explicit region list (`[TMD, ICD]`). The legacy
  `include` / `design_tm_surface_only` keys remain supported.
- `fix_ecd` (default true) is a hard guard: even `design_scope: all` cannot
  redesign ECD residues. `fixed_regions` adds ECD/ICD/TMD/LOOP hard fixes.
- `fixed_residues` + `protect_pad`, `tm_junction_keep`: hard constraints
  (e.g. heme-coordinating H175/H268 for STEAP1).
- `pockets`: list of `{name, residues, shell, fix_shell}`. Seed + CA shell are
  fixed; the seed/shell is also used by `pocket_rmsd`. Legacy
  `pocket_residues`/`pocket_shell_dist` are supported.
- `disulfides: auto` fixes SG-SG bonded cysteines only; `ptms: auto` fixes
  N-X-S/T sequons (off by default); `epitopes`, `pore_residues` and
  `pore_shell` add further hard constraints.
- `interface_mode`: `fixed` (default, all cross-chain contacts fixed), `core`
  (only residues with >= `interface_core_contacts` partners fixed) or `none`.
  `symmetric_mask` forces an identical designable set across homooligomer chains.
- `symmetry.groups`: explicit chain-equivalence groups for MPNN ties and the
  symmetry RMSD gate. Omit for heteromers; identical reference sequences are
  grouped automatically only when no explicit groups are given.
- `interface_pairs`: optional list of `[chain_a, chain_b]` pairs; per-pair
  `interface_recovery_<a>_<b>` metrics are reported and the global gate still
  applies. Without it, all chain pairs are scored.
- `cofactors`: optional residue-name strings (e.g. `HEM`) or mappings with
  `resnames`, `chain`, `residues`, `shell` and `fix_shell`; contacting protein
  residues are hard-fixed and `ligand_rmsd` is reported for matching hetero
  residues in the predicted model.
- `surface_sasa_threshold`: only residues with SASA above threshold are
  redesigned when the design scope is surface-only. With `relative_sasa: true`,
  `surface_rsa_threshold` (SASA divided by the residue-type maximum) is used
  instead, which transfers better across protein sizes and beta-barrels.
- `membrane`: optional lipid-facing slab. `normal: auto` uses a PCA heuristic;
  prefer an explicit normal/center or `membrane_opm_pdb` (OPM/PPM, z=0).
  `exclude_aromatic_belt` keeps interfacial Trp/Phe/Tyr native. Without a slab
  the mask falls back to SASA and the pipeline prints a warning.
- `assembly`: `as_is` (default), `prebuilt` + `assembly_pdb`, or `biological`
  (warns: automatic BIOMT expansion is intentionally not guessed).
- `design_mode`: `sequence` (template), `scaffold` (`backbone_dir`),
  `rf_diffusion`/`partial_diffusion` (run RFdiffusion via
  `rfdiffusion_script` + `rfdiffusion_contigs`, then ProteinMPNN).
- `interface_recovery_dist`: C-alpha cutoff (default 8 Å) for the
  `interface_recovery` gate. Dropped automatically for monomers.
- `filters`: per-design gates. Structure/geometry (TM-score, pLDDT, ECD RMSD,
  whole-TMD and redesigned-surface hydrophobic fraction, pocket RMSD, exposed
  aromatics, Cys, mutation budget), assembly gates (`iptm`, `chain_symmetry_rmsd`,
  `interface_mut`) plus stability gates (`soluble_score`, `residual_hydro`).
  A design passes only if every gate holds; failed gates are recorded per design
  in `fail_gates` and summarized in `report.md`.
  Geometry gates also include `ecd_centroid_shift` and `ecd_tm_angle`, which
  catch rigid-body ECD displacement that an ECD-only superposition hides;
  `ecd_internal_rmsd` is reported separately as a fold-quality metric.
  Stability gates include `abs_net_charge`, `isoelectric_point`,
  `over_polar`, `max_charge_run` and `max_polar_window`, so an over-polarized
  high-|charge| sequence cannot pass on solubility score alone. A condition
  may be `['>=', value]`, `[between, lo, hi]`, a list of conditions (all must
  hold), a key template such as `mutations_{chain}`, or a key list with
  `mode: any`.
  **`iptm`, `chain_symmetry_rmsd` and `interface_mut` are dropped automatically
  for monomers and for targets with no symmetry group; ECD gates are dropped
  when no ECD is defined**, so the same template works for monomers, heteromers
  and proteins without an ectodomain.
- `template_pdb`: template structure path (PDB **or mmCIF**); lets each config be
  self-contained for `--config-dir` batch runs (no `--pdb` needed).
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
- `predict_tool`: `colabfold` (default), `alphafold`, `boltz`, or an external
  marker (`af3`, `chai`, `external`, `none`). For external predictors, point
  `--results-dir` (or `results_dir`) at the returned PDB/mmCIF models.
- `strict_predict: true` (default) fails on a non-zero predictor exit;
  `require_predictions: true` also fails if any rank-1 model is missing.
- `colabfold_amber: true` appends `--amber` (relax predicted models).
- `predict_native: true` also predicts the unmutated sequence: its metrics are
  the achievable baseline for calibrating the gates (native ECD RMSD /
  surface hydrophobicity).

## Outputs (`<out>/outputs/`)

| file | produced by | content |
|------|-------------|---------|
| `<name>_designable.pdb` | `mask` / `report` | template with B-factor 1 on designable positions (visualize mask in PyMOL) |
| `design_inputs/<tag>/<name>*.jsonl` | `design` | ProteinMPNN fixed/ties/design-chains inputs (`tag` = template or backbone) |
| `designs/<tag>/runNN/seqs/<name>.fa` | ProteinMPNN | sampled sequences (one `runNN` per sweep combination) |
| `predict_in/<name>.fa` | `predict` | per-design multimer FASTA (one record per chain) |
| `predict/` | ColabFold/AlphaFold | predicted models (`*_rank_001_*.pdb`) |
| `report.csv` | `report` | per-design metrics + pass/fail + `fail_gates` |
| `report.md` | `report` | human-readable pass list + failed-gate summary when nothing passes |
| `summary.json` | `report` | schema version, counts, filters, top variants |
| `variants_top.fa`, `variants_all_passing.fa` | `report` | passing design sequences for downstream synthesis |

## Requirements

`pip install -r requirements.txt` (biopython, numpy, freesasa, PyYAML).
ProteinMPNN and ColabFold/AlphaFold are executed as external tools.
`python -m tmkit.pipeline --check-deps` reports exactly which libraries, paths
and executables are missing and returns a non-zero exit code when any required
item fails.

## Tests

`python -m pytest tests` (85 tests: generic topology, config validation,
constraints/pockets/interface modes, masks, membrane slab, reference numbering,
design inputs, MPNN command flags and sweeps, predictor command construction,
alignment-driven model evaluation, interface recovery, stability heuristics,
gates, external formats).