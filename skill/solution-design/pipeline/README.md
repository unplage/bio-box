# Soluble TM Redesign pipeline

Soluble redesign of transmembrane targets for expression/detergent approaches:
mutate only the solvent-exposed, transferable TM-facing residues to polar/charged
amino acids while preserving TM structure, ECD geometry and interaction surfaces,
then validate designs in silico by structure prediction.

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
| `predict` | assembles multimer FASTA, runs ColabFold/AlphaFold (or reads `results_dir` models) |
| `stability` | stability gates: clamped solubility propensity, residual hydrophobicity, over-polarization, interface-mutation check |
| `report` | CES evaluation (TM-score, ECD RMSD, pLDDT, surface-hydrophobicity, mutation count), CSV/Markdown/JSON + top variants FASTA |
| `native` | baseline row for the unmutated target (control metrics) |

`--dry-run` prints the exact external commands instead of executing them.

## Config (see `config.example.yaml`)

- `reference_sequences`: full-length sequence per chain; 1-based UCD positional numbering.
- `tm_segments` / `ecd_segments`: reference-numbered region bounds.
- `include`: `designable` (TMD only, default), `all`, `tm_plus_icd`.
- `fixed_residues` + `protect_pad`, `tm_junction_keep`: hard constraints.
- `fix_interface_contacts` (default on) hard-fixes residues contacting another
  chain (`interface_contact_dist`, 6 Å) to protect oligomer interfaces;
  `symmetric_mask` forces an identical designable set across chains for
  homooligomers.
- `surface_sasa_threshold`, `design_tm_surface_only`: only TM residues with
  solvent-accessible surface area above threshold are redesigned.
- `filters`: per-design gates. Structure/geometry (TM-score, pLDDT, ECD RMSD,
  surface hydrophobic fraction, mutation budget) plus stability gates
  (`soluble_score`, `residual_hydro`, `interface_mut`). A design passes only if
  every gate holds.
- `proteinmpnn_dir`: folder containing `protein_mpnn_run.py`
  (e.g. /opt/proteinmpnn). `soluble_model: true` uses the soluble weights.
- `predict_tool`: `colabfold` (default), `alphafold`, or `none`.
  `results_dir` may point at precomputed `*_rank_001_*.pdb` models.

## Outputs (`<out>/outputs/`)

| file | produced by | content |
|------|-------------|---------|
| `<name>_designable.pdb` | `mask` / `report` | template with B-factor 1 on designable positions (visualize mask in PyMOL) |
| `design_inputs/<name>*.jsonl` | `design` | ProteinMPNN fixed/ties/design-chains inputs |
| `designs/seqs/<name>.fa` | ProteinMPNN | sampled sequences |
| `predict_in/predict_input.fasta` | `predict` | multimer FASTA (chains joined by `:`) |
| `predict/` | ColabFold/AlphaFold | predicted models (`*_rank_001_*.pdb`) |
| `report.csv` | `report` | per-design metrics + pass/fail |
| `report.md` | `report` | human-readable pass list |
| `summary.json` | `report` | counts, filters, top variants |
| `variants_top.fa`, `variants_all_passing.fa` | `report` | passing design sequences for downstream synthesis |

## Requirements

`pip install -r requirements.txt` (biopython, numpy, freesasa, PyYAML).
ProteinMPNN and ColabFold/AlphaFold are executed as external tools.
`python -m tmkit.pipeline --check-deps` reports exactly which libraries, paths
and executables are missing and returns a non-zero exit code when any required
item fails.