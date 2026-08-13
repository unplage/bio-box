# Installation

Running order: set up the Python environment -> install ProteinMPNN ->
install a structure predictor (ColabFold recommended) -> wire config ->
`--check-deps` to verify -> run.

## 1. Python environment

Requires Python 3.8+.

```bash
cd skill/solution-design/pipeline
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt     # biopython, numpy, freesasa, PyYAML
```

Optional: `pip install scipy` for a faster inter-chain contact (KD-tree)
detection; without it the pipeline silently falls back to Biopython's bundled
KD-tree (no behaviour change).

freesasa needs its C library (`sudo apt install freesasa` or a conda install
`conda install -c bioconda freesasa`). Verify with:

```bash
python -m tmkit.pipeline --check-deps
```

## 2. ProteinMPNN (step `design`)

```bash
git clone https://github.com/dauparas/ProteinMPNN.git /opt/proteinmpnn
```

The directory must contain `protein_mpnn_run.py` and the model weight folders.
For soluble redesign the config sets `soluble_model: true`, which requires the
`soluble_model_weights` folder shipped with the repo exactly under
`/opt/proteinmpnn/soluble_model_weights`. Point the config at it:

```yaml
proteinmpnn_dir: /opt/proteinmpnn
soluble_model: true
```

The pipeline itself never calls the code directly — it only builds and runs
`python protein_mpnn_run.py ...`, so any standard ProteinMPNN checkout works.

## 3. Structure prediction (step `predict`)

Choose one strategy:

### 3a. ColabFold (recommended)

```bash
conda create -n colabfold python=3.10 -y
conda activate colabfold
conda install -c conda-forge -c bioconda colabfold -y
# offline (recommended for large-scale):
colabfold/bin/download_jackhmmer.sh .
```

`colabfold_batch` must be on `PATH`; GPU recommended (CPU with small
`--num-recycle` works for a handful of designs). First run builds the MSA —
allow network access or point at the local databases.

Config:
```yaml
predict_tool: colabfold
colabfold_model: alpha2_multimer_v3
num_recycles: 3
```

### 3b. AlphaFold (full DBs)

Install per official docs and fill every database path in config:

```yaml
predict_tool: alphafold
alphafold_script: /opt/alphafold/run_alphafold.py
data_dir: /databases
uniref90: /databases/uniref90/uniref90.fasta
mgnify:   /databases/mgnify/mgy_clusters_2022_05.fa
bfd:      /databases/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt
uniref30: /databases/uniref30/UniRef30_2021_06
pdb_seqres: /databases/pdb_seqres/pdb_seqres.txt
template_mmcif: /databases/pdb_mmcif/mmcif_files
obsolete: /databases/pdb_mmcif/obsolete.dat
```

### 3c. No predictor installed

Run the pipeline without GPU using externally computed models. Compute the
structures elsewhere (ColabFold notebook, another machine), place models in one
folder, and set:

```yaml
predict_tool: none
results_dir: /path/to/models
```

The pipeline reads `*_rank_001_*.pdb` (or `*_pred_0.pdb`) per design and runs
the full evaluation/report locally.

## 4. Config wiring

Copy `config.example.yaml` and edit:
- `reference_sequences` per chain (drives UCD positional numbering)
- `tm_segments` / `ecd_segments` (reference-numbered bounds)
- `proteinmpnn_dir`, predictor settings, `results_dir`
- optional `filters` incl. stability gates

## 5. Sanity checks

```bash
python -m tmkit.pipeline --config config.yaml --check-deps     # all paths/tools OK?
python -m tmkit.pipeline --config config.yaml --pdb pdb/8UCD.pdb \
    --steps mask,design,native,predict,stability,report --dry-run   # dry run
```

`--dry-run` prints every external command without executing it, and also works
without ProteinMPNN/ColabFold installed.

## 6. Troubleshooting

| symptom | cause / fix |
|---|---|
| `module 'freesasa' has no attribute ...` | wrong freesasa API version; `pip install --upgrade freesasa` or reinstall C lib |
| pipeline hangs on `mask` | scipy missing on huge structures; `pip install scipy` |
| `protein_mpnn_run.py` not found | `proteinmpnn_dir` wrong, or `--use_soluble_model` set without `soluble_model_weights` |
| `colabfold_batch: command not found` | not activated / not on PATH |
| empty `report.csv` after `predict` | model files missing; check naming `*_rank_001_*.pdb` and `results_dir` |
| must run as `python -m tmkit.pipeline` | package uses relative imports; run from inside `pipeline/` |
| `pass` all False | the strict multi-gate filters should be tuned per target, or design pool too small |