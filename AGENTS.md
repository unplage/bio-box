# AGENTS.md

## Repo shape

Not an app monorepo — a **biomedical tools collection** + one Python pipeline + one OpenCode skill:

| Path | What it is |
|:---|:---|
| `tools/*/` | Standalone tools; almost all **pure-frontend single HTML** (open in browser, no build/CDN) + a few tkinter/PyQt5 scripts. Each tool dir has its own `README.md` (convention: keep it updated when you add/changed the tool). |
| `clear.html` | Deliberately stays at **repo root** (not under `tools/`). |
| `skill/solution-design/pipeline/` | Only real package with tests: `tmkit` (soluble TM-protein redesign, STEAP1 reference target). |
| `skill/target_skill/` | OpenCode skill (`SKILL.md`) for patent research — not application code. |
| `database/igblast/` | IgBLAST data files (moved here from repo root). |
| Root `README.md` | Tool index with paths — **must stay in sync** with `tools/…` locations after the reorg. |

No CI, no linter/formatter config, no `package.json`. `.gitignore`: `.venv/`, `__pycache__/`, `*.pyc`.

## Pipeline (`skill/solution-design/pipeline/`)

```bash
cd skill/solution-design/pipeline
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt   # biopython numpy freesasa PyYAML (+ system lib freesasa)

# verify env (non-zero exit if missing)
python -m tmkit.pipeline --check-deps --config config.steap1.yaml

# dry-run external cmds (works without ProteinMPNN/ColabFold)
python -m tmkit.pipeline --config config.steap1.yaml --pdb pdb/8UCD.pdb \
  --steps annotate,mask,design,native,predict,stability,report --dry-run --out output

# tests (run from inside pipeline/ — required)
python -m pytest tests
```

**Gotchas**
- Must invoke as `python -m tmkit.pipeline` **from `pipeline/`** (relative imports; do not `python tmkit/pipeline.py`).
- Tests load `config.steap1.yaml` + `pdb/8UCD.pdb` via `tests/conftest.py` (session fixtures). Single test: `python -m pytest tests/test_mask.py::test_name`.
- Full `design`/`predict` steps need external ProteinMPNN (`proteinmpnn_dir`, `soluble_model: true`) and ColabFold/AlphaFold; without GPU use `predict_tool: none` + `results_dir`.
- `pipeline/.gitignore` ignores `output/`; root also ignores `.venv/`.
- Optional `scipy` only speeds interface contact detection; absence is silent fallback, not an error.

## HTML tools

- Single-file, offline, no npm/CDN. Open with `open tools/<name>/<file>.html`.
- **`tools/msa/`**: algorithm lives twice — `msa-worker.js` (source) and a full copy inside `msa.html` (`<script id="msaWorkerSource">`). The browser only runs the HTML copy; after editing the worker, re-inject it into the HTML or changes won’t take effect.
- Export conventions users rely on: FASTA sequences are **one line** (no 60-col wrap); MSA image export includes the position ruler (ticks every column, numeric labels every 10).
- Python GUI tools (`tools/mab-db`, `tools/Tm-calculate`, `tools/target-info/*.py`) need tkinter/PyQt5 + pandas/httpx as documented in each dir README — not installed by any repo-level setup.

## Docs / language

- User-facing docs and tool READMEs are **Chinese**; keep new tool docs in the same style (table of files, 使用方法, 导出说明).
- When moving/renaming files, update root `README.md` tool table and any relative links (e.g. `tools/target-info/target.html` hrefs).

## Git

- Do not commit secrets (PATs have been pasted in chat before — rotate, never write into files/remotes you commit).
