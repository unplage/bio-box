#!/usr/bin/env bash
# 场景 B：只有 Fab 结构和单独的靶点结构（无复合物）
# Fab 结构链：H（重链）、L（轻链）；靶点结构链：T
set -euo pipefail
cd "$(dirname "$0")/.."

: "${FAB:?请设置 FAB=/path/to/fab.pdb}"
: "${TARGET:?请设置 TARGET=/path/to/target.pdb}"

python run.py \
  --mode fab \
  --antibody-structure "$FAB" \
  --antibody-chains H,L \
  --target-structure "$TARGET" \
  --target-chains T \
  --target-name PD-L1 \
  --cdrs HCDR1,HCDR2,HCDR3,LCDR1,LCDR2,LCDR3 \
  --num-designs 8 \
  --offset 60 \
  --out-dir output/fab_pdl1
