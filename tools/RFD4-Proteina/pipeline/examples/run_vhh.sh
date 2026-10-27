#!/usr/bin/env bash
# 场景 A：已有 VHH–靶点复合物（最可靠）
# 复合物中的链：VHH = A，靶点 = T（按实际文件修改）
set -euo pipefail
cd "$(dirname "$0")/.."

: "${COMPLEX:?请设置 COMPLEX=/path/to/complex.pdb}"

python run.py \
  --mode vhh \
  --complex-structure "$COMPLEX" \
  --antibody-chains A \
  --target-chains T \
  --target-name HER2 \
  --epitope T:557-573,T:591-605 \
  --cdrs all \
  --num-designs 8 \
  --out-dir output/vhh_her2
