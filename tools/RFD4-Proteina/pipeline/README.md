# RFD4-Proteina 抗体亲和力成熟 Pipeline

用 NVIDIA 托管的 **RFD4-Proteina NIM** 做 VHH / Fab 的 CDR 重设计：

- 自动识别 CDR（ANARCI，IMGT/Kabat/Chothia/Martin）；
- 根据「有没有抗体–靶点复合物」自动选择条件策略；
- 生成服务端要求的 **annotated mmCIF**（PDB 会被 API 拒绝）；
- 提交 `antibody_design` 任务，下载结构；
- 输出中文 Markdown/HTML 报告、CSV/JSON 汇总和湿实验建议。

> 这是计算候选生成工具，**不能替代 SPR/BLI、稳定性与功能实验**。

## 目录

| 文件/目录 | 说明 |
|:---|:---|
| `run.py` | 命令行入口（等价于 `python -m rfd4affinity`） |
| `rfd4affinity/` | Python 包 |
| `examples/` | VHH/Fab 示例序列、表位写法、运行脚本 |
| `tests/` | pytest 测试（CDR、条件结构、报告） |
| `requirements.txt` | 唯一可选依赖 ANARCI（需系统 HMMER） |

## 工作流

```
序列/结构 ──► 链识别与重命名 ──► CDR 识别 ──► 条件策略
                                              │
                    ┌─────────────────────────┴─────────────────────────┐
                    │ A 复合物：靶点全原子固定 + 框架主链固定 + CDR 自由  │
                    │ B 有靶点无复合物：靶点全原子固定 + 框架组 1        │
                    │ C 无靶点：仅框架内部坐标组 1                        │
                    └─────────────────────────┬─────────────────────────┘
                                              ▼
                         annotated mmCIF ──► NIM ──► designs
                                              ▼
                    report.md / report.html / designs.csv / structures/
```

## 安装

```bash
cd tools/RFD4-Proteina/pipeline
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt     # 可选但强烈推荐：ANARCI
```

- Python 3.9+；核心代码只用标准库。
- ANARCI 需要系统里有 `hmmscan`（Ubuntu: `sudo apt install hmmer`）。
- 没有 ANARCI 时会退化为保守的基序启发式，并在日志/报告中标注 `heuristic`；**生产前必须人工核对 CDR 边界**。
- 不需要 GPU：所有计算都在 NVIDIA 托管 API 上完成。

## API Key

在 <https://build.nvidia.com> 申请试用 key（`nvapi-` 开头），然后：

```bash
export NVIDIA_API_KEY="nvapi-..."
```

也可以每次 `--api-key` 传入。**key 不会写入 `request.json`/`run.log`**。

托管端点是免费试用（约 40 RPM，官方声明不适合生产）。批量或生产请用自托管
NIM：`nvcr.io/nim/ipd/rfd4-proteina:1.0.0`。

## 快速开始

### 场景 A：已有抗体–靶点复合物（最可靠）

```bash
python run.py \
  --mode vhh \
  --complex-structure complex.pdb \
  --antibody-chains A \
  --target-chains T \
  --target-name HER2 \
  --epitope T:557-573,T:591-605 \
  --cdrs all \
  --num-designs 8 \
  --out-dir output/vhh_her2
```

复合物里的天然相对朝向被保留：靶点**全原子固定**，抗体框架主链固定、侧链自由重排；
选中的 CDR 序列与局部构象同时重设计。

### 场景 B：只有 VHH 序列 + 靶点结构（无复合物）

```bash
python run.py \
  --mode vhh \
  --vhh-seq examples/vhh_sequence.fasta \
  --target-structure 8UCD \        # 也可以是本地 .pdb/.cif
  --target-chains A \
  --target-name STEAP1 \
  --cdrs CDR3 \
  --num-designs 8 \
  --out-dir output/vhh_steap1
```

- 没有抗体结构时，`--fold auto` 会用 **ESMFold NIM** 折叠 VHH（需要 API key；该端点已标记 deprecated，失败时请改传 `--antibody-structure`）。
- 靶点**全原子**绝对固定；抗体框架放入 **内部坐标组 1**（保持框架内部几何），抗体–靶点
  相对朝向由模型采样。`--offset 60` 只控制初始质心距离。
- `--epitope` 里的链号是**输入靶点结构的原始链 ID**（例如 8UCD 的 `A`）；输出 CIF 会把靶点统一重命名为 `T/U/...`。
- **必须先做表位/结合模式验证**，不能把预测 pose 当真实结构。

### 场景 C：只有抗体、没有靶点（非靶向）

```bash
python run.py \
  --mode vhh \
  --vhh-seq examples/vhh_sequence.fasta \
  --cdrs all \
  --num-designs 8
```

只优化 CDR 序列/构象，结果只反映可开发性，不能解读为亲和力提升；报告会明确警告。

### Fab

Fab 需要 VH/VL 结构（ESMFold 不能可靠组装 VH–VL）：

```bash
python run.py \
  --mode fab \
  --antibody-structure fab.pdb \
  --antibody-chains H,L \
  --target-structure target.pdb \
  --target-chains T \
  --cdrs HCDR3,LCDR1,LCDR2,LCDR3 \
  --num-designs 8
```

也可以 `--complex-structure fab_complex.pdb --antibody-chains H,L --target-chains T`。

## 常用参数

| 参数 | 说明 |
|:---|:---|
| `--mode vhh\|fab` | 抗体形式，默认 `vhh` |
| `--vhh-seq / --vh-seq / --vl-seq` | 序列字符串或 FASTA 文件 |
| `--complex-structure` | 抗体–靶点复合物；与 `--target-structure` 互斥 |
| `--antibody-structure` | 仅抗体结构 |
| `--antibody-chains` | 抗体链；Fab 为 `H,L`（重,轻） |
| `--target-structure` | 靶点路径或 4 位 PDB ID（自动从 RCSB 下载） |
| `--target-chains` | 靶点链，默认全部 |
| `--cdrs` | `all`、`CDR3`、`HCDR1,LCDR3` 等；Fab 单独写 `CDR3` 会报歧义 |
| `--cdr-scheme` | `imgt`（默认）/`kabat`/`chothia`/`martin` |
| `--cdr-ranges` | 手动序列坐标，如 `H:99-109,L:89-97` |
| `--epitope` | 靶点热点，如 `557-573,591-605` 或 `T:557-573`（链号用**原始靶点链 ID**，不是输出中重命名后的 T/U） |
| `--fold` | `auto`（默认）/`esmfold`/`none` |
| `--offset` | docking 模式下靶点初始质心距离，默认 60 Å |
| `--fix-sidechains` | 抗体侧链也固定；靶点始终全原子固定，默认抗体只固定主链 |
| `--num-designs` | 设计数，默认 4 |
| `--num-steps` | 扩散步数；默认用 `--recycle` 预设（60 步 / 185 trunk pass） |
| `--recycle` | `thinking_low/medium/high` |
| `--seed` | 随机种子 |
| `--metrics-tags` | 逗号分隔，覆盖默认 `design,antibody,developability,...` |
| `--dry-run` | 只构建 `input_complex.cif` + `request.json`，不调用 API |
| `--replay response.json` | 复用已保存响应，重新出报告 |
| `--no-html` | 不生成 HTML |

`--help` 有完整列表和示例。

## CDR 识别

- ANARCI 把序列按选定 scheme 编号，再把 scheme 位置映射回 **输入序列的 1-based 位置**。
- CDR 区间：IMGT `27-38 / 56-65 / 105-117`；Kabat `31-35 / 50-65 / 95-102`；
  Chothia `26-32 / 52-56 / 95-102`；Martin `26-35 / 50-58 / 95-102`。
- 若输入的是结构 + 单独序列，序列与结构不一致时会做对齐并把 CDR 位置投影到
  结构的 `auth_seq_id`；一致性低于 80% 直接报错。
- `--cdr-ranges` 的坐标始终是**输入序列坐标**，优先级最高。
- `--dry-run` 后检查 `input_complex.cif` 的 `mask_sequence_1_residue` 列，
  确认 CDR 残基是 `False`（自由设计），框架是 `True`。

## 条件结构（annotated mmCIF）

NIM 只接受 UTF-8 mmCIF；本流程会把 PDB/普通 CIF 清洗后重写，并加上这些列：

| 列 | 含义 |
|:---|:---|
| `mask_coordinate_1_atom` | True = C_CRD 绝对坐标固定 |
| `condition_CoordinateGroup_1_atom` | C_CRDGRP：-1 未知、0 绝对、≥1 内部坐标组 |
| `mask_InternalCoordinate_1_atom` | True = 为组 ≥1 准备内部坐标 |
| `mask_sequence_1_residue` | True = C_SEQ 固定序列；False = 设计 |
| `mask_hotspot_1_atom` | True = C_HOT（只标侧链原子） |
| `condition_index-group_1_token` | 链/token 级坐标组编号 |
| `mask_index_1_residue` | True = 参与索引 |

默认策略：

- 靶点**全原子**绝对固定；抗体框架主链绝对/内部固定，**侧链自由重排**
  （`--fix-sidechains` 可把抗体侧链也固定）；
- 选中 CDR：不加坐标/序列约束 → 序列与局部构象一起重设计；
- 水、非标准残基、配体、多模型、altloc 会被丢弃；MSE→MET；
- `label_seq_id` 保留作者编号（插入码冲突时退化为 1..N），`auth_seq_id` 始终是原编号；
- 单请求上限 1024 tokens、32 MB；超限会给出明确错误。

> 也支持官方 `condition_spec` JSON，但本流程统一用 `annotated_cif`：
> 它自包含、可审计，且与官方 antibody 示例的列定义一致。

## 输出

```
output/<时间>_<靶点>_<模式>/
├── input_complex.cif      # 带条件的输入（可直接拿去调 API）
├── request.json           # 发给 NIM 的请求（不含 key）
├── response.json          # 原始响应
├── designs.csv            # 排名、CDR/框架突变、pI、GRAVY、风险位点
├── designs.json           # 完整解析结果
├── structures/            # 每个 design 的 mmCIF
├── report.md              # 中文报告
├── report.html            # 自包含单文件报告
└── run.log                # 运行日志（不含 key）
```

报告包含：运行摘要、输入链/CDR/热点、条件策略、候选总表、逐候选序列与突变、
pLDDT、pI/GRAVY/净电荷/风险基序、NIM 指标、聚合指标、湿实验建议、局限与文件清单。

**排名是启发式的**：通过标准 → pLDDT → 风险位点数，不等于真实亲和力排序。

## 湿实验建议（生成在报告第 6 节）

- VHH：pET-21b/22b 周质表达 + Ni-NTA + SEC；Fab：HEK293 共转染 + Protein L/G + SEC。
- 先 ELISA/流式初筛，再 BLI/SPR 测 kon/koff/KD；母本、无关抗体、缓冲液三组对照。
- 竞争实验/HDX-MS 验证表位；关键突变做丙氨酸扫描或回复突变。
- DSF Tm、SEC 单体率、冻融/40 °C 加速、AC-SINS/PSR/HIC。
- 组合突变时避免位点重叠与疏水斑块叠加，保留一个“只改 CDR3”的极简版本。

## 局限

- 计算预测 ≠ 实测亲和力；必须做 SPR/BLI 和功能实验。
- 场景 B/C 的结合模式/靶点缺失，不确定性最大；报告会显式标注。
- RFD4 可能改变 CDR 长度或框架边界；报告统计框架突变，优先推进框架突变少的候选。
- 托管 NIM 是试用端点，可能限流/变更；生产请自托管。
- 只处理 20 种标准氨基酸；抗体标签、His-tag、非天然残基请先去除。

## 测试

```bash
cd tools/RFD4-Proteina/pipeline
python -m pytest tests
```

测试不需要 API key，也不会联网；使用合成 PDB 覆盖三种场景、CDR 选择、
annotated CIF 列和报告生成。

## 常见问题

| 现象 | 处理 |
|:---|:---|
| `HTTP 401/403` | 检查 `NVIDIA_API_KEY`，试用 key 以 `nvapi-` 开头 |
| `HTTP 429` | 试用端点限流；降低 `--num-designs` 或稍后重试（客户端会自动退避） |
| `No standard protein atoms` | 输入可能是核酸/配体，或残基名不在标准表 |
| `ESMFold could not fold` | 端点已 deprecated；改用 `--antibody-structure` |
| `2048/1024 tokens` 超限 | 只保留靶点胞外域/单结构域，别把整个膜蛋白放进去 |
| CDR 边界可疑 | 用 `--cdr-scheme kabat` 或 `--cdr-ranges` 手动指定，并做 `--dry-run` 核对 |
| Fab 没有结构 | 先预测/实验解析 VH–VL；不要用单链 ESMFold 结果冒充 Fab |

## 参考

- 模型分析（本目录上级）：`../README.md`
- 官方文档：<https://rfd4-proteina.org/latest/>
- NIM 支持矩阵：<https://docs.nvidia.com/nim/bionemo/rfd4-proteina/latest/support-matrix.html>
- 论文：`../RFD4-Proteina-paper.pdf`
