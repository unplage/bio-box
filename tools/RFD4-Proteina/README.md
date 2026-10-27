# RFD4-Proteina 开源模型分析

> Baker Lab（华盛顿大学蛋白质设计研究所，IPD）× NVIDIA GenAIR 联合推出的新一代生物分子生成式设计模型。
> 2026 年 10 月 21 日开源（BSD-3-Clause），代码、权重、文档、在线 Demo 全部公开。

## 📄 文件

| 文件 | 说明 |
|:---|:---|
| `RFD4-Proteina-paper.pdf` | 原始预印本论文（54.8 MB，含正文与扩展数据） |
| `pipeline/` | 可运行的 VHH/Fab 抗体亲和力成熟 pipeline（调用 NVIDIA NIM，生成报告） |
| `README.md` | 本分析文档 |

## 📚 文献信息

| 项目 | 内容 |
|:---|:---|
| 标题 | *Biomolecular design through latent reasoning with RFD4-Proteina* |
| 作者 | Nathaniel Corley, Kieran Didi, Maxwell W. Kazman, Rafael I. Brent 等（共同一作，排名随机）；通讯作者 Kieran Didi、David Baker |
| 单位 | 华盛顿大学 IPD / 生物化学系；NVIDIA；牛津大学；清华大学；ETH Zurich；曼彻斯特大学等 |
| 年份 | 2026（预印本 Preprint） |
| 许可 | BSD 3-Clause（模型与权重），版权归 University of Washington |
| 项目页 | https://rfd4-proteina.org/ |
| 论文 PDF | https://rfd4-proteina.org/assets/paper/RFD4-Proteina.pdf |
| 代码 | https://github.com/RosettaCommons/RFD4-Proteina |
| 权重（原生） | https://huggingface.co/nvidia/RFD4-Proteina |
| 权重（Diffusers） | https://huggingface.co/nvidia/RFD4-Proteina-Diffusers |
| 文档 | https://rfd4-proteina.org/latest/ |
| PyPI | https://pypi.org/project/rfd4-proteina/ |

**BibTeX**

```bibtex
@article{corley2026rfd4proteina,
  title  = {Biomolecular design through latent reasoning with {RFD4-Proteina}},
  author = {Corley, Nathaniel and Didi, Kieran and Kazman, Maxwell W. and Brent, Rafael I. and others},
  year   = {2026},
  note   = {Preprint},
  url    = {https://rfd4-proteina.org/assets/paper/RFD4-Proteina.pdf}
}
```

## 🧠 模型核心

| 维度 | 内容 |
|:---|:---|
| 类型 | 流匹配（flow matching）生成模型，序列–结构联合生成（sequence–structure co-design） |
| 参数 | 约 4.356 亿参数（435,622,465），FP32 存储 1.74 GB；推理以 BF16 运行 |
| 核心创新 | 统一的 32 维「化学 + 构象」隐空间，覆盖标准/修饰氨基酸、糖类、DNA、RNA |
| 内在推理 | 循环精修（recurrent refinement），可按 low/medium/high「思考」预算自适应分配算力 |
| 外在推理 | 推理时搜索（reward-guided beam search），用几何/能量/结构预测 oracle 打分剪枝 |
| 物理约束 | 隐空间以力场能量监督，倾向物理合理的构象 |
| 输入 | 带注释的 CIF 结构，或 YAML/JSON 条件规格（固定 motif、hotspot、埋藏度、氢键等） |
| 输出 | CIF 格式设计结构（残基身份 + 全原子坐标） |
| 采样 | 已发布 `thinking_medium`：60 步、180 次 trunk 传递；可调 low/high |
| 加速 | 默认 `torch.compile`；可选 BioIR 后端（`inference_backend=bioir`，仅 CUDA 13） |

### 支持的 7 类设计任务

| 任务 | 作用 |
|:---|:---|
| `generation` | 无/有条件下从头生成（单体、蛋白/DNA binder） |
| `motif_scaffolding` | 固定 motif 位点，在其周围生成支架 |
| `enzyme_design` | 围绕未编号催化位点设计酶支架 |
| `small_molecule_binder` | 围绕小分子配体设计结合口袋 |
| `antibody_design` | 在固定 Ig 框架上重设计 CDR 环，可针对固定抗原；也覆盖 TCR |
| `partial_diffusion` | 部分加噪再精修，用于优化或增加多样性 |
| `inverse_folding` | 固定骨架，重设计序列与侧链 |

### 关键性能（论文 / 官网数据）

- 在蛋白结合、小分子结合、酶支架、DNA 结合四类基准上，**通过率与吞吐量同时优于所有对比方法**。
- 蛋白结合：high-thinking 单次生成 clustered pass rate **10.6%**（次优 BindCraft 5.9%），吞吐高约 **42 倍**；低→高思考，通过率 7.7%→10.6%，吞吐 2388→986 designs/GPU-hour。
- 小分子结合：14.3% vs 2.4%（6.1 倍）；酶支架 12.1% vs 5.7%；DNA 结合 3.5% vs 0.8%。
- 目标覆盖：24 个蛋白靶点解出 22 个，77 个小分子靶点全解，37 个酶任务解出 35 个。
- 无需额外训练即可胜任侧链打包、固定骨架序列设计、构象采样、设计精修等专门任务。

### 湿实验验证（论文）

- 钙依赖钙调蛋白 binder：18 个中 12 个 < 1 µM，无钙时不结合。
- 肽–MHC binder：对靶肽 >100 倍特异性，可区分仅差 1–2 个残基的脱靶肽。
- 内在无序肽 binder：YAP1 112 pM、Tau 527 pM（据称为该表位的首个从头设计 binder）。
- TrkA binder：可区分物种与旁系同源（hTrkA 188 个中 158 个亚微摩尔）。
- 糖基化屏蔽：9 个中 8 个 >90% oligomannose。
- 磷酸化特异 binder：组蛋白 H3 pS10 / pT11。
- 双金属磷酸三酯酶：催化效率较此前设计高至 1000 倍（kcat/KM 达 3600 M⁻¹s⁻¹）。
- 对映选择性连接酶：用于 lenacapavir 合成的酰胺成键步骤，最佳设计转化 79% 的 (S)-胺——**首个用于 API 合成的从头设计酶概念验证**。

## 🧬 在抗体药物领域能做什么

> 说明：RFD4-Proteina 是「生成 + 打分 + 搜索」的早期发现工具，输出为计算设计，**必须经湿实验验证**；它不是完整的抗体药物开发流水线，也不直接产出可申报的 IgG 分子。

### 1. VHH / 纳米抗体 CDR 设计（最直接）

- 固定 VHH 框架与靶点，模型**逐候选重采样 CDR1/CDR2/CDR3 长度与序列**，框架在空间上可相对靶点移动。
- 每个候选自动采样靶点附近的 hotspot patch，适合快速生成大量纳米抗体候选。
- 适合做纳米抗体发现、表位偏好筛选、多靶点并行探索。

### 2. Fab / IgG 的 CDR 重设计

- 在**固定 Ig 框架 + 固定抗原**条件下重设计 CDR（如 HER2 CDR H3 示例），用于：
  - 亲和力成熟 / 表位改换（epitope resurfacing）；
  - 现有抗体骨架上的 CDR 多样化，降低对天然免疫库的依赖；
  - 抗原特异性优化（可结合 `partial_diffusion` 对已有设计做局部精修）。

### 3. 从头蛋白 binder 作为抗体替代/补充

- 针对 PD-L1、TrkA 等治疗相关靶点设计 50–120 残基的 mini-binder：
  - 体积小、可穿组织、易表达、成本低，可作纳米抗体/抗体片段替代；
  - 支持**表位特异**设计（hotspot 指定），可绕开保守表位、追求功能阻断或别构。
- 论文中的 pMHC binder、磷酸化特异 binder、IDR binder 展示了「抗体难做的靶点」方向：
  - **肽–MHC**（肿瘤/感染新抗原）——TCR 样结合；
  - **磷酸化肽**（如 H3 pS10/T11）——翻译后修饰特异；
  - **内在无序蛋白**（YAP1、Tau 等）——传统抗体难以稳定识别。

### 4. 糖基化与可开发性相关设计

- 糖基屏蔽（glycan shielding）设计可调控 N-糖基化加工（保留 Man9/Glc1Man9），与**疫苗佐剂、抗体 Fc 糖型、免疫原性调控**相关。
- `inverse_folding` / 固定骨架序列设计可用于**人源化/序列优化/稳定性改造**的候选重采样（非专门人源化模块，需配合其他工具评估）。

### 5. 亲和力成熟与设计精修

- `partial_diffusion`：从已有抗体/结合体结构出发局部加噪再生成，可在保留关键接触的同时优化界面、提升亲和力或增加多样性。
- 结合 reward-guided search，用 AF3/OpenFold3 打分或自定义能量指标在生成轨迹上做 beam search，直接朝「像天然复合物、界面互补好」的方向搜索。

### 6. 酶/ADC 相关（抗体偶联与工艺）

- 论文的酶设计（磷酸三酯酶、连接酶）可用于**抗体偶联化学、linker 合成、API 工艺**等场景的从头酶开发——例如 lenacapavir 案例说明可做药物合成步骤的酶催化替代。

### 7. 局限与注意

- `antibody_design` 是**在已有框架上重设计 CDR**，不是从零生成完整 IgG；需要提供框架与抗原结构。
- 输出通过率是「计算标准」，**结合亲和力、免疫原性、可开发性、稳定性仍需实验与专门工具评估**。
- 在线 Studio 与生成流程会经 **SecureDNA** 生物安全筛查；涉及受控序列需遵守相关法规。
- 模型仅支持 Linux + NVIDIA GPU，无 Windows/macOS 支持。

### 抗体领域建议工作流

```
靶点结构/表位选择 (CIF, hotspot)
   → rfd4-generate 生成候选 (VHH/Fab/mini-binder, best-of-N)
   → AF3 / OpenFold3 refold + 抗体专用指标 (rfd4-evaluate)
   → 聚类、选优、partial_diffusion 精修
   → 湿实验 (BLI/SPR/酵母展示) 验证
```

## 🖥️ 部署硬件要求

### 最低软件/系统

| 项目 | 要求 |
|:---|:---|
| 操作系统 | **Linux x86_64**（不支持 Windows / macOS） |
| Python | **3.12** |
| GPU | **NVIDIA GPU**，驱动兼容 **CUDA 12 或 CUDA 13** |
| CPU/内存（NIM 生产容器） | ≥ 12 核 CPU，≥ 16 GB RAM |
| 存储 | ≥ 60 GB NVMe（模型缓存 + 环境） |
| GPU 数量 | **每个设计请求用单卡**；多卡可并行服务多个请求 |

### 官方测试/验证过的 GPU

| 架构 | 数据中心卡 | 工作站/桌面卡 |
|:---|:---|:---|
| Blackwell | B200、B300、GB200、GB300、RTX PRO 6000 Blackwell Server Edition | RTX PRO 6000 Blackwell Workstation Edition、GB10（DGX Spark） |
| Hopper | H100 SXM / PCIe / NVL、H200 / NVL、GH200 144GB | — |
| Ada Lovelace | L40S | RTX 6000 Ada Generation |
| Ampere | A100 SXM4 40GB/80GB、A100 PCIe 40GB/80GB | — |

### 显存（VRAM）要点

- 权重仅 **1.74 GB（FP32）**，但**峰值显存随 token 数近似平方增长、并随 batch size 线性增长**，与权重关系不大。
- `sampling.batch_size` 默认等于 `sampling.num_replicates`：多 replicate 请求会变成一次大前向。
- 官方测试卡最低为 **A100 40GB**；生产建议 **48–80 GB**（L40S / RTX 6000 Ada 48GB、A100 80GB、H100/H200）。
- 消费级 24GB 卡未列入官方支持清单，短序列/小 batch 或可运行，但不受支持、长设计易 OOM。
- 若启用 refolding oracle（AF2/AF3/OpenFold3），需额外独立环境/权重与显存。

### 存储与下载体积

| 内容 | 体积 |
|:---|:---|
| Python 环境（PyTorch + CUDA 库） | 约 7 GB |
| 模型权重 bundle（含 autoencoder） | 约 1.7 GB |
| AlphaFold2 参数（可选，refold/search） | 约 5.3 GB |
| OpenFold3 运行时镜像（可选） | 约 4 GB |
| OpenFold3 权重（可选） | 约 1.4 GB |
| 验证基准集（可选） | 约 63 MB |

### 部署方式

| 方式 | 命令/说明 |
|:---|:---|
| 源码（推荐） | `git clone --depth 1 … && cd RFD4-Proteina && make setup` |
| PyPI | `uv pip install --torch-backend cu132 "rfd4-proteina[all-cuda13]"`（CUDA 12 用 `all-cuda12` / `cu126`） |
| Docker | `docker build -t rfd4-proteina:cu13 .`；`--build-arg RFD4_ACCELERATOR=cu12|cpu` |
| Apptainer/Singularity | `apptainer build rfd4.sif docker-daemon://rfd4-proteina:cu13`，运行时加 `--nv` |
| SLURM 集群 | `./scripts/launch_generate.sh --njobs 4 --gpus-per-job 1 …` |
| NVIDIA NIM（生产 API） | NGC 拉取 `nvcr.io/nim/ipd/rfd4-proteina:latest`，单卡服务多个请求 |

### 快速上手（源码）

```bash
git clone --depth 1 https://github.com/RosettaCommons/RFD4-Proteina.git
cd RFD4-Proteina
make setup                 # 创建 .venv、缓存权重、配置 DATA_PATH
source .venv/bin/activate

# 一个 PD-L1 binder 示例（50-120 残基）
cat > binder.yaml <<'EOF'
binder:
  input: docs/examples/inputs/5o45_cropped.cif.zst
  contig: "50-120, /0, A17-131"
  conditions: {C_CRDGRP: {1: [{select: "sel('A/*/*/[N,CA,C,O]')"}]},
               C_SEQ: {true: [{select: "sel('A')"}]}}
EOF
rfd4-generate inputs=binder.yaml num_replicates=2 out_dir=./runs/first '~metrics'
```

> 首次冒烟测试可加 `inference_experiment=quickstart`（关指标、跳过 `torch.compile` 预热、只生成 1 个）。
> 设计输出在 `runs/first/inference_samples/`。

## 📌 结论（针对抗体药物）

- **值得跟进**：RFD4-Proteina 在 VHH/Fab CDR 设计、表位特异 mini-binder、难靶点（pMHC/磷酸化/无序蛋白）上能力突出，且通过率与速度同时领先，适合作为抗体/纳米抗体早期发现的**候选生成器**。
- **落地门槛**：必须是 Linux + NVIDIA 环境；建议 48GB 及以上单卡（A100 80GB / L40S / RTX 6000 Ada / H100 起步），配 ≥16GB 内存、≥60GB 存储；生产可用 NIM 容器，单卡服务多请求。
- **定位**：生成与计算筛选工具，**不能替代湿实验**；完整抗体开发仍需亲和力、特异性、稳定性、免疫原性、可开发性与工艺评估。

## 🛠️ 抗体亲和力成熟 Pipeline

`pipeline/` 把上面的模型能力落到可直接运行的抗体工程流程：

- 输入 VHH/Fab 序列或结构，自动用 ANARCI 识别 CDR；
- 按「有复合物 / 有靶点无复合物 / 无靶点」三种情形构造 annotated mmCIF；
- 调用 NVIDIA 托管 NIM（也可换成自托管端点）做 `antibody_design`；
- 输出中文 Markdown/HTML 报告、`designs.csv`、结构文件和湿实验建议。

```bash
cd tools/RFD4-Proteina/pipeline
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
export NVIDIA_API_KEY="nvapi-..."

# 只构建请求，不调用 API
python run.py --complex-structure complex.pdb --antibody-chains A \
  --target-chains T --epitope T:557-573 --dry-run

# 真正提交（VHH + 靶点结构）
python run.py --vhh-seq examples/vhh_sequence.fasta \
  --target-structure 8UCD --target-chains A --cdrs all --num-designs 8
```

详见 [`pipeline/README.md`](pipeline/README.md)。

## 🔗 参考链接

- 项目主页：https://rfd4-proteina.org/
- 文档：https://rfd4-proteina.org/latest/
- GitHub：https://github.com/RosettaCommons/RFD4-Proteina
- Hugging Face 权重：https://huggingface.co/nvidia/RFD4-Proteina
- NVIDIA NIM 支持矩阵：https://docs.nvidia.com/nim/bionemo/rfd4-proteina/latest/support-matrix.html
- NVIDIA NIM 性能：https://docs.nvidia.com/nim/bionemo/rfd4-proteina/latest/performance.html
- VHH 设计教程：https://rfd4-proteina.org/latest/examples/antibody_design_tutorial.html
- UW Medicine 新闻：https://newsroom.uw.edu/news-releases/researchers-achieve-major-advance-in-protein-design
