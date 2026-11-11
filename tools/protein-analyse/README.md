# 🧬 Protein-Analyse 蛋白质分析工具集

单文件纯前端 PWA，浏览器直接打开即可离线使用，无需构建、无外部依赖。

| 文件 | 说明 |
|:---|:---|
| `index.html` | 主页面（PWA，全部算法内联） |
| `manifest.json` | PWA manifest，可安装到桌面/手机 |
| `sw.js` | Service Worker，离线缓存 |
| `icon.svg` | 应用图标 |

## 功能模块

| 模块 | 功能 | 算法来源 |
|:---|:---|:---|
| ⚖️ 分子量 / pI | 分子量、等电点、氨基酸组成表、**多 pKa 参数集对比**、Cys 二硫键处理、PTM 提示 | ExPASy 残基平均质量；Henderson–Hasselbalch + 二分法求净电荷=0；可选 Bjellqvist/IPC_protein/EMBOSS/Sillero/Nozaki–Tanford |
| 🌡️ KD 疏水性 | 平均 KD、滑动窗口剖面图、残基热图 | Kyte & Doolittle (1982) J Mol Biol |
| 💧 抗体粘度 | 输入 VH/VL，预测 180 mg/mL 高浓度粘度、HI、Fv 净电荷与电荷不对称性，并给出风险分层 | Sharma et al. (2014) PNAS 序列模型；Eisenberg (1984) 疏水标度 |
| 🔬 A280 浓度 | ε₂₈₀、mg/mL、mM/μM、1 mg/mL 吸光度 | Pace et al. (1995)：`ε = 5500·W + 1490·Y + 125·Cystine`；Beer-Lambert |
| 🧪 摩尔换算 | 分子的量 (nmol/pmol)、摩尔浓度、质量/浓度互转 | `n = m/MW`，`c = n/V`，1 g/L ≡ 1 mg/mL |
| 🧮 Tm 计算 | 引物 Tm、GC%、退火温度建议；支持单序列/批量 FASTA/每行一条序列与 CSV 导出 | 完全复刻 `tools/Tm-calculate/primer_Tm-calculator.py`：SantaLucia 最近邻参数 + 校准补偿 + Na⁺ 盐校正 |
| 📚 参考手册 | 遗传密码子表、20 种氨基酸性质表（含分子式/结构式/残基质量）、IUPAC 简并碱基表、核酸定量换算、缓冲液/培养基、抗生素/母液、蛋白标签/蛋白酶、限制性内切酶、单位换算、荧光蛋白/报告基因；支持搜索与分类筛选 | NCBI 标准遗传密码（transl_table=1）；IUPAC 碱基歧义代码；ExPASy 残基质量；常用实验配方与试剂参数 |
| 🧬 密码子优化 | 按 7 个物种的密码子偏好生成编码序列，排除酶切位点、控制 GC 窗口/聚聚体/CpG | Sharp & Li CAI + beam search（K=60）+ 局部同义修复；Kazusa CUD 频率表 |

## 等电点（pI）算法与可靠性

**计算框架**：Henderson–Hasselbalch 方程 + 二分法。对序列中可离子化基团逐类计算电荷：

```
Q(pH) = Σ_正 1/(1+10^(pH−pKa)) − Σ_负 1/(1+10^(pKa−pH))
pI    = Q(pH) = 0 的根（Q 随 pH 单调递减，二分区间 0–14）
```

正电荷基团：N 端、His、Lys、Arg；负电荷基团：C 端、Asp、Glu、Cys、Tyr。所有 pKa 参数集共用这一框架，差异只来自 pKa 取值。

**可选 pKa 参数集**（页面下拉框切换，结果对比表同时列出全部）：

| 参数集 | 来源 | 说明 |
|:---|:---|:---|
| Bjellqvist / ExPASy（默认） | Bjellqvist et al. 1993/1994 | 与 ExPASy ProtParam、Biopython `IsoelectricPoint` 逐位一致；含首/末残基末端修正（如 N 端 M→7.0、C 端 D→4.55） |
| IPC_protein | Kozlowski 2016, Biology Direct | 基于 2324 个实验 pI 优化的九参数集；经典 pKa 集中蛋白基准最优（RMSE ≈ 0.87） |
| EMBOSS | EMBOSS `iep` | 经典九参数集 |
| Sillero | Sillero & Ribeiro 1989 | 经典九参数集 |
| Nozaki–Tanford | Nozaki & Tanford 1967 | 经典九参数集 |

**Cys 二硫键处理**：可选①全部按游离巯基计入；②每 2 个 Cys 形成 1 对二硫键，仅未配对的奇余 Cys 带电；③全部视为氧化，忽略 Cys 电荷。成键 Cys 的巯基不再电离。

**PTM 提示**：结果面板会提示潜在 N-糖基化序列（N-X-S/T，X≠P）与可磷酸化的 S/T/Y 数量，但**不参与计算**。

**可靠性**：本工具给出的是序列层面的理论估算，未考虑三维结构微环境（埋藏、盐桥、氢键）、金属离子、翻译后修饰及非标准残基。与实验 pI 的典型偏差约 **±0.5–1.5 pH**，不能替代等电聚焦等实验测定。蛋白基准 RMSE 参考（Kozlowski 2016/2021）：IPC_protein ≈ 0.87、Bjellqvist ≈ 0.94、Nozaki ≈ 1.02、Sillero ≈ 1.05；机器学习方法 IPC 2.0 可降至 ≈ 0.85，但需联网/外部服务，本离线工具未采用。

## 抗体粘度（Sharma 模型）

**适用对象**：单克隆抗体 Fv/Fab（分别输入 VH 与 VL 可变区）。预测条件为 **180 mg/mL、25°C、pH 5.5、200 mM 精氨酸盐酸盐**。

**模型**（Sharma et al. 2014, *PNAS* 111:18601）：

```
η (cP) = 10^[0.15 + 1.26·HI − 0.043·q − 0.020·q_sym]
```

| 特征 | 定义 |
|:---|:---|
| `q` | Fv 在 pH 5.5 的净电荷 = q_VH + q_VL；仅计 Asp/Glu/His/Lys/Arg，pKa 取 3.9/4.3/6.0/10.8/12.5（Berg 等），按 Henderson–Hasselbalch 计算 |
| `q_sym` | 电荷不对称性（FvCSP）= q_VH × q_VL |
| `HI` | 疏水指数 = −Σ_疏水 E_i / Σ_亲水 E_j；疏水残基 A/C/F/I/L/P/V/W/Y，亲水残基 D/E/G/H/K/M/N/Q/R/S/T |
| `E` | Eisenberg 归一化共识疏水标度（Eisenberg et al. 1984）：A 0.62、C 0.29、F 1.19、G 0.48、I 1.38、L 1.06、M 0.64、P 0.12、V 1.08、W 0.81、Y 0.26、D −0.90、E −0.74、H −0.40、K −1.50、N −0.78、Q −0.85、R −2.53、S −0.18、T −0.05 |

**风险分层**（仅工程参考）：`<10 cP` 低、`10–20 cP` 中、`>20 cP` 偏高。20 cP 阈值借自 DeepViscosity 的 150 mg/mL 判据，本模型为 180 mg/mL，二者不可直接等同。

**示例**：曲妥珠单抗 Fv（VH 120 aa / VL 107 aa）→ HI ≈ 1.127、q ≈ 5.15、q_sym ≈ 6.43、η ≈ 16.6 cP；文献实测高浓度约 8.6 cP，说明本模型存在低估/偏差（Sharma 原文 MAE 约 7 ± 9 cP）。

**与 DeepViscosity 的关系**：DeepViscosity（Kalejaye/Lai 2025）是 102 个 ANN 集成 + 30 个 DeepSP CNN 描述符的数据驱动分类器，需 ANARCI/IMGT 比对和外部模型权重，输出 150 mg/mL 下低（≤20 cP）/高（>20 cP）分类，准确率约 87.5–89.5%。它**无法在离线单文件网页中复现**，本工具明确不实现、也不冒充其结果；这里采用可解释、可离线复算的 Sharma 序列模型。

**局限**：未考虑三维结构、盐/辅料、聚集、翻译后修饰；训练集有限，对部分抗体偏差较大（曲妥珠单抗实测约 8.6 cP，本模型约 16.6 cP）。高浓度制剂决策仍需实验流变学验证。

## 密码子优化

**数据**：Kazusa Codon Usage Database（Nakamura et al. *NAR* 2003，GenBank R160）——
大肠杆菌(8087 CDS)、人(93487)、小鼠(53036)、酿酒酵母(14411)、枯草芽孢杆菌(2529)、CHO(331)、毕赤酵母(137)；
每张表的 64 密码子频率、编码区 GC%、CDS 数均已与线上数据逐位核对。

**算法**（详见页面「算法说明」卡片）：

1. 目标函数 `score = Σ ln w(cᵢ) − P_位点 − P_聚聚体 − P_GC窗口 − P_CpG`，其中 `w(c)=f(c)/max f(同义)`（Sharp & Li 1987 CAI 相对适应度）；
2. beam search（宽度 60）逐密码子扩展同义候选并按累计得分剪枝——无约束时可证明达到贪心最优（测试验证与逐位最高频序列 CAI 完全一致）；
3. 限制性内切酶位点为**硬约束**（当步禁止生成，自动含反向互补链，如 BsaI 同时排除 GGTCTC/GAGACC）；
4. 结束后对残余违规做局部同义替换修复（位点 → 聚聚体 → GC 窗口，迭代多轮）；
5. 输出前强制**翻译回原文校验**，结果面板逐项展示校验状态。

**选项**：42 种常用酶（含 Golden Gate 的 BsaI/BbsI/BsmBI，以及 **SfiI/XmnI/BglI/BstEII/HinfI/DdeI 等含 N 简并碱基的位点**）+ 自定义位点（支持 IUPAC 简并码 RYSWKMBDHVN）；
GC 窗口上下限/窗长；避免 ≥5 聚聚体、避免 CpG；终止密码子自动（取该物种最常用）或手动 TAA/TAG/TGA/不追加。

**位点排除机制**：
1. 位点编译为 IUPAC 正则（`N→[ACGT]`、`R→[AG]` 等），如 SfiI `GGCCNNNNNGGCC` 一次匹配全部 4⁵ 种中间变体；同时自动追加反向互补链；
2. beam search 硬约束：每放置一个密码子，检测序列末段是否**新完成**任一位点（含简并变体），命中即拒绝该同义候选；
3. 所有同义候选均冲突时降级计罚（-300）保证可继续；
4. 生成后局部同义替换修复残余位点 → 输出前全串复扫校验并在结果面板展示。

**评价指标**：CAI（不含终止子）、整体 GC%、窗口 GC 范围、最长聚聚体、位点检出列表。

## 参考手册（查询与学习）

- **标准遗传密码子表**：16 行 × 4 列矩阵——行 = 第 1/2 位碱基，列 = 第 3 位碱基；每格同时显示密码子、三字母简写与单字母，浅色底纹区分 20 种氨基酸；3 个终止密码子红色高亮、起始密码子 AUG★ 橙色描边；
- **20 种氨基酸性质表**：1/3 字母简写、中英文名、**分子式（数字转下标，如 C₃H₇NO₂）**、**结构式（等宽字体绘制主链与侧链键连）**、**残基质量 (Da)**、**pI (等电点)**、极性/酸碱分类（非极性、极性非带电、酸性 D/E、碱性 H/K/R）与性质说明，按疏水→亲水→酸→碱排序；
- **简并碱基符号表**：IUPAC 16 个符号（A/C/G/T/U + R/Y/S/W/K/M/B/D/H/V/N）及所含碱基、中文名、引物设计示例；
- **核酸与定量换算**：碱基配对/氢键数、A260→浓度（dsDNA 50 / ssDNA 33 / RNA 40 µg/mL）、纯度比值、质量↔pmol↔拷贝数、DNA/RNA 分子量、oligo Tm；
- **缓冲液与培养基**：50× TAE、10× TBE/PBS/TBS、TE、SDS-PAGE/转膜/上样缓冲液、EDTA/Tris/NaCl/NaOAc/SDS 母液；LB、**2×YT**、TB、SOC、M9 培养基；
- **抗生素与母液**：Amp/Kan/Cm/Tet/Spec/Hyg/Puro/G418/Zeocin/Blasticidin 的母液、工作浓度、储存；IPTG/X-gal/DTT/β-ME/PMSF/EDTA、蛋白酶抑制剂、RNase/DNase/溶菌酶/蛋白酶 K；
- **蛋白标签与蛋白酶**：His6/GST/MBP/SUMO/FLAG/HA/Myc/V5/Strep/CBP/GFP 的纯化与洗脱；TEV/PreScission/Thrombin/Factor Xa/肠激酶/SUMO 酶切位点与缓冲液；
- **限制性内切酶**：33 种常用酶（含 BsaI/BbsI/BsmBI Golden Gate、SfiI/XmnI/BglI/BstEII/HinfI/DdeI 简并位点）的识别序列、切点、末端与同尾酶备注；
- **单位换算与实验计算**：SI 词头、质量/体积/浓度/稀释、% w/v、ppm、rpm↔rcf、A260、温度/压力/酶活、DNA→蛋白估算、母液稀释表；
- **荧光蛋白与报告基因**：常用荧光蛋白 Ex/Em、荧光素酶/β-gal/SEAP 报告基因、常用抗性基因；
- **🧮 引物 Tm 计算（第 7 面板）**：完全复刻 `tools/Tm-calculate/primer_Tm-calculator.py`——最近邻热力学模型 + 校准参数（56 条序列，RMSE 1.16°C）、Na⁺ 盐校正、末端 AT 修正、退火模式对称因子；输出 Tm、GC%、标准酶/高保真酶/梯度 PCR 退火建议；支持 FASTA 批量、每行一条纯序列与 CSV 导出；
- 顶部搜索框实时过滤当前子表；分类 chip 仅在氨基酸表显示；密码子表搜索时非匹配项淡显。
- 酶切位点/缓冲液/荧光蛋白等为常用速查数据，厂商、批次和实验条件可能存在差异，关键实验请核对试剂说明书、FPbase 或最新文献。

## 使用方法

```bash
# 本地直接打开
open tools/protein-analyse/index.html   # macOS
# 或部署到 GitHub Pages 后作为 PWA 安装
```

1. 粘贴 FASTA 或纯氨基酸序列（自动过滤非标准字符）；
2. 切换标签页，点「⚡ 计算」得到结果；
3. 「⚖️ 分子量 / pI」可选择 pKa 参数集与 Cys 二硫键处理方式；结果下方同时给出全部参数集的 pI 对比、Cys/PTM 提示与可靠性说明；
4. A280 模块可配置光程与二硫键假设（氧化 Cystine ×125 / 还原不计）；
5. 摩尔换算优先取「质量浓度」，否则用「质量 + 体积」；
6. 「📚 参考手册」无需输入序列，可直接查询密码子/氨基酸/简并碱基，以及核酸定量、缓冲液/培养基、抗生素/母液、蛋白标签/蛋白酶、限制性内切酶、单位换算、荧光蛋白/报告基因；
7. 「🧬 密码子优化」读取同一序列输入：选物种 → 点选要排除的酶切位点 → ⚡ 开始优化 → 下载 FASTA。
8. 「💧 抗体粘度」在独立输入框分别粘贴 VH/VL 可变区（FASTA 或纯序列）→ ⚡ 计算粘度；示例按钮可载入曲妥珠单抗 Fv 对照。

## 导出说明

- 分子量/pI 结果支持一键复制；
- 全部计算在浏览器本地完成，序列不上传、不联网（Service Worker 离线可用）。
- Service Worker 使用独立缓存（`protein-analyse-*`），更新时只清理本工具缓存，不影响同域下其他 PWA 的离线缓存。

## 核心公式

```
MW      = Σ(残基平均质量) + 18.0153
pI      : 二分法求 netCharge(pH) = 0（Henderson-Hasselbalch）
Q(pH)   = Σ_正 1/(1+10^(pH−pKa)) − Σ_负 1/(1+10^(pKa−pH))
          正: N端/His/Lys/Arg；负: C端/Asp/Glu/Cys/Tyr
ε280    = 5500·n(Trp) + 1490·n(Tyr) + 125·n(Cystine)   [M⁻¹cm⁻¹]
c(M)    = A280 / (ε · l)
c(mg/mL)= A280 · MW / (ε · l)      # 因 1 g/L ≡ 1 mg/mL
n(mol)  = m(g) / MW
c(M)    = n / V(L)
η(cP)   = 10^[0.15 + 1.26·HI − 0.043·q − 0.020·q_sym]
          q = q_VH + q_VL（pH 5.5）；q_sym = q_VH·q_VL
          HI = −Σ疏水 E / Σ亲水 E（Eisenberg 1984）
```
