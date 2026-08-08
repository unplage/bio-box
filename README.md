# 🧬 Bio-Box

生物医学工具合集 | Biomedical Tools Collection

<img width="1917" height="1039" alt="image" src="https://github.com/user-attachments/assets/21d8e786-63b5-4ead-b87f-ee59498457c9" />
<img width="1102" height="744" alt="image" src="https://github.com/user-attachments/assets/d7dfaf35-be91-4053-94a3-4d96751945c1" />

## 📦 工具列表

| 工具 | 文件 | 描述 |
|:---|:---|:---|
| 🎯 靶点调研报告生成器 | `target.html` | PWA Web 应用，一键生成靶点调研报告 |
| 🧹 精准数据清理工具 | `clear.html` | Web 应用，用于清理和标准化数据 |
| 🧬 BsaI 酶切位点重复分析 | `bsaitest.html` | 分析4-mer在 BsaI 位点数据库中的出现频率 |
| 🧪 NNK 引物生成器 | `NNKprimer.html` | 输入编码区序列，生成 NNK 突变引物组合 |
| 🔬 2ED Maturation 反向引物 | `2edmaturation.html` | 截取 NNK 前24bp 计算 Tm，生成反向互补引物 |
| 📊 CDR 序列提取工具 | `CDR-extract-V1.2.1-260618.py` | 从 Excel/CSV 提取 CDR 序列 |
| 🔍 CDR3/CDR2/CDR1 搜索工具 | `igblast-solo-SQLite260618-V1.3.2.py` | SQLite 数据库搜索 CDR 区域 |
| 🧮 引物 Tm 值计算器 | `primer_Tm-calculator.py` | 批量计算引物 Tm 值 |
| 🎯 靶点调研（Python 版） | `target_info.py` | PyQt5 桌面版靶点调研工具 |
| 🦠 菌落计数 | [cfu_count](https://github.com/unplage/cfu_count) | 菌落计数工具（独立仓库） |

---

## 🎯 靶点调研报告生成器 (PWA)

**target.html** 是一个纯前端 PWA Web 应用，可直接在浏览器中运行，无需安装任何软件。

### ✨ 功能特性

- **多数据源整合**：PubMed、ClinicalTrials.gov、Open Targets、UniProt、ChEMBL、Human Protein Atlas、USPTO/Lens.org 专利
- **AI 智能分析**：支持 DeepSeek、小米 MiMo、智谱 GLM、自定义 OpenAI 兼容接口
- **联网情报分析**：智谱 GLM 支持实时联网搜索获取最新行业情报
- **一键生成报告**：Markdown、JSON、PPT 多格式导出
- **PWA 支持**：可安装到手机/桌面，支持离线使用

### 🚀 快速开始

#### 方式一：直接访问（推荐）

如果仓库已部署到 GitHub Pages，直接访问：
```
https://unplage.github.io/bio-box/target.html
```

#### 方式二：本地使用

1. 克隆仓库
```bash
git clone https://github.com/unplage/bio-box.git
cd bio-box
```

2. 用浏览器打开 `target.html`
```bash
# macOS
open target.html

# Linux
xdg-open target.html

# Windows
start target.html
```

#### 方式三：部署到 GitHub Pages

1. Fork 本仓库
2. 进入仓库 Settings → Pages
3. Source 选择 `main` 分支，目录选择 `/ (root)`
4. 点击 Save，等待部署完成

### 📱 安装为 PWA

#### 手机端（iOS/Android）

1. 用 Chrome/Safari 打开 `target.html`
2. 点击浏览器菜单中的「添加到主屏幕」或「安装应用」
3. 确认安装，即可像 App 一样使用

#### 桌面端（Windows/Mac/Linux）

1. 用 Chrome/Edge 打开 `target.html`
2. 点击地址栏右侧的安装图标 📥
3. 或进入菜单 → 更多工具 → 创建快捷方式
4. 即可在桌面创建独立应用

### ⚙️ 配置说明

#### 1. AI 供应商设置

| 供应商 | API Key 获取地址 | 联网搜索 |
|:---|:---|:---:|
| DeepSeek | [api.deepseek.com](https://api.deepseek.com) | ❌ |
| 小米 MiMo | [api.xiaomimimo.com](https://api.xiaomimimo.com) | ❌ |
| 智谱 GLM | [open.bigmodel.cn](https://open.bigmodel.cn) | ✅ |
| 自定义 | OpenAI 兼容接口 | ❌ |

> 💡 **提示**：选择智谱 GLM 时，AI 综合情报研判会自动联网搜索最新行业信息

#### 2. CORS 代理（可选）

如果遇到跨域请求失败，可配置 CORS 代理：
- 推荐代理：`https://corsproxy.io/?`
- 在设置中填入代理前缀即可

#### 3. 专利数据源

| 数据源 | 说明 | 需要 API Key |
|:---|:---|:---:|
| USPTO | 美国专利商标局 | ✅ |
| Lens.org | 全球专利检索 | ✅ |

### 📋 使用流程

```
1. 输入靶点名称（如：PD-1、EGFR、HER2）
   ↓
2. 选择数据源（文献、临床、靶点、药物、专利）
   ↓
3. 配置 AI 供应商和 API Key
   ↓
4. 点击「🚀 开始生成」
   ↓
5. 等待数据检索和 AI 分析（约 1-2 分钟）
   ↓
6. 查看报告并导出（Markdown / JSON / PPT）
```

### 📊 报告内容

生成的报告包含以下章节：

1. **🧬 靶点概述与背景** - 基因符号、蛋白类别、功能描述、相关疾病
2. **🧬 表达谱与抗体可及性** - HPA 数据、蛋白定位、抗体信息
3. **📈 研究进展** - 文献数量、年份趋势、核心文献
4. **📚 文献列表** - PubMed 检索结果详情
5. **🧪 临床试验** - ClinicalTrials.gov 试验数据
6. **💊 药物研发管线** - Open Targets 药物信息
7. **📜 专利调研** - USPTO/Lens.org 专利数据
8. **🧠 AI 综合情报研判** - AI 分析与建议
9. **🔬 靶向活性分子** - ChEMBL 化合物线索

### 🎨 PPT 导出功能

支持将报告导出为 PowerPoint 格式，包含：

- 封面页（靶点名称、生成日期）
- 靶点概述页
- 表达谱数据表格
- 研究进展与文献趋势
- 临床试验概况
- 药物管线列表
- 专利调研与 AI 解读
- 数据来源致谢

> 💡 PPT 使用 [PptxGenJS](https://github.com/gitbrent/PptxGenJS) 库生成，兼容 PowerPoint、Keynote、LibreOffice、Google Slides

---

## 🧬 BsaI 酶切位点重复分析

`bsaitest.html` 分析 DNA 序列中每个4-mer在 BsaI 位点数据库中的出现频率。

### 功能

- 输入 DNA 序列，滑动4字符窗口提取所有4-mer
- 统计每个4-mer及其反向互补在数据库中的出现次数
- 显示所有位点（含0次），支持下载结果

### 使用

直接用浏览器打开 `bsaitest.html`，输入序列点击「分析」即可。

---

## 🧪 NNK 引物生成器

`NNKprimer.html` 输入编码区序列，自动生成 NNK 突变引物组合。

### 功能

- 输入 Upper Prefix、编码区序列、突变位点数、Later Suffix
- 使用组合数学生成所有可能的 NNK 替换位置
- 支持下载 `primer_list.txt`

### 使用

直接用浏览器打开 `NNKprimer.html`，填写参数点击「Generate」即可。

---

## 🔬 2ED Maturation 反向引物生成器

`2edmaturation.html` 处理含 NNK 的序列，自动生成反向引物并计算 Tm 值。

### 功能

- 输入含 NNK 的 DNA 序列（每行一条）
- 截取 NNK 前24个碱基，使用 SantaLucia 2004 热力学参数计算 Tm
- 生成从截取起始位到序列末尾的反向互补序列
- 支持下载 `primer_list.txt` 和 `tm_values.txt`

### Tm 计算参数

- 校准基于56条序列（RMSE: 1.16°C）
- 默认条件：Na⁺ 50mM，引物浓度 0.25µM

---

## 🧹 精准数据清理工具

`clear.html` 是一个 Web 数据清理工具，用于：

- 数据去重
- 格式标准化
- 缺失值处理
- 异常值检测

---

## 🐍 Python 工具

### CDR 序列提取工具

```bash
python CDR-extract-V1.2.1-260618.py
```

功能：
- 从 Excel/CSV 文件提取 CDR1/CDR2/CDR3 序列
- 支持 .xlsx 和 .csv 格式
- 自动识别编码和分隔符

### CDR3/CDR2/CDR1 搜索工具

```bash
python igblast-solo-SQLite260618-V1.3.2.py
```

功能：
- SQLite 数据库存储和搜索
- 支持多区域组合搜索
- FASTA 批量分析

### 引物 Tm 值计算器

```bash
python primer_Tm-calculator.py
```

功能：
- 基于校准参数的 Tm 值计算
- 单序列和批量处理
- GC 含量和退火温度建议

---

## 🛠️ 技术栈

### target.html (PWA)

- **前端**：纯 HTML/CSS/JavaScript，无框架依赖
- **PWA**：支持安装到设备，可离线使用
- **API**：PubMed、ClinicalTrials.gov、Open Targets、UniProt、ChEMBL、HPA
- **AI**：OpenAI 兼容接口（DeepSeek/MiMo/智谱）
- **导出**：Markdown、JSON、PPT (PptxGenJS)

### Python 工具

- **GUI**：tkinter / PyQt5
- **数据处理**：pandas
- **异步**：asyncio + httpx

---

## 📝 更新日志

### v2.1 (2026-08-08)

- ✨ 新增 `bsaitest.html` - BsaI 酶切位点重复分析
- ✨ 新增 `NNKprimer.html` - NNK 引物生成器
- ✨ 新增 `2edmaturation.html` - 2ED Maturation 反向引物生成器
- 🐛 修复 `target.html` PDB 搜索查询逻辑

### v2.0 (2026-08-07)

- ✨ 新增 PPT 导出功能
- ✨ 新增 AI 联网情报分析（智谱 GLM）
- 🐛 修复专利搜索 403 错误
- 🐛 修复 MiMo 供应商联网搜索问题

### v1.0

- 🎉 初始版本发布
- 📊 整合多数据源靶点调研
- 🤖 支持多 AI 供应商

---

## 📄 许可证

MIT License

---

## 🤝 贡献

欢迎提交 Issue 和 Pull Request！

## 📧 联系方式

- GitHub: [unplage](https://github.com/unplage)
