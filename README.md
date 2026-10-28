# 🧬 Bio-Box

生物医学工具合集 | Biomedical Tools Collection

## 📦 工具列表

| 工具 | 路径 | 描述 |
|:---|:---|:---|
| 🎯 靶点调研报告生成器 | `tools/target-info/target.html` | PWA Web 应用，多数据源整合，AI 智能分析，一键生成报告 |
| 🧹 精准数据清理工具 | `clear.html` | Web 应用，数据清理与标准化 |
| 🧬 BsaI 酶切位点重复分析 | `tools/maturation/bsaitest.html` | 4-mer 在 BsaI 位点数据库中的出现频率分析 |
| 🧪 NNK 引物生成器 | `tools/maturation/NNKprimer.html` | 输入编码区序列，生成 NNK 突变引物组合 |
| 🧬 His-mut 多点组氨酸突变引物 | `tools/his-mut/Hismut.html` | 从 FASTA 定位突变点，生成组氨酸突变引物 |
| 🔬 2ED Maturation 反向引物 | `tools/maturation/2edmaturation.html` | 截取 NNK 前24bp 计算 Tm，生成反向互补引物 |
| 📊 CDR 序列提取工具 | `tools/mab-db/CDR-extract-V1.2.1-260618.py` | 从 Excel/CSV 提取 CDR 序列 |
| 🔍 CDR 搜索工具 | `tools/mab-db/igblast-solo-SQLite260618-V1.3.2.py` | SQLite 数据库搜索 CDR 区域 |
| 🧮 引物 Tm 值计算器 | `tools/Tm-calculate/primer_Tm-calculator.py` | 批量计算引物 Tm 值 |
| 🎯 靶点调研（Python 版） | `tools/target-info/target_info.py` | PyQt5 桌面版靶点调研工具 |
| 🤖 **RFD4-Proteina 亲和力成熟 Pipeline** | `tools/RFD4-Proteina/` | CDR 识别、NIM 调用、报告生成、湿实验建议 |
| 🤖 **NVIDIA Bio Tools** | `tools/nvidia-bio/` | ProteinMPNN / Boltz-2 统一 CLI 与 Python API |

---

## 🎯 靶点调研报告生成器 (PWA)

`tools/target-info/target.html` - 纯前端 PWA，浏览器直接运行，支持离线。

**核心能力**：
- 多数据源：PubMed、ClinicalTrials、Open Targets、UniProt、ChEMBL、HPA、专利
- AI 分析：DeepSeek、MiMo、智谱 GLM、自定义 OpenAI 兼容接口
- 联网搜索：智谱 GLM 支持实时联网
- 多格式导出：Markdown、JSON、PPT、Word、HTML、PDF
- PWA 支持：可安装到手机/桌面，离线可用

**快速开始**：
```bash
# 本地使用
git clone https://github.com/unplage/bio-box.git
cd bio-box
open tools/target-info/target.html  # macOS
# 或直接部署到 GitHub Pages
```

---

## 🧬 RFD4-Proteina 亲和力成熟 Pipeline

`tools/RFD4-Proteina/` - CDR 识别、NIM 调用、结构生成、报告生成。

**核心功能**：
- 自动 CDR 识别 (ANARCI/IMGT/Kabat/Chothia/Martin)
- 三种模式：复合物 / 有靶点无复合物 / 仅抗体
- NVIDIA RFD4-Proteina NIM 调用 (托管/自托管)
- 中文 Markdown/HTML 报告 + CSV/JSON + 结构文件
- 湿实验建议 + 开发性风险评估

**关键修复**：
- ✅ CDR 编号位置→序列位置修正 (修复 CDR3 吃掉 FR4)
- ✅ Expandable segment 编码 (官方 antibody 示例风格，CDR 可变长)
- ✅ 默认 metric tags 仅 `design` + 条件类，避免 422
- ✅ CDR 长度范围控制 (`--cdr-length-range`)
- ✅ API key 自动脱敏，不落盘

---

## 🤖 NVIDIA Bio Tools

`tools/nvidia-bio/` - 统一 CLI 与 Python API。

```bash
# ProteinMPNN
python main.py proteinmpnn input.pdb --api-key $KEY -o out.fa
python main.py proteinmpnn input.pdb --self-hosted  # 自托管无需 key

# Boltz-2
python main.py boltz2 --sequence "MALWMRLL..." --api-key $KEY -o out/
python main.py boltz2 --chains chains.json --api-key $KEY -o out/
python main.py boltz2 --sequence "..." --self-hosted -o out/
```

⚠️ **托管 API 限制**：仅返回序列/结构，**无置信度/概率/熵**，生产建议自托管。

---

## 🔬 其他工具

- **BsaI 酶切分析** / **NNK 引物** / **His-mut** / **2ED Maturation** - 浏览器打开即用
- **CDR 提取/搜索** / **Tm 计算** - Python 脚本运行
- **精准数据清理** - `clear.html` 浏览器打开

---

## 🛠️ 技术栈

| 类别 | 工具 |
|------|------|
| Web/PWA | 纯 HTML/CSS/JS，PWA 支持，IndexedDB |
| Python GUI | tkinter / PyQt5 |
| 数据处理 | pandas, biotite |
| NIM 客户端 | requests, httpx, asyncio |
| 结构处理 | biotite, mmCIF/PDB 解析 |

---

## 📄 许可证

MIT License

---

## 🤝 贡献

欢迎提交 Issue 和 Pull Request！

---

## 📧 联系方式

- GitHub: [unplage](https://github.com/unplage)