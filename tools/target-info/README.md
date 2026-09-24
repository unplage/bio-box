# 靶点调研报告生成器 (TargetInfo)

一键检索多数据源并生成靶点调研报告，含三个版本。

## 文件说明

| 文件 | 类型 | 说明 |
|:---|:---|:---|
| `target.html` | 纯前端 PWA | 推荐入口，浏览器直接打开即可使用，可安装到桌面/手机 |
| `target.py` | PyQt5 桌面版 | 功能最全的桌面版，无 CORS 限制，集成更多专业数据库 |
| `target_info.py` | PyQt5 桌面版 | 早期桌面版 |

## 使用方法

### target.html（推荐）

```bash
# macOS
open tools/target-info/target.html
# Linux
xdg-open tools/target-info/target.html
# Windows
start tools\target-info\target.html
```

1. 输入靶点名称（如 `PD-1`、`EGFR`、`HER2`，也支持 `EGFR T790M` 突变特异检索）
2. 勾选数据源（文献 / 临床 / 靶点 / 药物 / 专利）
3. 配置 AI 供应商与 API Key（DeepSeek / 小米 MiMo / 智谱 GLM / 自定义 OpenAI 兼容接口，Key 仅存本地浏览器）
4. 点击「🚀 开始生成」，等待约 1–2 分钟
5. 导出报告：Markdown / JSON / PPT / Word / HTML / 打印 PDF

其他功能：

- **历史报告库**：基于 IndexedDB 自动保存，支持查看、下载、删除、多记录对比
- **追问问答**：基于当前报告的本地 BM25 检索 + AI 作答并标注来源
- **PWA**：可「添加到主屏幕」安装，支持离线使用
- **CORS 代理**：专利检索内置 r.jina.ai / allorigins 自动兜底，可填自定义代理

GitHub Pages 部署后可直接访问：
`https://unplage.github.io/bio-box/tools/target-info/target.html`

### target.py（桌面版）

```bash
pip install PyQt5 httpx python-pptx python-dotenv pydantic
python tools/target-info/target.py
```

数据源覆盖：PubMed / OpenAlex / Semantic Scholar、ClinicalTrials.gov / ISRCTN / ChiCTR、Open Targets / UniProt / ClinVar / KEGG、PDB / AlphaFold / STRING、HPA / GTEx、ChEMBL / PubChem、专利与 AI 分析（同网页版）。支持导出 Markdown / JSON / PPT 等。

### target_info.py（早期桌面版）

```bash
pip install PyQt5 httpx python-pptx python-dotenv pydantic matplotlib
python tools/target-info/target_info.py
```

依赖 `.env` 可配置 `PUBMED_EMAIL`（PubMed 检索要求）。

## 配置说明

| 配置项 | 说明 |
|:---|:---|
| AI API Key | 页面内填写，存浏览器 localStorage / 桌面版本地配置，不会上传 |
| NCBI API Key | 可选，提升 PubMed 请求配额 |
| CORS 代理 | 可选，专利源跨域时使用，内置免费代理链自动兜底 |
| 专利数据源 | Google Patents（免 Key）、USPTO、Lens.org、Espacenet、MCP |

详细功能与报告章节见仓库根目录 [README.md](../../README.md)。
