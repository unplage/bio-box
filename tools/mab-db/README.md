# 单抗数据库工具 (mab-db)

抗体 CDR 序列提取与 SQLite 检索工具（tkinter 图形界面）。

## 文件说明

| 文件 | 说明 |
|:---|:---|
| `CDR-extract-V1.2.1-260618.py` | 从 Excel/CSV 批量提取 CDR1/CDR2/CDR3 氨基酸序列，导出 FASTA |
| `igblast-solo-SQLite260618-V1.3.2.py` | 将 CDR 数据存入 SQLite，支持多区域组合搜索与 FASTA 批量分析 |

配套 IgBLAST 数据文件位于仓库 `database/igblast/`（`alpaca_gl.aux`、`alpaca.ndm.imgt`、`alpaca.ll` 等）。

## 依赖

```bash
pip install pandas
# 需要带 tkinter 的 Python 3（Linux: sudo apt install python3-tk）
```

## 使用方法

### 1. CDR 序列提取

```bash
python tools/mab-db/CDR-extract-V1.2.1-260618.py
```

1. 在界面中选择输入文件夹（递归扫描其中所有 `.xlsx` / `.csv`）与输出目录
2. 点击开始提取
3. 生成三个 FASTA 文件：`cdr1_sequences.fasta`、`cdr2_sequences.fasta`、`cdr3_sequences.fasta`

说明：

- 列名匹配只识别明确含 `aa` 的列（如 `cdr3_aa`、`CDR3.aa`），避免误匹配 `cdr3` 等简写
- 自动识别 CSV 编码与分隔符

### 2. CDR 搜索（SQLite）

```bash
python tools/mab-db/igblast-solo-SQLite260618-V1.3.2.py
```

1. **导入数据**：选择 `.xlsx` / `.csv` 文件导入，数据写入当前目录的 `cdr3_search.db`
   - 文件名以 `VHH` / `VH` / `VL` 开头时自动标注来源类型
   - 也支持手工添加单条记录
2. **搜索**：可按 CDR1 / CDR2 / CDR3 多区域组合检索
3. **FASTA 批量分析**：粘贴或加载 FASTA 批量查询
4. 结果详情会解析 `full_row_json` 展示原始行全部列

说明：数据库文件默认在运行目录下生成（`cdr3_search.db`），建议在固定目录运行以便复用。
