# His-mut 多点组氨酸突变引物生成器 (his-mut)

纯前端 HTML 工具：从测序 FASTA 序列中按标记序列定位多个突变点，生成组氨酸（His/CAC）突变引物并去重计数。

## 文件说明

| 文件 | 说明 |
|:---|:---|
| `Hismut.html` | `His-mut.py` 的纯前端 Web 版，浏览器直接打开即可使用 |

## 使用方法

```bash
# macOS
open tools/his-mut/Hismut.html
# Linux
xdg-open tools/his-mut/Hismut.html
# Windows
start tools\his-mut\Hismut.html
```

1. **输入序列**
   - 粘贴多序列 FASTA 到文本框，或点击「选择 .seq 文件」读取多个文件
   - 可点击「加载示例数据」快速试用
2. **配置突变点**
   - 点击「＋ 添加突变点」可配置多个，每个点在完整序列中独立定位
   - 每个突变点参数：标签、标记序列 motif、替换起始位、替换长度、突变碱基、上游臂长、下游臂长
   - 默认预填 L13 突变点（motif `AGCCCCTAAGCTCCTGATCTAT`，偏移 15→CAC，臂长 15/15），与原程序一致；可「重置为默认」
3. **生成**：点击「⚡ 生成引物」
   - 生成正向引物及反向互补引物，按原程序逻辑子串去重并计数（COUNT / INCLUDE）
4. **下载**：点击「📥 下载结果文件」，包含：
   - `primer_unique.txt`
   - `realign_primer_aa.txt`
   - `primer_unique.csv`
   - `de-primer_seq.txt`
   - `de-error.txt`
