# 多序列比对工具 (msa)

纯前端多序列比对（MSA）：核酸 / 蛋白质 FASTA 比对、保守性着色、consensus、导出。数据不离开本机。

## 文件说明

| 文件 | 说明 |
|:---|:---|
| `msa.html` | 主界面（已内嵌 Worker 算法），浏览器直接打开 |
| `msa-worker.js` | 比对算法源码（构建时注入 `msa.html`，一般无需单独打开） |

## 依赖

无。浏览器直接打开即可（建议 Chrome / Edge / Firefox 较新版本，需支持 Web Worker）。

## 使用方法

```bash
# macOS
open tools/msa/msa.html
# Linux
xdg-open tools/msa/msa.html
# Windows
start tools\msa\msa.html
```

1. **上传序列**：拖拽或选择一个或多个 FASTA 文件（`.fa` / `.fasta` / `.seq` / `.txt`），自动合并记录；重名 ID 自动加 `_2` 后缀
2. **参数（可选，⚙️ 展开）**：
   - 序列类型：自动识别 / 手动指定核酸、蛋白质
   - 打分：自动（核酸=匹配分，蛋白=BLOSUM62）或 BLOSUM62 / PAM250 / 自定义匹配错配
   - gapOpen / gapExtend（代价 = gapOpen + (k−1)×gapExtend）
   - 是否惩罚两端缺口（默认自由末端）
   - 着色方式、consensus 保守阈值
3. **开始比对**：进度条显示两两比对 / 建树 / 渐进比对阶段
4. **查看结果**：左侧序列名 + consensus 行 + 保守度条 + 残基着色
5. **导出**：
   - 比对 FASTA（序列单行，不折行）
   - consensus FASTA：所有序列在每一列的多数残基拼成的一条「一致序列」，用于概括比对的共有序列
   - Newick 引导树（UPGMA）
   - 比对图片 PNG（canvas 截图下载）

## 算法说明

- 两两全局比对：Needleman–Wunsch，仿射缺口，Web Worker 计算
- 距离：1 − 比对一致率（p-distance）
- 引导树：UPGMA
- 合并：Clustal 风格渐进 profile 比对（列间平均打分）
- 适用规模：约 ≤50 条 × ≤500 bp/aa；更大规模请用 MUSCLE / MAFFT 等专业工具
