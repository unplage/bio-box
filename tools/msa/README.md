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
2. **粘贴序列（可选）**：在拖拽区下方的「或直接粘贴 FASTA / 序列文本」文本框中粘贴，**输入/粘贴后自动解析**，无需额外按钮；支持带 `>` header 的 FASTA（多行序列自动合并），或每行一条裸序列（自动命名 seq1、seq2…）。可与文件上传混合使用，重名自动去重。数据不足 2 条时点「开始比对」会提示补充输入
3. **参数（可选，⚙️ 展开）**：
   - 序列类型：自动识别 / 手动指定核酸、蛋白质
   - 打分：自动（核酸=匹配分，蛋白=BLOSUM62）或 BLOSUM62 / PAM250 / 自定义匹配错配
   - gapOpen / gapExtend（代价 = gapOpen + (k−1)×gapExtend）
   - 是否惩罚两端缺口（默认自由末端）
   - 着色方式、consensus 保守阈值
4. **开始比对**：进度条显示两两比对 / 建树 / 渐进比对阶段
5. **查看结果**：支持「图形视图 / 文本视图」切换；**按引导树聚类顺序展示**——相近序列相邻（如抗体 VH/VL 分簇），不保留输入文件中的原始顺序，导出 FASTA 同样按此顺序
   - 图形：顶部位点标尺（每 10 列标注一次）+ 保守度条 + consensus 行 + 残基着色，gap 显示为短横线（导出 PNG 同样包含）
   - 文本：与图形视图一致的着色渲染（单行紧凑位点标尺：主刻度 `|` 紧跟位置数字、次刻度 `·`，规格贴近图形视图 16px 标尺；保守度条 + cons 行 + 残基颜色 + gap 短横线），等宽对齐，可鼠标选中，**滚动时标尺以不透明单行固定在顶部、与容器无缝隙**；或点「📋 复制比对文本」一键复制纯文本（含同样的刻度行）
   - **✂️ 编辑序列**：勾选参与比对的序列子集（面板始终列出本次会话的全部序列，可随时勾选/取消），点「🔁 重新比对所选」仅对勾选序列用当前参数重跑——适合剔除远缘、污染或无关序列后再看结果。未勾选的序列不会消失：以灰色原始序列显示在比对下方的「未参与比对」区，不计入共识 / 引导树 / 比对 FASTA 导出（复制文本与 PNG 包含该区）。清空或重新「开始比对」会重置会话
6. **导出**：
   - 比对 FASTA（序列单行，不折行）
   - consensus FASTA：所有序列在每一列的多数残基拼成的一条「一致序列」，用于概括比对的共有序列
   - Newick 引导树（Neighbor-Joining，含枝长）
   - 比对图片 PNG（canvas 截图下载）

## 算法说明

- 两两全局比对：Needleman–Wunsch，仿射缺口，Web Worker 计算
- 距离：1 − 比对一致率（p-distance）
- 引导树：Neighbor-Joining（NJ）
- 合并：Clustal 风格渐进 profile 比对（列间平均打分）
- 迭代精修：2 轮逐序列摘除重插（MUSCLE 风格），改善早期比对错误
- 定位：适合小规模、近缘序列的日常分析；远缘序列、结构域混搭或需要发表级质量时，请用 MAFFT / MUSCLE5 等专业工具
- 适用规模：约 ≤50 条 × ≤500 bp/aa

## 测试

零依赖回归测试（Node ≥ 16）：

```bash
node tools/msa/tests/run.js            # 全部（约 0.3 s）
node tools/msa/tests/run.js oracle     # 按名称过滤
MSA_WORKER=path/to/worker.js node tools/msa/tests/run.js   # 测试指定 worker 副本
```

覆盖范围：
- Gotoh 打分/回溯与**穷举 oracle 对拍**（小序列全路径枚举取最优，捕获回溯指针类 bug）
- 历史 pY 指针互换 bug 回归夹具（E51H2/h82H 对 + 抗体 6 序列全链路）
- NJ 拓扑、枝长与路径距离还原（可加矩阵）
- `runMSA` 不变量：无全 gap 列、列长一致、consensus 一致、newick 叶子齐全、确定性
- `detectType` / `makeScorer` 矩阵锁定

修改 `msa-worker.js` 后请先跑测试，再同步内嵌副本到 `msa.html`。
