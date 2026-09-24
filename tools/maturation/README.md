# 定向进化/成熟化引物工具 (maturation)

纯前端 HTML 工具集，浏览器直接打开即可使用，无需安装依赖。

## 文件说明

| 文件 | 说明 |
|:---|:---|
| `bsaitest.html` | BsaI 酶切位点重复分析：统计 DNA 序列中各 4-mer 在 BsaI 位点数据库中的出现频率 |
| `NNKprimer.html` | NNK 引物生成器：由编码区序列生成所有可能的 NNK 突变引物组合 |
| `2edmaturation.html` | 2ED Maturation 反向引物：截取 NNK 前 24 bp 计算 Tm 并生成反向互补引物 |

## 使用方法

### BsaI 酶切位点重复分析

```bash
# macOS
open tools/maturation/bsaitest.html
# Linux
xdg-open tools/maturation/bsaitest.html
```

1. 输入 DNA 序列（如 `ATGTGGTCAGCGATCGATCGATCG`）
2. 点击「分析」
3. 查看每个 4-mer 及其反向互补在数据库中的出现次数（含 0 次），支持下载结果

### NNK 引物生成器

```bash
open tools/maturation/NNKprimer.html
```

1. 填写参数：
   - **Upper Prefix**：上游前缀（如 `gccATG`）
   - **编码区序列**：长度须为 3 的倍数（如 `ATGTGGTCAGCG`）
   - **突变位点数**：如 `2`
   - **Later Suffix**：下游后缀（如 `tta`）
2. 点击「Generate」
3. 点击「下载 primer_list.txt」保存结果

### 2ED Maturation 反向引物

```bash
open tools/maturation/2edmaturation.html
```

1. 粘贴含 NNK 的 DNA 序列，每行一条
2. 点击「Process」
3. 下载 `primers.txt`（反向引物）与 `tm_values.txt`（Tm 数据）

计算说明：截取 NNK 前 24 个碱基，使用 SantaLucia 2004 热力学参数计算 Tm；校准基于 56 条序列（RMSE 1.16°C），默认 Na⁺ 50 mM、引物浓度 0.25 µM。
