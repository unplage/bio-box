# 引物 Tm 值计算器 (Tm-calculate)

基于校准参数的 DNA 引物 Tm 值批量计算工具（tkinter 图形界面）。

## 文件说明

| 文件 | 说明 |
|:---|:---|
| `primer_Tm-calculator.py` | Tm / GC 含量 / 退火温度建议计算，支持单条与批量 |

校准基于 56 条序列，RMSE 1.16°C；默认条件 Na⁺ 50 mM、引物浓度 0.25 µM。

## 依赖

```bash
# 标准库即可，需带 tkinter 的 Python 3（Linux: sudo apt install python3-tk）
```

## 使用方法

```bash
python tools/Tm-calculate/primer_Tm-calculator.py
```

界面提供三种模式：

1. **单序列计算**
   - 手动输入或粘贴 FASTA 格式序列
   - 显示 Tm、GC 含量、退火温度建议等详情

2. **批量 FASTA 处理**
   - 上传 FASTA 文件或粘贴 FASTA 文本
   - 批量计算并导出结果

3. **批量序列处理（非 FASTA）**
   - 每行一条 DNA 序列（无标题行）
   - 批量计算并导出（支持 CSV）

结果可导出 CSV 保存。
