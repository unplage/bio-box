# FASTA 格式转换工具 (format-convert)

纯前端 HTML 工具：多个 FASTA 文件与单个 FASTA 文件之间互转，浏览器本地运行，不上传数据。

## 文件说明

| 文件 | 说明 |
|:---|:---|
| `fasta-convert.html` | 合并（多 → 单）与拆分（单 → 多），支持 `.fa` / `.fasta` / `.seq` / `.txt` |

## 使用方法

```bash
# macOS
open tools/format-convert/fasta-convert.html
# Linux
xdg-open tools/format-convert/fasta-convert.html
# Windows
start tools\format-convert\fasta-convert.html
```

### 合并：多个文件 → 单个文件

1. 切换到「⤵️ 合并：多 → 单」标签
2. 拖拽或点击选择多个文件（可混选 `.fa` `.fasta` `.seq` `.txt`）
3. 可选配置：
   - **输出扩展名**：`.fa` / `.fasta` / `.seq` / `.txt`
   - **标题前添加文件名前缀**：避免不同文件中同名 ID 冲突（如 `file1|seq001`）
   - **ID 冲突时自动重编号**：重复 ID 追加 `_2`、`_3` …
4. 用 ↑↓ 调整合并顺序，✕ 移除个别文件
5. 点击「预览」查看结果，再点「⬇️ 合并并下载」

说明：无 `>` 标题的纯序列文件会以文件名（去掉扩展名）作为标题。

### 拆分：单个文件 → 多个文件

1. 切换到「⤴️ 拆分：单 → 多」标签
2. 拖拽或选择单个 FASTA 文件
3. 可选配置：
   - **输出扩展名**
   - **文件名取自**：序列 ID 或序号（`001`, `002`, …）
   - **保留完整标题**：关闭则只保留 ID
4. 下载方式二选一：
   - **逐个下载**：浏览器可能询问是否允许多文件下载，请选择允许
   - **打包为 ZIP 下载**：原生实现（不压缩），一次下载全部记录

说明：文件名会清洗非法字符；ID 重复时自动加 `_2`、`_3` 后缀。
