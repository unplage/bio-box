# -*- coding: utf-8 -*-
"""
CDR序列提取工具 V1.2.1
- 修正列名匹配逻辑：只匹配明确包含 'aa' 的列名（如 cdr3_aa, CDR3.aa），避免误匹配 cdr3 等简写列
- 支持 .xlsx 和 .csv 文件（自动识别编码/分隔符）
- 增加更多错误处理
"""

import os
import sys
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
from pathlib import Path
import pandas as pd
import traceback
import threading
import queue
import csv

# ==================== 核心提取函数 ====================
def extract_cdrs_from_folder(folder_path, output_dir, log_callback=None):
    """
    扫描文件夹，提取所有匹配的 Excel/CSV 文件的 CDR 序列，生成 FASTA 文件
    返回 (成功计数, 失败计数, 错误信息列表)
    """
    root_path = Path(folder_path)
    if not root_path.exists():
        return 0, 0, [f"目录不存在: {root_path}"]

    # 创建输出目录
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # 准备三个 FASTA 文件，键名使用 cdr1_aa, cdr2_aa, cdr3_aa
    fasta_files = {
        'cdr1_aa': out_dir / 'cdr1_sequences.fasta',
        'cdr2_aa': out_dir / 'cdr2_sequences.fasta',
        'cdr3_aa': out_dir / 'cdr3_sequences.fasta'
    }
    for f in fasta_files.values():
        f.write_text('', encoding='utf-8')

    # 打开文件句柄（追加模式）
    handles = {key: open(path, 'a', encoding='utf-8') for key, path in fasta_files.items()}

    total_success = 0
    total_fail = 0
    errors = []

    # 同时搜索 xlsx 和 csv 文件
    all_files = list(root_path.rglob("*.xlsx")) + list(root_path.rglob("*.csv"))
    if not all_files:
        errors.append("未找到任何 .xlsx 或 .csv 文件")
        for h in handles.values():
            h.close()
        return 0, 0, errors

    # 遍历文件
    for file_idx, file_path in enumerate(all_files):
        if log_callback:
            log_callback(f"处理文件: {file_path.name} ({file_idx+1}/{len(all_files)})")

        # 确定前缀（文件名开头）
        base_name = file_path.name
        prefix = None
        if base_name.startswith("VHH"):
            prefix = "VHH"
        elif base_name.startswith("VH"):
            prefix = "VH"
        elif base_name.startswith("VL"):
            prefix = "VL"
        else:
            errors.append(f"跳过不支持的文件: {file_path.name}（前缀应为 VHH/VH/VL）")
            total_fail += 1
            continue

        # 根据扩展名读取文件
        ext = file_path.suffix.lower()
        try:
            if ext == '.xlsx':
                df = pd.read_excel(file_path)
            else:  # .csv
                # 尝试多种编码
                try:
                    df = pd.read_csv(file_path, encoding='utf-8')
                except UnicodeDecodeError:
                    try:
                        df = pd.read_csv(file_path, encoding='gbk')
                    except UnicodeDecodeError:
                        df = pd.read_csv(file_path, encoding='latin-1')
                # 如果只有一列且包含分隔符，自动检测分隔符
                if len(df.columns) == 1 and (',' in df.iloc[0, 0] or '\t' in df.iloc[0, 0]):
                    with open(file_path, 'r', encoding='utf-8') as f:
                        sample = f.read(1024)
                        sniffer = csv.Sniffer()
                        delimiter = sniffer.sniff(sample).delimiter
                    df = pd.read_csv(file_path, encoding='utf-8', sep=delimiter)
        except Exception as e:
            errors.append(f"读取失败 {file_path.name}: {e}")
            total_fail += 1
            continue

        if df.empty:
            errors.append(f"文件为空: {file_path.name}")
            total_fail += 1
            continue

        # 第一列作为原始 ID
        first_col = df.columns[0]
        df.rename(columns={first_col: 'original_id'}, inplace=True)

        # 识别各 CDR 列 - 只匹配明确包含 'aa' 的列名，避免误匹配 cdr1/cdr2/cdr3 等简写
        cdr_cols = {}
        cdr_map = {
            'cdr1_aa': ['cdr1_aa', 'CDR1_aa', 'CDR1.aa', 'cdr1aa'],
            'cdr2_aa': ['cdr2_aa', 'CDR2_aa', 'CDR2.aa', 'cdr2aa'],
            'cdr3_aa': ['cdr3_aa', 'CDR3_aa', 'CDR3.aa', 'cdr3aa']
        }
        for target, aliases in cdr_map.items():
            found = None
            for col in df.columns:
                if col.lower() in [a.lower() for a in aliases]:
                    found = col
                    break
            if found:
                cdr_cols[target] = found
            else:
                errors.append(f"文件 {file_path.name} 缺少列: {target}，该CDR将跳过")

        if not cdr_cols:
            errors.append(f"文件 {file_path.name} 没有任何CDR列，跳过")
            total_fail += 1
            continue

        # 遍历每一行
        file_success = 0
        for idx, row in df.iterrows():
            orig_id = str(row.get('original_id', f'row_{idx}'))
            new_id = f"{orig_id}_{prefix}"

            for cdr_type, col_name in cdr_cols.items():
                seq = row.get(col_name)
                if pd.isna(seq) or not seq:
                    continue
                seq = str(seq).strip()
                if len(seq) == 0:
                    continue
                handle = handles[cdr_type]
                handle.write(f">{new_id}\n{seq}\n")
                file_success += 1

        total_success += file_success
        if log_callback:
            log_callback(f"  从 {file_path.name} 提取 {file_success} 条序列")

    # 关闭所有句柄
    for h in handles.values():
        h.close()

    if log_callback:
        log_callback("提取完成")
    return total_success, total_fail, errors


# ==================== GUI 界面 ====================
class CDRExtractorGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("CDR序列提取工具 V1.2.1")
        self.root.geometry("700x600")
        self.root.minsize(600, 500)

        self.input_folder = tk.StringVar(value="")
        self.output_folder = tk.StringVar(value="./Extracted_FASTA")
        self.log_queue = queue.Queue()

        self.setup_ui()
        self.process_log_queue()

    def setup_ui(self):
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.pack(fill=tk.BOTH, expand=True)

        ttk.Label(main_frame, text="CDR序列提取工具 V1.2.1（仅匹配含 _aa 的CDR列）",
                  font=("Microsoft YaHei", 16, "bold")).pack(pady=(0,10))

        # 输入文件夹
        input_frame = ttk.LabelFrame(main_frame, text="输入文件夹（含子文件夹）", padding="8")
        input_frame.pack(fill=tk.X, pady=5)
        row1 = ttk.Frame(input_frame)
        row1.pack(fill=tk.X, pady=2)
        ttk.Label(row1, text="文件夹:").pack(side=tk.LEFT, padx=5)
        entry_in = ttk.Entry(row1, textvariable=self.input_folder, width=50)
        entry_in.pack(side=tk.LEFT, padx=5)
        ttk.Button(row1, text="浏览", command=self.browse_input).pack(side=tk.LEFT, padx=5)

        # 输出文件夹
        output_frame = ttk.LabelFrame(main_frame, text="输出文件夹", padding="8")
        output_frame.pack(fill=tk.X, pady=5)
        row2 = ttk.Frame(output_frame)
        row2.pack(fill=tk.X, pady=2)
        ttk.Label(row2, text="文件夹:").pack(side=tk.LEFT, padx=5)
        entry_out = ttk.Entry(row2, textvariable=self.output_folder, width=50)
        entry_out.pack(side=tk.LEFT, padx=5)
        ttk.Button(row2, text="浏览", command=self.browse_output).pack(side=tk.LEFT, padx=5)

        action_frame = ttk.Frame(main_frame)
        action_frame.pack(fill=tk.X, pady=10)
        self.start_btn = ttk.Button(action_frame, text="开始提取", command=self.start_extract, width=15)
        self.start_btn.pack(side=tk.LEFT, padx=5)
        self.progress_bar = ttk.Progressbar(action_frame, mode='indeterminate', length=200)
        self.progress_bar.pack(side=tk.LEFT, padx=20)

        log_frame = ttk.LabelFrame(main_frame, text="运行日志", padding="8")
        log_frame.pack(fill=tk.BOTH, expand=True, pady=5)
        self.log_text = scrolledtext.ScrolledText(log_frame, height=15, font=("Consolas", 9), wrap=tk.WORD)
        self.log_text.pack(fill=tk.BOTH, expand=True)
        self.log_text.config(state=tk.DISABLED)

        status_frame = ttk.Frame(main_frame)
        status_frame.pack(fill=tk.X, pady=5)
        self.status_label = ttk.Label(status_frame, text="就绪", foreground="#7f8c8d")
        self.status_label.pack(side=tk.LEFT)
        ttk.Button(status_frame, text="清空日志", command=self.clear_log).pack(side=tk.RIGHT, padx=5)

    def browse_input(self):
        folder = filedialog.askdirectory()
        if folder:
            self.input_folder.set(folder)

    def browse_output(self):
        folder = filedialog.askdirectory()
        if folder:
            self.output_folder.set(folder)

    def log(self, msg):
        self.log_queue.put(msg)
        self.process_log_queue()

    def clear_log(self):
        self.log_text.config(state=tk.NORMAL)
        self.log_text.delete(1.0, tk.END)
        self.log_text.config(state=tk.DISABLED)

    def process_log_queue(self):
        try:
            while True:
                msg = self.log_queue.get_nowait()
                self.log_text.config(state=tk.NORMAL)
                self.log_text.insert(tk.END, msg + "\n")
                self.log_text.see(tk.END)
                self.log_text.config(state=tk.DISABLED)
        except queue.Empty:
            pass
        self.root.after(100, self.process_log_queue)

    def start_extract(self):
        input_path = self.input_folder.get().strip()
        output_path = self.output_folder.get().strip()
        if not input_path:
            messagebox.showwarning("警告", "请选择输入文件夹")
            return
        if not output_path:
            output_path = "./Extracted_FASTA"

        self.start_btn.config(state=tk.DISABLED, text="提取中...")
        self.progress_bar.start(10)
        self.status_label.config(text="正在提取...")
        self.log("开始提取CDR序列...")

        def worker():
            try:
                success, fail, errors = extract_cdrs_from_folder(
                    input_path,
                    output_path,
                    log_callback=self.log
                )
                self.root.after(0, lambda: self.finish_extract(success, fail, errors))
            except Exception as e:
                err_msg = traceback.format_exc()
                self.log(f"发生错误:\n{err_msg}")
                self.root.after(0, lambda: self.finish_extract(0, 0, [str(e)]))

        threading.Thread(target=worker, daemon=True).start()

    def finish_extract(self, success, fail, errors):
        self.progress_bar.stop()
        self.start_btn.config(state=tk.NORMAL, text="开始提取")
        self.status_label.config(text="提取完成")
        msg = f"提取完成: 成功 {success} 条序列"
        if fail > 0:
            msg += f", 失败 {fail} 个文件"
        self.log(msg)
        if errors:
            self.log("--- 错误/警告信息 ---")
            for err in errors[:20]:
                self.log(f"  {err}")
            if len(errors) > 20:
                self.log(f"  ... 还有 {len(errors)-20} 条信息")
        messagebox.showinfo("完成", msg)


def main():
    root = tk.Tk()
    app = CDRExtractorGUI(root)
    root.mainloop()

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print("程序异常:", str(e))
        traceback.print_exc()
        input("按回车退出...")