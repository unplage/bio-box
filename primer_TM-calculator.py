#V2版-增加非fasta格式批量处理
#!/usr/bin/env python3
#Tm批量计算-GUI版
"""
DNA引物Tm值计算器 - 图形界面版本
基于校准参数: 56条序列 | RMSE: 1.16°C

功能：
1. 单序列计算：手动输入或粘贴FASTA格式序列，显示Tm值、GC含量、退火温度建议等
2. 批量FASTA处理：上传FASTA文件或粘贴FASTA文本，批量计算并导出
3. 批量序列处理：每行一条DNA序列（无标题），批量计算并导出
"""

import math
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
from typing import Union, List, Dict, Tuple
import csv
import os

# ==================== 核心计算模块 ====================

CALIBRATION_PARAMS = {
    "salt_correction_factor": 1.6,
    "terminal_AT_penalty": 2.497083789430314,
    "entropy_scaling_factor": 0.9675869088668494,
    "initiation_adjustment": -0.0007103390016273732,
    "calibration_rmse": 1.1553994766937066,
    "calibration_points": 56
}

NN_DH = {
    'AA': -7.6, 'TT': -7.6, 'AT': -7.2, 'TA': -7.2,
    'CA': -8.5, 'TG': -8.5, 'GT': -8.4, 'AC': -8.4,
    'CT': -7.8, 'AG': -7.8, 'GA': -8.2, 'TC': -8.2,
    'CG': -10.6, 'GC': -10.6, 'GG': -8.0, 'CC': -8.0
}

NN_DS = {
    'AA': -21.3, 'TT': -21.3, 'AT': -20.4, 'TA': -20.4,
    'CA': -22.7, 'TG': -22.7, 'GT': -22.4, 'AC': -22.4,
    'CT': -21.0, 'AG': -21.0, 'GA': -22.2, 'TC': -22.2,
    'CG': -27.2, 'GC': -27.2, 'GG': -19.9, 'CC': -19.9
}


def calculate_tm(sequence: str, na_conc: float = 50.0,
                 primer_conc: float = 0.25, annealing_mode: bool = False,
                 return_details: bool = False) -> Union[float, Dict]:
    """计算单条DNA序列的Tm值"""
    if na_conc <= 0:
        raise ValueError(f"Na⁺浓度必须为正数，当前为 {na_conc} mM")
    if primer_conc <= 0:
        raise ValueError(f"引物浓度必须为正数，当前为 {primer_conc} μM")

    seq = sequence.upper().strip()
    n = len(seq)

    if n < 4:
        raise ValueError(f"序列太短 ({n} bp)，至少需要4bp")

    dH = 0.2 + CALIBRATION_PARAMS["initiation_adjustment"] * 0.05
    dS = -5.7 * CALIBRATION_PARAMS["entropy_scaling_factor"]

    for i in range(n - 1):
        nn = seq[i:i+2]
        dH += NN_DH[nn]
        dS += NN_DS[nn] * CALIBRATION_PARAMS["entropy_scaling_factor"]

    if seq[-1] in ('A', 'T'):
        dH += CALIBRATION_PARAMS["terminal_AT_penalty"]
        dS += 6.9

    monovalent_M = na_conc / 1000.0
    salt_effect = 0.368 * (n - 1) * math.log(monovalent_M) * CALIBRATION_PARAMS["salt_correction_factor"]
    dS_corrected = dS + salt_effect

    R = 1.987
    Ct = primer_conc * 1e-6

    if annealing_mode:
        sym_factor = 1
    else:
        sym_factor = 4

    denominator = dS_corrected + R * math.log(Ct / sym_factor)

    if denominator >= 0:
        gc_count = seq.count('G') + seq.count('C')
        tm = 2 * (n - gc_count) + 4 * gc_count
    else:
        Tm_K = (dH * 1000) / denominator
        tm = Tm_K - 273.15

    if not return_details:
        return round(tm, 2)

    gc_count = seq.count('G') + seq.count('C')
    at_count = n - gc_count

    return {
        "tm": round(tm, 2),
        "sequence": seq,
        "length": n,
        "gc_percent": round(gc_count / n * 100, 2),
        "gc_count": gc_count,
        "at_count": at_count,
        "conditions": {
            "na_conc": na_conc,
            "primer_conc": primer_conc,
            "annealing_mode": annealing_mode
        },
        "annealing_suggestions": {
            "standard_enzyme": [round(tm - 3, 1), round(tm - 1, 1)],
            "high_fidelity": [round(tm + 3, 1), round(tm + 5, 1)],
            "gradient_pcr": [round(tm - 5, 1), round(tm + 5, 1)]
        },
        "calibration_info": {
            "rmse": CALIBRATION_PARAMS["calibration_rmse"],
            "data_points": CALIBRATION_PARAMS["calibration_points"]
        }
    }


def calculate_tm_batch(sequences: List[str], na_conc: float = 50.0,
                       primer_conc: float = 0.25, annealing_mode: bool = False,
                       return_format: str = "list") -> Union[List[float], Dict]:
    results = []
    errors = []

    for i, seq in enumerate(sequences):
        try:
            if return_format == "dict":
                result = calculate_tm(seq, na_conc, primer_conc, annealing_mode, return_details=True)
                result["index"] = i
                results.append(result)
            else:
                tm = calculate_tm(seq, na_conc, primer_conc, annealing_mode, return_details=False)
                results.append(tm)
        except Exception as e:
            errors.append({
                "index": i,
                "sequence": seq[:50] + "..." if len(seq) > 50 else seq,
                "error": str(e)
            })

    if return_format == "dict":
        return {
            "total": len(sequences),
            "successful": len(results),
            "failed": len(errors),
            "results": results,
            "errors": errors
        }
    return results


def validate_sequence(sequence: str) -> Tuple[bool, str]:
    seq = sequence.upper().strip()
    if not seq:
        return False, "序列不能为空"
    if len(seq) < 4:
        return False, f"序列太短 ({len(seq)} bp)，至少需要4bp"
    if len(seq) > 5000:
        return False, f"序列太长 ({len(seq)} bp)，最大支持5000bp"
    invalid_chars = [c for c in seq if c not in "ATCG"]
    if invalid_chars:
        return False, f"序列包含非法字符: {', '.join(set(invalid_chars))}"
    return True, seq


def parse_fasta(text: str) -> List[Tuple[str, str]]:
    entries = []
    current_header = None
    current_seq = []

    for line in text.strip().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if current_header is not None:
                full_seq = ''.join(current_seq)
                if full_seq:
                    entries.append((current_header, full_seq))
            current_header = line[1:].strip()
            current_seq = []
        else:
            cleaned = ''.join(ch for ch in line.upper() if ch in 'ATCG')
            if not cleaned:
                continue
            if current_header is None:
                current_header = "unnamed"
                current_seq = [cleaned]
            else:
                current_seq.append(cleaned)

    if current_header is not None and current_seq:
        full_seq = ''.join(current_seq)
        if full_seq:
            entries.append((current_header, full_seq))

    if not entries:
        cleaned = ''.join(ch for ch in text.upper() if ch in 'ATCG')
        if cleaned:
            entries.append(("unnamed", cleaned))

    return entries


def parse_plain_sequences(text: str) -> List[str]:
    """解析纯文本序列：每行一条序列，去除空行和空格"""
    sequences = []
    for line in text.strip().splitlines():
        line = line.strip()
        if not line:
            continue
        # 只保留ATCG字符
        cleaned = ''.join(ch for ch in line.upper() if ch in 'ATCG')
        if cleaned:
            sequences.append(cleaned)
    return sequences


# ==================== GUI 应用程序 ====================

class TmCalculatorApp:
    def __init__(self, root):
        self.root = root
        self.root.title("DNA引物Tm值计算器 v1.0")
        self.root.geometry("1100x700")
        self.root.resizable(True, True)

        style = ttk.Style()
        style.theme_use('clam')

        self.na_conc = tk.DoubleVar(value=50.0)
        self.primer_conc = tk.DoubleVar(value=0.25)
        self.annealing_mode = tk.BooleanVar(value=False)

        self.create_menu()
        self.create_main_layout()

    def create_menu(self):
        menubar = tk.Menu(self.root)
        self.root.config(menu=menubar)
        file_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="文件", menu=file_menu)
        file_menu.add_command(label="退出", command=self.root.quit)
        help_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="帮助", menu=help_menu)
        help_menu.add_command(label="使用说明", command=self.show_help)

    def create_main_layout(self):
        param_frame = ttk.LabelFrame(self.root, text="计算参数", padding=10)
        param_frame.pack(fill=tk.X, padx=10, pady=5)

        row1 = ttk.Frame(param_frame)
        row1.pack(fill=tk.X, pady=2)

        ttk.Label(row1, text="Na⁺浓度 (mM):").pack(side=tk.LEFT, padx=5)
        na_entry = ttk.Entry(row1, textvariable=self.na_conc, width=10)
        na_entry.pack(side=tk.LEFT, padx=5)
        ttk.Label(row1, text=" 引物浓度 (μM):").pack(side=tk.LEFT, padx=5)
        primer_entry = ttk.Entry(row1, textvariable=self.primer_conc, width=10)
        primer_entry.pack(side=tk.LEFT, padx=5)

        cb = ttk.Checkbutton(row1, text="退火模式 (Annealing Mode)", variable=self.annealing_mode)
        cb.pack(side=tk.LEFT, padx=20)
        ttk.Label(row1, text="PCR模式(默认)使用对称因子4，退火模式使用对称因子1",
                  foreground="gray", font=("Arial", 8)).pack(side=tk.LEFT, padx=5)

        info_label = ttk.Label(row1, text=f"校准RMSE: {CALIBRATION_PARAMS['calibration_rmse']}°C (基于{CALIBRATION_PARAMS['calibration_points']}条序列)")
        info_label.pack(side=tk.RIGHT, padx=10)

        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # 三个选项卡
        self.single_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.single_frame, text="单序列计算")

        self.batch_fasta_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.batch_fasta_frame, text="批量FASTA处理")

        self.batch_plain_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.batch_plain_frame, text="批量序列处理")

        self.create_single_tab()
        self.create_batch_fasta_tab()
        self.create_batch_plain_tab()

    def create_single_tab(self):
        left_frame = ttk.Frame(self.single_frame)
        left_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=5, pady=5)

        ttk.Label(left_frame, text="DNA序列 (支持FASTA格式，将自动提取第一条序列):", font=("Arial", 10, "bold")).pack(anchor=tk.W)
        self.seq_text = scrolledtext.ScrolledText(left_frame, height=15, width=50, font=("Courier", 10))
        self.seq_text.pack(fill=tk.BOTH, expand=True, pady=5)

        btn_frame = ttk.Frame(left_frame)
        btn_frame.pack(fill=tk.X, pady=5)
        ttk.Button(btn_frame, text="计算Tm", command=self.calculate_single).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="清除", command=self.clear_single).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="粘贴示例", command=self.paste_example).pack(side=tk.LEFT, padx=5)

        right_frame = ttk.LabelFrame(self.single_frame, text="计算结果", padding=10)
        right_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=5, pady=5)

        self.result_text = scrolledtext.ScrolledText(right_frame, height=20, width=50, font=("Courier", 10), wrap=tk.WORD)
        self.result_text.pack(fill=tk.BOTH, expand=True)

    def create_batch_fasta_tab(self):
        # 顶部按钮
        top_frame = ttk.Frame(self.batch_fasta_frame)
        top_frame.pack(fill=tk.X, padx=5, pady=5)

        ttk.Button(top_frame, text="打开FASTA文件", command=self.load_fasta_file).pack(side=tk.LEFT, padx=5)
        ttk.Button(top_frame, text="从文本框解析", command=self.parse_batch_text).pack(side=tk.LEFT, padx=5)
        ttk.Button(top_frame, text="批量计算", command=self.calculate_batch_fasta).pack(side=tk.LEFT, padx=5)
        ttk.Button(top_frame, text="导出结果", command=self.export_batch_fasta_results).pack(side=tk.LEFT, padx=5)
        ttk.Button(top_frame, text="清除表格", command=self.clear_batch_fasta_table).pack(side=tk.LEFT, padx=5)

        # FASTA文本区域
        ttk.Label(self.batch_fasta_frame, text="FASTA文本 (可直接粘贴):", font=("Arial", 9, "bold")).pack(anchor=tk.W, padx=5)
        self.batch_fasta_text = scrolledtext.ScrolledText(self.batch_fasta_frame, height=8, font=("Courier", 9))
        self.batch_fasta_text.pack(fill=tk.X, padx=5, pady=2)

        # 结果表格
        table_frame = ttk.LabelFrame(self.batch_fasta_frame, text="批量计算结果", padding=5)
        table_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        columns = ("ID", "序列(前50bp)", "长度(bp)", "GC%", "Tm(°C)", "标准酶退火范围")
        self.tree_fasta = ttk.Treeview(table_frame, columns=columns, show="headings", height=12)
        for col in columns:
            self.tree_fasta.heading(col, text=col)
            if col == "ID":
                self.tree_fasta.column(col, width=150)
            elif col == "序列(前50bp)":
                self.tree_fasta.column(col, width=250)
            else:
                self.tree_fasta.column(col, width=100)

        vsb = ttk.Scrollbar(table_frame, orient="vertical", command=self.tree_fasta.yview)
        self.tree_fasta.configure(yscrollcommand=vsb.set)
        self.tree_fasta.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        vsb.pack(side=tk.RIGHT, fill=tk.Y)

        self.batch_fasta_results_data = []

    def create_batch_plain_tab(self):
        """新增：批量序列处理（每行一条序列）"""
        # 顶部按钮
        top_frame = ttk.Frame(self.batch_plain_frame)
        top_frame.pack(fill=tk.X, padx=5, pady=5)

        ttk.Button(top_frame, text="从文件加载", command=self.load_plain_file).pack(side=tk.LEFT, padx=5)
        ttk.Button(top_frame, text="批量计算", command=self.calculate_batch_plain).pack(side=tk.LEFT, padx=5)
        ttk.Button(top_frame, text="导出结果", command=self.export_batch_plain_results).pack(side=tk.LEFT, padx=5)
        ttk.Button(top_frame, text="清除表格", command=self.clear_batch_plain_table).pack(side=tk.LEFT, padx=5)

        # 序列文本区域
        ttk.Label(self.batch_plain_frame, text="序列列表 (每行一条DNA序列，仅含ATCG):", font=("Arial", 9, "bold")).pack(anchor=tk.W, padx=5)
        self.batch_plain_text = scrolledtext.ScrolledText(self.batch_plain_frame, height=8, font=("Courier", 9))
        self.batch_plain_text.pack(fill=tk.X, padx=5, pady=2)

        # 结果表格
        table_frame = ttk.LabelFrame(self.batch_plain_frame, text="批量计算结果", padding=5)
        table_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        columns = ("序号", "序列(前50bp)", "长度(bp)", "GC%", "Tm(°C)", "标准酶退火范围")
        self.tree_plain = ttk.Treeview(table_frame, columns=columns, show="headings", height=12)
        for col in columns:
            self.tree_plain.heading(col, text=col)
            if col == "序号":
                self.tree_plain.column(col, width=60)
            elif col == "序列(前50bp)":
                self.tree_plain.column(col, width=300)
            else:
                self.tree_plain.column(col, width=100)

        vsb = ttk.Scrollbar(table_frame, orient="vertical", command=self.tree_plain.yview)
        self.tree_plain.configure(yscrollcommand=vsb.set)
        self.tree_plain.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        vsb.pack(side=tk.RIGHT, fill=tk.Y)

        self.batch_plain_results_data = []

    # ------------------ 单序列计算 ------------------
    def calculate_single(self):
        raw_text = self.seq_text.get("1.0", tk.END).strip()
        if not raw_text:
            messagebox.showwarning("警告", "请输入或粘贴DNA序列")
            return

        entries = parse_fasta(raw_text)
        if not entries:
            messagebox.showerror("错误", "无法识别序列，请检查输入内容")
            return

        header, seq = entries[0]
        valid, msg = validate_sequence(seq)
        if not valid:
            messagebox.showerror("序列无效", msg)
            return

        try:
            na = self.na_conc.get()
            primer = self.primer_conc.get()
            anneal_mode = self.annealing_mode.get()

            details = calculate_tm(seq, na, primer, anneal_mode, return_details=True)

            result_str = f"""
═══════════════════════════════════════
  Tm计算结果 (基于最近邻热力学模型)
═══════════════════════════════════════

序列ID: {header}
序列长度: {details['length']} bp
GC含量: {details['gc_percent']}%
GC数/AT数: {details['gc_count']}/{details['at_count']}

计算条件:
  Na⁺浓度: {na} mM
  引物浓度: {primer} μM
  模式: {'退火模式' if anneal_mode else 'PCR模式'}

▶ Tm值: {details['tm']} °C

推荐退火温度 (标准Taq酶): {details['annealing_suggestions']['standard_enzyme'][0]} - {details['annealing_suggestions']['standard_enzyme'][1]} °C
推荐退火温度 (高保真酶):   {details['annealing_suggestions']['high_fidelity'][0]} - {details['annealing_suggestions']['high_fidelity'][1]} °C
梯度PCR建议范围:           {details['annealing_suggestions']['gradient_pcr'][0]} - {details['annealing_suggestions']['gradient_pcr'][1]} °C

校准信息:
  RMSE: {details['calibration_info']['rmse']}°C (基于{details['calibration_info']['data_points']}条实验序列)
═══════════════════════════════════════
"""
            self.result_text.delete("1.0", tk.END)
            self.result_text.insert("1.0", result_str)

        except Exception as e:
            messagebox.showerror("计算错误", f"计算Tm值时出错:\n{str(e)}")

    def clear_single(self):
        self.seq_text.delete("1.0", tk.END)
        self.result_text.delete("1.0", tk.END)

    def paste_example(self):
        example = """>Example_primer_F
ATCGATCGATCGATCG
>Example_primer_R
GCCAGTGCCAAGCTTGCA
"""
        self.seq_text.delete("1.0", tk.END)
        self.seq_text.insert("1.0", example)

    # ------------------ 批量FASTA处理 ------------------
    def load_fasta_file(self):
        filepath = filedialog.askopenfilename(
            title="选择FASTA文件",
            filetypes=[("FASTA文件", "*.fa *.fasta *.txt"), ("所有文件", "*.*")]
        )
        if not filepath:
            return
        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                content = f.read()
            self.batch_fasta_text.delete("1.0", tk.END)
            self.batch_fasta_text.insert("1.0", content)
            messagebox.showinfo("成功", f"已加载文件: {os.path.basename(filepath)}\n请点击「批量计算」按钮")
        except Exception as e:
            messagebox.showerror("错误", f"读取文件失败:\n{str(e)}")

    def parse_batch_text(self):
        text = self.batch_fasta_text.get("1.0", tk.END).strip()
        if not text:
            messagebox.showwarning("警告", "请先在文本框中粘贴FASTA格式序列")
            return

        entries = parse_fasta(text)
        if not entries:
            messagebox.showerror("错误", "未找到有效序列，请检查格式")
            return

        self.clear_batch_fasta_table()
        for header, seq in entries:
            short_seq = seq[:50] + "..." if len(seq) > 50 else seq
            self.tree_fasta.insert("", tk.END, values=(header, short_seq, len(seq), "待计算", "待计算", "待计算"))
        messagebox.showinfo("解析完成", f"成功解析 {len(entries)} 条序列，请点击「批量计算」")

    def clear_batch_fasta_table(self):
        for item in self.tree_fasta.get_children():
            self.tree_fasta.delete(item)
        # 不重置数据变量，避免影响导出（将在计算时重新赋值）

    def calculate_batch_fasta(self):
        text = self.batch_fasta_text.get("1.0", tk.END).strip()
        if not text:
            messagebox.showwarning("警告", "请先粘贴FASTA序列或加载文件")
            return

        entries = parse_fasta(text)
        if not entries:
            messagebox.showerror("错误", "未找到有效序列")
            return

        na = self.na_conc.get()
        primer = self.primer_conc.get()
        anneal_mode = self.annealing_mode.get()

        # 验证序列
        valid_entries = []
        invalid_details = []
        for header, seq in entries:
            valid, msg = validate_sequence(seq)
            if valid:
                valid_entries.append((header, seq))
            else:
                preview = seq[:50] + "..." if len(seq) > 50 else seq
                invalid_details.append((header, preview, msg))

        if not valid_entries:
            err_lines = ["所有序列均无效，原因如下："]
            for h, preview, err in invalid_details[:10]:
                err_lines.append(f"• {h}: {err} (片段: {preview})")
            messagebox.showerror("序列验证失败", "\n".join(err_lines))
            return

        if invalid_details:
            warn_lines = [f"以下 {len(invalid_details)} 条序列无效，将被跳过："]
            for h, preview, err in invalid_details[:5]:
                warn_lines.append(f"• {h}: {err}")
            messagebox.showwarning("部分序列无效", "\n".join(warn_lines))

        headers = [h for h, _ in valid_entries]
        sequences = [seq for _, seq in valid_entries]

        try:
            batch_result = calculate_tm_batch(sequences, na, primer, anneal_mode, return_format="dict")
            self.batch_fasta_results_data = batch_result["results"]
            self.clear_batch_fasta_table()

            if not self.batch_fasta_results_data:
                if batch_result.get("failed", 0) > 0:
                    err_msg = "所有序列计算失败:\n"
                    for err in batch_result["errors"][:10]:
                        err_msg += f"  {err['error']}\n"
                    messagebox.showerror("计算失败", err_msg)
                else:
                    messagebox.showwarning("无结果", "未产生任何有效结果")
                return

            for res in self.batch_fasta_results_data:
                idx = res["index"]
                header = headers[idx] if idx < len(headers) else f"Seq_{idx}"
                seq_short = res["sequence"][:50] + "..." if len(res["sequence"]) > 50 else res["sequence"]
                std_anneal = f"{res['annealing_suggestions']['standard_enzyme'][0]}-{res['annealing_suggestions']['standard_enzyme'][1]}"
                self.tree_fasta.insert("", tk.END, values=(
                    header, seq_short, res["length"],
                    f"{res['gc_percent']}%", f"{res['tm']}", std_anneal
                ))

            self.tree_fasta.update_idletasks()

            if batch_result.get("failed", 0) > 0:
                err_msg = f"成功: {batch_result['successful']}, 失败: {batch_result['failed']}\n"
                for err in batch_result["errors"]:
                    err_msg += f"  {err['error']}\n"
                messagebox.showwarning("部分失败", err_msg)
            else:
                messagebox.showinfo("计算完成", f"成功计算 {batch_result['successful']} 条序列")

        except Exception as e:
            messagebox.showerror("计算错误", str(e))

    def export_batch_fasta_results(self):
        if not self.batch_fasta_results_data:
            messagebox.showwarning("警告", "没有可导出的结果，请先进行批量计算")
            return
        filepath = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV文件", "*.csv")])
        if not filepath:
            return
        try:
            with open(filepath, 'w', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                writer.writerow(["ID", "序列(完整)", "长度(bp)", "GC%", "Tm(°C)",
                                 "标准酶退火范围", "高保真酶退火范围", "梯度PCR范围",
                                 "Na+浓度(mM)", "引物浓度(μM)", "计算模式"])
                for res in self.batch_fasta_results_data:
                    writer.writerow([
                        res.get("index", ""), res["sequence"], res["length"], res["gc_percent"], res["tm"],
                        f"{res['annealing_suggestions']['standard_enzyme'][0]}-{res['annealing_suggestions']['standard_enzyme'][1]}",
                        f"{res['annealing_suggestions']['high_fidelity'][0]}-{res['annealing_suggestions']['high_fidelity'][1]}",
                        f"{res['annealing_suggestions']['gradient_pcr'][0]}-{res['annealing_suggestions']['gradient_pcr'][1]}",
                        res["conditions"]["na_conc"], res["conditions"]["primer_conc"],
                        "退火模式" if res["conditions"]["annealing_mode"] else "PCR模式"
                    ])
            messagebox.showinfo("导出成功", f"已保存至 {filepath}")
        except Exception as e:
            messagebox.showerror("导出失败", str(e))

    # ------------------ 批量序列处理（每行一条） ------------------
    def load_plain_file(self):
        filepath = filedialog.askopenfilename(
            title="选择文本文件（每行一条序列）",
            filetypes=[("文本文件", "*.txt"), ("所有文件", "*.*")]
        )
        if not filepath:
            return
        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                content = f.read()
            self.batch_plain_text.delete("1.0", tk.END)
            self.batch_plain_text.insert("1.0", content)
            messagebox.showinfo("成功", f"已加载文件: {os.path.basename(filepath)}")
        except Exception as e:
            messagebox.showerror("错误", f"读取文件失败:\n{str(e)}")

    def clear_batch_plain_table(self):
        for item in self.tree_plain.get_children():
            self.tree_plain.delete(item)
        # 不重置数据变量

    def calculate_batch_plain(self):
        text = self.batch_plain_text.get("1.0", tk.END).strip()
        if not text:
            messagebox.showwarning("警告", "请在文本框中粘贴序列（每行一条）")
            return

        # 解析每行序列
        raw_sequences = parse_plain_sequences(text)
        if not raw_sequences:
            messagebox.showerror("错误", "未找到有效序列（至少需要4个ATCG字母）")
            return

        na = self.na_conc.get()
        primer = self.primer_conc.get()
        anneal_mode = self.annealing_mode.get()

        # 验证每条序列
        valid_seqs = []
        invalid_details = []
        for i, seq in enumerate(raw_sequences):
            valid, msg = validate_sequence(seq)
            if valid:
                valid_seqs.append(seq)
            else:
                preview = seq[:50] + "..." if len(seq) > 50 else seq
                invalid_details.append((i+1, preview, msg))

        if not valid_seqs:
            err_lines = ["所有序列均无效："]
            for idx, preview, err in invalid_details[:10]:
                err_lines.append(f"行{idx}: {err} (片段: {preview})")
            messagebox.showerror("序列验证失败", "\n".join(err_lines))
            return

        if invalid_details:
            warn_lines = [f"以下 {len(invalid_details)} 条序列无效，将被跳过："]
            for idx, preview, err in invalid_details[:5]:
                warn_lines.append(f"行{idx}: {err}")
            messagebox.showwarning("部分序列无效", "\n".join(warn_lines))

        try:
            batch_result = calculate_tm_batch(valid_seqs, na, primer, anneal_mode, return_format="dict")
            self.batch_plain_results_data = batch_result["results"]
            self.clear_batch_plain_table()

            if not self.batch_plain_results_data:
                if batch_result.get("failed", 0) > 0:
                    err_msg = "所有序列计算失败:\n"
                    for err in batch_result["errors"][:10]:
                        err_msg += f"  {err['error']}\n"
                    messagebox.showerror("计算失败", err_msg)
                else:
                    messagebox.showwarning("无结果", "未产生任何有效结果")
                return

            for res in self.batch_plain_results_data:
                seq_short = res["sequence"][:50] + "..." if len(res["sequence"]) > 50 else res["sequence"]
                std_anneal = f"{res['annealing_suggestions']['standard_enzyme'][0]}-{res['annealing_suggestions']['standard_enzyme'][1]}"
                self.tree_plain.insert("", tk.END, values=(
                    res["index"] + 1,  # 序号从1开始
                    seq_short,
                    res["length"],
                    f"{res['gc_percent']}%",
                    f"{res['tm']}",
                    std_anneal
                ))

            self.tree_plain.update_idletasks()

            if batch_result.get("failed", 0) > 0:
                err_msg = f"成功: {batch_result['successful']}, 失败: {batch_result['failed']}\n"
                for err in batch_result["errors"]:
                    err_msg += f"  {err['error']}\n"
                messagebox.showwarning("部分失败", err_msg)
            else:
                messagebox.showinfo("计算完成", f"成功计算 {batch_result['successful']} 条序列")

        except Exception as e:
            messagebox.showerror("计算错误", str(e))

    def export_batch_plain_results(self):
        if not self.batch_plain_results_data:
            messagebox.showwarning("警告", "没有可导出的结果，请先进行批量计算")
            return
        filepath = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV文件", "*.csv")])
        if not filepath:
            return
        try:
            with open(filepath, 'w', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                writer.writerow(["序号", "序列(完整)", "长度(bp)", "GC%", "Tm(°C)",
                                 "标准酶退火范围", "高保真酶退火范围", "梯度PCR范围",
                                 "Na+浓度(mM)", "引物浓度(μM)", "计算模式"])
                for res in self.batch_plain_results_data:
                    writer.writerow([
                        res["index"] + 1,
                        res["sequence"],
                        res["length"],
                        res["gc_percent"],
                        res["tm"],
                        f"{res['annealing_suggestions']['standard_enzyme'][0]}-{res['annealing_suggestions']['standard_enzyme'][1]}",
                        f"{res['annealing_suggestions']['high_fidelity'][0]}-{res['annealing_suggestions']['high_fidelity'][1]}",
                        f"{res['annealing_suggestions']['gradient_pcr'][0]}-{res['annealing_suggestions']['gradient_pcr'][1]}",
                        res["conditions"]["na_conc"],
                        res["conditions"]["primer_conc"],
                        "退火模式" if res["conditions"]["annealing_mode"] else "PCR模式"
                    ])
            messagebox.showinfo("导出成功", f"已保存至 {filepath}")
        except Exception as e:
            messagebox.showerror("导出失败", str(e))

    def show_help(self):
        help_text = """使用说明：

1. 单序列计算：
   - 输入或粘贴FASTA格式序列（自动取第一条），显示详细Tm结果和退火建议。

2. 批量FASTA处理：
   - 支持FASTA格式（>开头），可打开文件或直接粘贴文本。
   - 点击“批量计算”后显示表格，可导出CSV。

3. 批量序列处理（新功能）：
   - 每行一条DNA序列，无需标题行。
   - 支持从文本文件加载或直接粘贴。
   - 自动过滤非ATCG字符，计算Tm并显示结果。

参数说明：
- Na⁺浓度 (mM)：默认50 mM
- 引物浓度 (μM)：默认0.25 μM
- 退火模式：PCR模式（对称因子4）或退火模式（对称因子1）

校准精度：RMSE = 1.16°C (基于56条序列)
"""
        messagebox.showinfo("帮助", help_text)


def main():
    root = tk.Tk()
    app = TmCalculatorApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
