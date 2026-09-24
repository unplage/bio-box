# -*- coding: utf-8 -*-
"""
CDR3/CDR2/CDR1序列搜索工具 V1.3.2
- 导入支持 .xlsx 和 .csv（自动识别编码/分隔符）
- 修正 full_row_json 为单行JSON，大幅减小数据库体积
- 详情显示完整原始列（full_row_json 解析）
- 列名匹配仅识别含 'aa' 的CDR列
- 手工添加、FASTA批量分析、多区域搜索
"""

import os
import sys
import sqlite3
import json
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
from pathlib import Path
import pandas as pd
import traceback
import csv

# ==================== 配置 ====================
DEFAULT_DIR = Path("./Analyse_file")
DB_PATH = "cdr3_search.db"

def get_file_prefix(filename):
    """根据文件名返回前缀（VHH/VH/VL），否则返回None"""
    base = os.path.basename(filename)
    if base.startswith("VHH"):
        return "VHH"
    elif base.startswith("VH"):
        return "VH"
    elif base.startswith("VL"):
        return "VL"
    else:
        return None

# ==================== 安全显示辅助 ====================
def safe_str(value):
    return '' if value is None else str(value)

# ==================== 数据库操作 ====================
def init_database():
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS cdr3_data (
            sequence_id TEXT PRIMARY KEY,
            source_type TEXT,
            original_id TEXT,
            cdr3_aa TEXT,
            v_call TEXT,
            d_call TEXT,
            j_call TEXT,
            junction_aa TEXT,
            junction_aa_length INTEGER,
            v_identity REAL,
            d_identity REAL,
            j_identity REAL,
            productive TEXT,
            sequence_aa TEXT,
            cdr1_aa TEXT,
            cdr2_aa TEXT,
            full_row_json TEXT,
            created_time TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    ''')
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_cdr3 ON cdr3_data(cdr3_aa)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_cdr2 ON cdr3_data(cdr2_aa)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_cdr1 ON cdr3_data(cdr1_aa)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_seqid ON cdr3_data(sequence_id)")
    conn.commit()
    conn.close()

def import_excel_files(folder_path):
    """
    扫描文件夹（含子文件夹）中的所有 .xlsx 和 .csv 文件，
    导入数据（追加，按 sequence_id 去重）
    返回: (新增总数, 错误信息列表)
    """
    init_database()
    root_path = Path(folder_path)
    if not root_path.exists():
        return 0, f"目录不存在: {root_path}"

    all_files = list(root_path.rglob("*.xlsx")) + list(root_path.rglob("*.csv"))
    if not all_files:
        return 0, "未找到任何 .xlsx 或 .csv 文件"

    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()
    total_new = 0
    errors = []

    for file_path in all_files:
        prefix = get_file_prefix(file_path.name)
        if not prefix:
            errors.append(f"跳过不支持的文件: {file_path.name}")
            continue

        # 读取文件
        ext = file_path.suffix.lower()
        try:
            if ext == '.xlsx':
                df = pd.read_excel(file_path)
            else:  # .csv
                try:
                    df = pd.read_csv(file_path, encoding='utf-8')
                except UnicodeDecodeError:
                    try:
                        df = pd.read_csv(file_path, encoding='gbk')
                    except UnicodeDecodeError:
                        df = pd.read_csv(file_path, encoding='latin-1')
                # 如果只有一列，自动检测分隔符
                if len(df.columns) == 1 and (',' in df.iloc[0, 0] or '\t' in df.iloc[0, 0]):
                    with open(file_path, 'r', encoding='utf-8') as f:
                        sample = f.read(1024)
                        sniffer = csv.Sniffer()
                        delimiter = sniffer.sniff(sample).delimiter
                    df = pd.read_csv(file_path, encoding='utf-8', sep=delimiter)
        except Exception as e:
            errors.append(f"读取失败 {file_path.name}: {e}")
            continue

        if df.empty:
            errors.append(f"文件为空: {file_path.name}")
            continue

        # 第一列作为原始ID
        df.rename(columns={df.columns[0]: 'original_id'}, inplace=True)

        # 查找 cdr3_aa 列（仅匹配含 'aa' 的列）
        if 'cdr3_aa' not in df.columns:
            possible = [c for c in df.columns if 'cdr3' in c.lower() and 'aa' in c.lower()]
            if possible:
                df.rename(columns={possible[0]: 'cdr3_aa'}, inplace=True)
            else:
                errors.append(f"文件 {file_path.name} 缺少CDR3列，跳过")
                continue

        # 其他列映射（仅匹配含 'aa' 的列名变体）
        col_map = {
            'sequence_aa': ['sequence_aa', 'seq_aa', 'Sequence_aa', 'Sequence AA'],
            'cdr1_aa': ['cdr1_aa', 'CDR1_aa', 'CDR1.aa', 'cdr1aa'],
            'cdr2_aa': ['cdr2_aa', 'CDR2_aa', 'CDR2.aa', 'cdr2aa'],
            'v_call': ['v_call', 'V_call'],
            'd_call': ['d_call', 'D_call'],
            'j_call': ['j_call', 'J_call'],
            'junction_aa': ['junction_aa', 'Junction_aa'],
            'junction_aa_length': ['junction_aa_length', 'junction_length'],
            'v_identity': ['v_identity', 'V_identity'],
            'd_identity': ['d_identity', 'D_identity'],
            'j_identity': ['j_identity', 'J_identity'],
            'productive': ['productive', 'Productive']
        }
        for target, candidates in col_map.items():
            if target not in df.columns:
                found = None
                for cand in candidates:
                    if cand in df.columns:
                        found = cand
                        break
                if found:
                    df.rename(columns={found: target}, inplace=True)
                else:
                    df[target] = ''

        # 生成新ID
        df['sequence_id'] = df['original_id'].astype(str) + '_' + prefix
        df['source_type'] = prefix

        # 选择需要的列
        base_cols = ['sequence_id', 'source_type', 'original_id', 'cdr3_aa', 'v_call', 'd_call', 'j_call',
                     'junction_aa', 'junction_aa_length', 'v_identity', 'd_identity', 'j_identity',
                     'productive', 'sequence_aa', 'cdr1_aa', 'cdr2_aa']
        cols = [c for c in base_cols if c in df.columns]
        df_sub = df[cols].copy()

        # ===== 修正：每行存储自身的 JSON（仅该行） =====
        df_sub['full_row_json'] = df.apply(lambda row: row.to_json(force_ascii=False), axis=1)
        # ==============================================

        # 获取现有ID集合
        existing_ids = set()
        cursor.execute("SELECT sequence_id FROM cdr3_data")
        for row in cursor.fetchall():
            existing_ids.add(row[0])

        # 过滤新记录
        df_new = df_sub[~df_sub['sequence_id'].isin(existing_ids)]
        new_count = len(df_new)
        if new_count > 0:
            df_new.to_sql('cdr3_data', conn, if_exists='append', index=False, chunksize=500)
            total_new += new_count

    conn.commit()
    return total_new, "\n".join(errors) if errors else None

def add_manual_record(seq_id, cdr3_aa, cdr1_aa='', cdr2_aa='', v_call='', d_call='', j_call='',
                     junction_aa='', junction_aa_length='', v_identity='', d_identity='',
                     j_identity='', productive='', sequence_aa='', source_type='manual'):
    init_database()
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()
    try:
        cursor.execute("SELECT 1 FROM cdr3_data WHERE sequence_id = ?", (seq_id,))
        exists = cursor.fetchone() is not None

        fields = ['sequence_id', 'source_type', 'cdr3_aa', 'cdr1_aa', 'cdr2_aa', 'v_call', 'd_call', 'j_call',
                  'junction_aa', 'junction_aa_length', 'v_identity', 'd_identity', 'j_identity', 'productive', 'sequence_aa']
        values = [seq_id, source_type, cdr3_aa, cdr1_aa, cdr2_aa, v_call, d_call, j_call,
                  junction_aa, junction_aa_length, v_identity, d_identity, j_identity, productive, sequence_aa]
        values = [None if v == '' else v for v in values]

        if exists:
            set_clause = ", ".join([f"{f} = ?" for f in fields[1:]])
            sql = f"UPDATE cdr3_data SET {set_clause} WHERE sequence_id = ?"
            cursor.execute(sql, values[1:] + [seq_id])
            conn.commit()
            return True, f"记录已更新: {seq_id}"
        else:
            placeholders = ", ".join(["?"] * len(fields))
            sql = f"INSERT INTO cdr3_data ({', '.join(fields)}) VALUES ({placeholders})"
            cursor.execute(sql, values)
            conn.commit()
            return True, f"记录已添加: {seq_id}"
    except Exception as e:
        return False, str(e)
    finally:
        conn.close()

def search_region(region, seq, exact_match=True):
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()
    if exact_match:
        sql = f"SELECT * FROM cdr3_data WHERE {region} = ? ORDER BY sequence_id"
        params = (seq,)
    else:
        sql = f"SELECT * FROM cdr3_data WHERE {region} LIKE ? ORDER BY sequence_id"
        params = (f'%{seq}%',)
    cursor.execute(sql, params)
    rows = cursor.fetchall()
    results = [dict(row) for row in rows]
    conn.close()
    return results

def get_db_stats():
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()
    cursor.execute("SELECT COUNT(*) FROM cdr3_data")
    count = cursor.fetchone()[0]
    conn.close()
    return count

def parse_fasta(filepath):
    seqs = {}
    current_id = None
    current_seq = []
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_id is not None:
                    seqs[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id is not None:
            seqs[current_id] = ''.join(current_seq)
    return seqs

def analyze_fasta_file(filepath):
    seqs = parse_fasta(filepath)
    results = []
    for seq_id, seq in seqs.items():
        seq_upper = seq.upper()
        matches = search_region('cdr3_aa', seq_upper, exact_match=True)
        exists = len(matches) > 0
        matched_ids = [m['sequence_id'] for m in matches] if exists else []
        results.append({
            'seq_id': seq_id,
            'cdr3_seq': seq_upper,
            'exists': exists,
            'matched_ids': ', '.join(matched_ids)
        })
    return results

# ==================== GUI ====================
class CDR3SearchApp:
    def __init__(self, root):
        self.root = root
        self.root.title("CDR3/CDR2/CDR1 搜索工具 V1.3.2")
        self.root.geometry("1150x850")
        self.root.minsize(1000, 750)

        self.notebook = ttk.Notebook(root)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        self.create_search_tab()
        self.create_management_tab()

        try:
            init_database()
            self.update_stats()
        except Exception as e:
            messagebox.showerror("错误", f"数据库初始化失败: {e}")

    # ================== 搜索标签页 ==================
    def create_search_tab(self):
        self.search_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.search_tab, text="序列搜索")

        tool_frame = ttk.Frame(self.search_tab)
        tool_frame.pack(fill=tk.X, pady=5)
        self.stats_label = ttk.Label(tool_frame, text="数据库: 0 条记录", foreground="#2ecc71")
        self.stats_label.pack(side=tk.LEFT, padx=10)
        ttk.Button(tool_frame, text="刷新统计", command=self.update_stats).pack(side=tk.LEFT, padx=5)
        ttk.Button(tool_frame, text="清空数据库", command=self.clear_db).pack(side=tk.LEFT, padx=5)
        ttk.Button(tool_frame, text="导出全部", command=self.export_all).pack(side=tk.LEFT, padx=5)

        search_frame = ttk.LabelFrame(self.search_tab, text="搜索条件", padding="10")
        search_frame.pack(fill=tk.X, pady=10, padx=5)

        row0 = ttk.Frame(search_frame)
        row0.pack(fill=tk.X, pady=2)
        ttk.Label(row0, text="搜索区域:").pack(side=tk.LEFT, padx=5)
        self.region_var = tk.StringVar(value="cdr3_aa")
        region_combo = ttk.Combobox(row0, textvariable=self.region_var,
                                    values=["cdr3_aa", "cdr2_aa", "cdr1_aa"], state="readonly", width=10)
        region_combo.pack(side=tk.LEFT, padx=5)
        ttk.Label(row0, text="序列:").pack(side=tk.LEFT, padx=5)
        self.search_entry = ttk.Entry(row0, width=40, font=("Courier", 10))
        self.search_entry.pack(side=tk.LEFT, padx=5)
        self.search_entry.bind("<Return>", lambda e: self.do_search())

        row1 = ttk.Frame(search_frame)
        row1.pack(fill=tk.X, pady=5)
        self.exact_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(row1, text="精确匹配", variable=self.exact_var).pack(side=tk.LEFT, padx=5)
        ttk.Button(row1, text="搜索", command=self.do_search, width=10).pack(side=tk.LEFT, padx=5)
        ttk.Button(row1, text="清空结果", command=self.clear_results, width=10).pack(side=tk.LEFT, padx=5)
        ttk.Button(row1, text="导出搜索结果", command=self.export_search_results, width=15).pack(side=tk.LEFT, padx=5)

        result_paned = ttk.PanedWindow(self.search_tab, orient=tk.VERTICAL)
        result_paned.pack(fill=tk.BOTH, expand=True, pady=5, padx=5)

        table_frame = ttk.LabelFrame(result_paned, text="搜索结果", padding="5")
        result_paned.add(table_frame, weight=3)
        columns = ("sequence_id", "cdr1_aa", "cdr2_aa", "cdr3_aa", "v_call", "j_call", "sequence_aa")
        self.search_tree = ttk.Treeview(table_frame, columns=columns, show="headings", height=12)
        widths = {"sequence_id": 180, "cdr1_aa": 100, "cdr2_aa": 100, "cdr3_aa": 120,
                  "v_call": 120, "j_call": 80, "sequence_aa": 250}
        for col in columns:
            self.search_tree.heading(col, text=col)
            self.search_tree.column(col, width=widths.get(col, 100))
        scroll_y = ttk.Scrollbar(table_frame, orient=tk.VERTICAL, command=self.search_tree.yview)
        scroll_x = ttk.Scrollbar(table_frame, orient=tk.HORIZONTAL, command=self.search_tree.xview)
        self.search_tree.configure(yscrollcommand=scroll_y.set, xscrollcommand=scroll_x.set)
        self.search_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scroll_y.pack(side=tk.RIGHT, fill=tk.Y)
        scroll_x.pack(side=tk.BOTTOM, fill=tk.X)

        self.search_tree.bind("<<TreeviewSelect>>", self.on_search_select)
        self.search_tree.bind("<Double-1>", self.on_search_double)

        detail_frame = ttk.LabelFrame(result_paned, text="选中行详情 (鼠标选择文字后按 Ctrl+C 复制)", padding="5")
        result_paned.add(detail_frame, weight=1)
        self.detail_text = scrolledtext.ScrolledText(detail_frame, height=6, bg="#fafafa",
                                                     font=("Consolas", 9), wrap=tk.WORD)
        self.detail_text.pack(fill=tk.BOTH, expand=True)
        self.detail_text.config(state=tk.DISABLED)

        self.search_info = ttk.Label(table_frame, text="", foreground="#7f8c8d")
        self.search_info.pack(side=tk.BOTTOM, pady=2)

    # ================== 数据管理标签页 ==================
    def create_management_tab(self):
        self.manage_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.manage_tab, text="数据管理")

        # 导入区域
        import_frame = ttk.LabelFrame(self.manage_tab, text="导入Excel/CSV（支持选择文件夹，含子文件夹）", padding="10")
        import_frame.pack(fill=tk.X, padx=5, pady=5)

        path_frame = ttk.Frame(import_frame)
        path_frame.pack(fill=tk.X, pady=5)
        ttk.Label(path_frame, text="目标文件夹:").pack(side=tk.LEFT, padx=5)
        self.folder_path_var = tk.StringVar(value=str(DEFAULT_DIR))
        path_entry = ttk.Entry(path_frame, textvariable=self.folder_path_var, width=50)
        path_entry.pack(side=tk.LEFT, padx=5)
        ttk.Button(path_frame, text="浏览", command=self.browse_folder).pack(side=tk.LEFT, padx=5)

        btn_frame = ttk.Frame(import_frame)
        btn_frame.pack(fill=tk.X, pady=5)
        self.import_btn = ttk.Button(btn_frame, text="开始导入", command=self.do_import)
        self.import_btn.pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="重置为默认 (Analyse_file)", command=self.reset_default_path).pack(side=tk.LEFT, padx=5)
        self.import_status = ttk.Label(btn_frame, text="", foreground="#7f8c8d")
        self.import_status.pack(side=tk.LEFT, padx=20)

        # 手工添加
        add_frame = ttk.LabelFrame(self.manage_tab, text="手工添加记录", padding="10")
        add_frame.pack(fill=tk.X, padx=5, pady=5)

        fields = [
            ("Sequence ID *", "add_id"),
            ("CDR3_aa *", "add_cdr3"),
            ("CDR1_aa", "add_cdr1"),
            ("CDR2_aa", "add_cdr2"),
            ("V_call", "add_v"),
            ("D_call", "add_d"),
            ("J_call", "add_j"),
            ("Junction_aa", "add_junc"),
            ("Junction_aa_length", "add_junc_len"),
            ("V_identity", "add_v_id"),
            ("D_identity", "add_d_id"),
            ("J_identity", "add_j_id"),
            ("Productive", "add_prod"),
            ("Sequence_aa", "add_seq_aa"),
        ]
        self.add_entries = {}
        for i, (label, key) in enumerate(fields):
            ttk.Label(add_frame, text=label).grid(row=i//2, column=(i%2)*2, sticky=tk.W, padx=5, pady=2)
            entry = ttk.Entry(add_frame, width=30)
            entry.grid(row=i//2, column=(i%2)*2+1, padx=5, pady=2)
            self.add_entries[key] = entry

        add_btn_frame = ttk.Frame(add_frame)
        add_btn_frame.grid(row=len(fields)//2 + 1, column=0, columnspan=4, pady=10)
        ttk.Button(add_btn_frame, text="添加/更新记录", command=self.do_add_record).pack(side=tk.LEFT, padx=5)
        ttk.Button(add_btn_frame, text="清空字段", command=self.clear_add_fields).pack(side=tk.LEFT, padx=5)

        # 批量分析FASTA
        fasta_frame = ttk.LabelFrame(self.manage_tab, text="批量分析FASTA文件（检查CDR3是否存在于数据库）", padding="10")
        fasta_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        fasta_btn_frame = ttk.Frame(fasta_frame)
        fasta_btn_frame.pack(fill=tk.X, pady=5)
        ttk.Button(fasta_btn_frame, text="选择FASTA文件并分析", command=self.do_analyze_fasta).pack(side=tk.LEFT, padx=5)
        self.fasta_status = ttk.Label(fasta_btn_frame, text="", foreground="#7f8c8d")
        self.fasta_status.pack(side=tk.LEFT, padx=20)
        ttk.Button(fasta_btn_frame, text="导出分析结果", command=self.export_fasta_results).pack(side=tk.LEFT, padx=5)

        result_frame = ttk.Frame(fasta_frame)
        result_frame.pack(fill=tk.BOTH, expand=True, pady=5)
        cols = ("FASTA_ID", "CDR3序列", "是否存在", "匹配的数据库ID")
        self.fasta_tree = ttk.Treeview(result_frame, columns=cols, show="headings", height=8)
        for col in cols:
            self.fasta_tree.heading(col, text=col)
            self.fasta_tree.column(col, width=150)
        scroll_y = ttk.Scrollbar(result_frame, orient=tk.VERTICAL, command=self.fasta_tree.yview)
        scroll_x = ttk.Scrollbar(result_frame, orient=tk.HORIZONTAL, command=self.fasta_tree.xview)
        self.fasta_tree.configure(yscrollcommand=scroll_y.set, xscrollcommand=scroll_x.set)
        self.fasta_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scroll_y.pack(side=tk.RIGHT, fill=tk.Y)
        scroll_x.pack(side=tk.BOTTOM, fill=tk.X)

        self.fasta_results = []

    # ================== 通用方法 ==================
    def update_stats(self):
        try:
            count = get_db_stats()
            self.stats_label.config(text=f"数据库: {count} 条记录")
        except:
            self.stats_label.config(text="数据库: 未连接")

    def clear_db(self):
        if messagebox.askyesno("确认清空", "确定要清空数据库中的所有数据吗？此操作不可恢复！"):
            try:
                conn = sqlite3.connect(DB_PATH)
                conn.execute("DELETE FROM cdr3_data")
                conn.commit()
                conn.close()
                self.update_stats()
                self.clear_results()
                messagebox.showinfo("成功", "数据库已清空")
            except Exception as e:
                messagebox.showerror("错误", f"清空失败: {e}")

    def export_all(self):
        file_path = filedialog.asksaveasfilename(defaultextension=".xlsx", filetypes=[("Excel files", "*.xlsx")])
        if not file_path:
            return
        try:
            conn = sqlite3.connect(DB_PATH)
            df = pd.read_sql("SELECT * FROM cdr3_data", conn)
            conn.close()
            if df.empty:
                messagebox.showwarning("警告", "数据库为空")
                return
            if 'full_row_json' in df.columns:
                df.drop('full_row_json', axis=1, inplace=True)
            df.to_excel(file_path, index=False)
            messagebox.showinfo("成功", f"导出成功，共 {len(df)} 条记录")
        except Exception as e:
            messagebox.showerror("错误", f"导出失败: {e}")

    # ================== 导入相关 ==================
    def browse_folder(self):
        folder = filedialog.askdirectory()
        if folder:
            self.folder_path_var.set(folder)

    def reset_default_path(self):
        self.folder_path_var.set(str(DEFAULT_DIR))

    def do_import(self):
        folder_path = self.folder_path_var.get().strip()
        if not folder_path:
            messagebox.showwarning("警告", "请先选择文件夹")
            return
        self.import_btn.config(state=tk.DISABLED, text="导入中...")
        self.root.update()
        try:
            new_count, err = import_excel_files(folder_path)
            self.update_stats()
            if err:
                messagebox.showwarning("导入完成", f"新增 {new_count} 条记录\n\n警告:\n{err}")
            else:
                messagebox.showinfo("导入完成", f"新增 {new_count} 条记录")
            self.import_status.config(text=f"上次导入: {new_count} 条新增")
        except Exception as e:
            messagebox.showerror("导入失败", str(e))
            self.import_status.config(text="导入失败")
        finally:
            self.import_btn.config(state=tk.NORMAL, text="开始导入")

    # ================== 搜索相关 ==================
    def do_search(self):
        seq = self.search_entry.get().strip().upper()
        if not seq:
            messagebox.showwarning("警告", "请输入序列")
            return
        region = self.region_var.get()
        exact = self.exact_var.get()
        try:
            results = search_region(region, seq, exact)
        except Exception as e:
            messagebox.showerror("搜索失败", str(e))
            return
        for item in self.search_tree.get_children():
            self.search_tree.delete(item)
        self.detail_text.config(state=tk.NORMAL)
        self.detail_text.delete(1.0, tk.END)
        self.detail_text.config(state=tk.DISABLED)
        if not results:
            self.search_info.config(text="未找到匹配记录")
            return
        for row in results:
            seq_aa = safe_str(row.get('sequence_aa'))
            display_seq_aa = seq_aa[:50] + '...' if len(seq_aa) > 50 else seq_aa
            values = (
                safe_str(row.get('sequence_id')),
                safe_str(row.get('cdr1_aa')),
                safe_str(row.get('cdr2_aa')),
                safe_str(row.get('cdr3_aa')),
                safe_str(row.get('v_call')),
                safe_str(row.get('j_call')),
                display_seq_aa
            )
            self.search_tree.insert("", tk.END, values=values, tags=(json.dumps(row),))
        self.search_info.config(text=f"找到 {len(results)} 条记录")

    def clear_results(self):
        for item in self.search_tree.get_children():
            self.search_tree.delete(item)
        self.search_info.config(text="")
        self.search_entry.delete(0, tk.END)
        self.detail_text.config(state=tk.NORMAL)
        self.detail_text.delete(1.0, tk.END)
        self.detail_text.config(state=tk.DISABLED)

    def on_search_select(self, event):
        """单击行时下方详情显示，包含 full_row_json 解析"""
        selection = self.search_tree.selection()
        if not selection:
            return
        item = selection[0]
        tags = self.search_tree.item(item, "tags")
        if not tags:
            return
        try:
            row_dict = json.loads(tags[0])
        except:
            return
        lines = []
        for key, value in row_dict.items():
            if key.startswith('created_time'):
                continue
            if key == 'full_row_json' and value:
                try:
                    extra_data = json.loads(value)
                    lines.append("--- 原始完整行 (full_row_json) ---")
                    for k, v in extra_data.items():
                        lines.append(f"  {k}: {safe_str(v)}")
                    continue
                except:
                    pass
            lines.append(f"{key}: {safe_str(value)}")
        detail = "\n".join(lines)
        self.detail_text.config(state=tk.NORMAL)
        self.detail_text.delete(1.0, tk.END)
        self.detail_text.insert(1.0, detail)
        self.detail_text.config(state=tk.DISABLED)

    def on_search_double(self, event):
        """双击弹出可滚动、可复制的详情窗口，包含 full_row_json"""
        selection = self.search_tree.selection()
        if not selection:
            return
        item = selection[0]
        tags = self.search_tree.item(item, "tags")
        if tags:
            try:
                row_dict = json.loads(tags[0])
                self.show_detail_dialog(row_dict)
            except:
                pass

    def show_detail_dialog(self, row_dict):
        """自定义详情窗口，显示所有字段（包括 full_row_json 解析）"""
        lines = []
        for k, v in row_dict.items():
            if k.startswith('created_time'):
                continue
            if k == 'full_row_json' and v:
                try:
                    extra_data = json.loads(v)
                    lines.append("--- 原始完整行 (full_row_json) ---")
                    for sub_k, sub_v in extra_data.items():
                        lines.append(f"  {sub_k}: {safe_str(sub_v)}")
                    continue
                except:
                    pass
            lines.append(f"{k}: {safe_str(v)}")
        detail_text = "\n".join(lines)

        win = tk.Toplevel(self.root)
        win.title("详细信息")
        win.geometry("800x600")
        win.minsize(600, 400)

        frame = ttk.Frame(win)
        frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        text_widget = scrolledtext.ScrolledText(frame, wrap=tk.NONE, font=("Consolas", 10))
        text_widget.pack(fill=tk.BOTH, expand=True)

        text_widget.insert(tk.END, detail_text)
        text_widget.config(state=tk.DISABLED)

        btn_frame = ttk.Frame(win)
        btn_frame.pack(fill=tk.X, pady=5)
        ttk.Button(btn_frame, text="关闭", command=win.destroy).pack(side=tk.RIGHT, padx=5)
        ttk.Button(btn_frame, text="复制全部",
                   command=lambda: (win.clipboard_clear(),
                                    win.clipboard_append(detail_text),
                                    messagebox.showinfo("提示", "已复制全部内容"))).pack(side=tk.RIGHT, padx=5)

    def export_search_results(self):
        items = self.search_tree.get_children()
        if not items:
            messagebox.showwarning("警告", "没有搜索结果可导出")
            return
        data = []
        for item in items:
            tags = self.search_tree.item(item, "tags")
            if tags:
                try:
                    data.append(json.loads(tags[0]))
                except:
                    continue
        if not data:
            return
        df = pd.DataFrame(data)
        file_path = filedialog.asksaveasfilename(defaultextension=".xlsx", filetypes=[("Excel files", "*.xlsx")])
        if file_path:
            df.to_excel(file_path, index=False)
            messagebox.showinfo("成功", f"已导出 {len(df)} 条记录")

    # ================== 数据管理 ==================
    def clear_add_fields(self):
        for entry in self.add_entries.values():
            entry.delete(0, tk.END)

    def do_add_record(self):
        seq_id = self.add_entries['add_id'].get().strip()
        cdr3 = self.add_entries['add_cdr3'].get().strip().upper()
        if not seq_id or not cdr3:
            messagebox.showwarning("警告", "Sequence ID 和 CDR3_aa 为必填字段")
            return

        conn = sqlite3.connect(DB_PATH)
        cursor = conn.cursor()
        cursor.execute("SELECT 1 FROM cdr3_data WHERE sequence_id = ?", (seq_id,))
        exists = cursor.fetchone() is not None
        conn.close()

        if exists:
            if not messagebox.askyesno("确认覆盖", f"ID '{seq_id}' 已存在，是否覆盖该记录？"):
                return

        fields = {
            'cdr1_aa': self.add_entries['add_cdr1'].get().strip(),
            'cdr2_aa': self.add_entries['add_cdr2'].get().strip(),
            'v_call': self.add_entries['add_v'].get().strip(),
            'd_call': self.add_entries['add_d'].get().strip(),
            'j_call': self.add_entries['add_j'].get().strip(),
            'junction_aa': self.add_entries['add_junc'].get().strip(),
            'junction_aa_length': self.add_entries['add_junc_len'].get().strip(),
            'v_identity': self.add_entries['add_v_id'].get().strip(),
            'd_identity': self.add_entries['add_d_id'].get().strip(),
            'j_identity': self.add_entries['add_j_id'].get().strip(),
            'productive': self.add_entries['add_prod'].get().strip(),
            'sequence_aa': self.add_entries['add_seq_aa'].get().strip()
        }
        try:
            if fields['junction_aa_length']:
                fields['junction_aa_length'] = int(fields['junction_aa_length'])
            else:
                fields['junction_aa_length'] = None
        except ValueError:
            messagebox.showwarning("警告", "Junction_aa_length 必须是整数")
            return
        for key in ['v_identity', 'd_identity', 'j_identity']:
            if fields[key]:
                try:
                    fields[key] = float(fields[key])
                except ValueError:
                    messagebox.showwarning("警告", f"{key} 必须是数字")
                    return
            else:
                fields[key] = None

        success, msg = add_manual_record(
            seq_id=seq_id,
            cdr3_aa=cdr3,
            cdr1_aa=fields['cdr1_aa'],
            cdr2_aa=fields['cdr2_aa'],
            v_call=fields['v_call'],
            d_call=fields['d_call'],
            j_call=fields['j_call'],
            junction_aa=fields['junction_aa'],
            junction_aa_length=fields['junction_aa_length'],
            v_identity=fields['v_identity'],
            d_identity=fields['d_identity'],
            j_identity=fields['j_identity'],
            productive=fields['productive'],
            sequence_aa=fields['sequence_aa']
        )
        if success:
            messagebox.showinfo("成功", msg)
            self.update_stats()
            self.clear_add_fields()
        else:
            messagebox.showerror("错误", f"添加失败: {msg}")

    def do_analyze_fasta(self):
        file_path = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fasta *.fa *.faa"), ("All files", "*.*")])
        if not file_path:
            return
        self.fasta_status.config(text="分析中...")
        self.root.update()
        try:
            results = analyze_fasta_file(file_path)
            self.fasta_results = results
            for item in self.fasta_tree.get_children():
                self.fasta_tree.delete(item)
            exists_count = sum(1 for r in results if r['exists'])
            for r in results:
                seq_display = r['cdr3_seq'][:50] + ('...' if len(r['cdr3_seq']) > 50 else '')
                self.fasta_tree.insert("", tk.END, values=(
                    r['seq_id'],
                    seq_display,
                    "是" if r['exists'] else "否",
                    r['matched_ids']
                ))
            self.fasta_status.config(text=f"分析完成: 共 {len(results)} 条，存在 {exists_count} 条，不存在 {len(results)-exists_count} 条")
        except Exception as e:
            messagebox.showerror("分析失败", str(e))
            self.fasta_status.config(text="分析失败")

    def export_fasta_results(self):
        if not self.fasta_results:
            messagebox.showwarning("警告", "请先执行FASTA分析")
            return
        df = pd.DataFrame(self.fasta_results)
        file_path = filedialog.asksaveasfilename(defaultextension=".xlsx", filetypes=[("Excel files", "*.xlsx")])
        if file_path:
            df.to_excel(file_path, index=False)
            messagebox.showinfo("成功", f"已导出 {len(df)} 条分析结果")

def main():
    root = tk.Tk()
    app = CDR3SearchApp(root)
    root.mainloop()

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print("程序异常:", str(e))
        traceback.print_exc()
        input("按回车退出...")