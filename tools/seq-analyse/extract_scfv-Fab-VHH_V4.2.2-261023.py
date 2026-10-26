# -*- coding: utf-8 -*-
"""
抗体测序数据处理系统 V4.2.2
支持 scFv / Fab / VHH 三种分析模式，提供序列提取、翻译、去重、
多序列比对、IgBlast 注释、重复序列整理与结果预览。
"""
# 修订记录：
# V1.1：优化 if 查找，使之定义更精确
# V1.2：精简 if 条件代码
# V1.3：增加翻译功能，生成 7 个文件，单独列出有问题的序列
# V2.0：防止测序文件是多行 FASTA 格式
# V3.0：调用本地 IgBlast 处理生成的核酸序列
# V3.X：scfv_aa 去重并计算个数，生成比对序列以及进化树（可用 MEGA 打开）
# V3.Y：增加 rabbit、rat 数据库；根据氨基酸序列去重后，对 scfv_aa 去重后的核酸序列进行 IgBlast
# V3.Y.1（230307）：优化序列切片，尽量避免 GS 区过多引起程序异常终止
# V3.AI（20260120）：去重引用 seqkit 提高速度，整体代码使用 DS 优化
# V4.0（20260615）：重写代码提高性能；IgBlast 升级为 1.22 以匹配 J 基因模糊匹配；适配 Python 3.14
# V4.1（20261026）：增加 Fab 分析模式；集成 VHH 分析模块，支持 scFv/Fab/VHH 三种模式切换
# V4.1.1（20261026）：删除冗余备份代码，清理未使用导入与变量；重构 GUI 布局与视觉样式（模式单选、文件树、结果双击预览、统一日志样式），提升专业性与易用性
# V4.1.2（20261026）：新增分隔序列自定义功能；默认沿用各模式标准序列，支持按模式手工修改、恢复默认并校验输入
# V4.2.0（20261026）：稳健性重构：修复重复明细/重复ID/部分链翻译/去重返回值/输出顺序/位点匹配/IgBlast解析；统一绝对路径与工具参数列表；增加任务快照、步骤状态、旧结果归档、原子写入、日志限流、预检、运行清单、输入格式扩展与翻译质控
# V4.2.1（20261026）：预检按BLAST索引前缀识别数据库；参数设置新增输入格式多选，默认 .seq，支持 .txt/.fa/.fas/.fasta/.fna 及对应 gzip
# V4.2.2（20261026）：修复seqkit版本检测命令；工具解析优先bin/程序目录并统一返回绝对路径

import sys
import os
import subprocess
import threading
import queue
import shlex
import shutil
import json
import gzip
import fnmatch
from pathlib import Path
from datetime import datetime
import traceback

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

# tkinter imports
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
from tkinter import font as tkfont


def get_base_dir():
    """获取可执行文件或脚本所在目录（兼容打包后）"""
    if getattr(sys, 'frozen', False):
        # 打包后的 exe
        return os.path.dirname(sys.executable)
    else:
        # 普通脚本
        return os.path.dirname(os.path.abspath(__file__))


def resolve_clustalw(base_dir=None):
    """按程序目录、bin 目录和 PATH 查找 clustalw"""
    candidates = []
    if base_dir:
        candidates.append(Path(base_dir) / "clustalw2.exe")
    base_path = Path(get_base_dir())
    candidates.extend([
        base_path / "clustalw2.exe",
        base_path / "bin" / "clustalw2.exe",
        base_path / "clustalw2",
        base_path / "bin" / "clustalw2",
    ])
    try:
        candidates.extend([
            Config.BIN_DIR / "clustalw2.exe",
            Config.BIN_DIR / "clustalw2",
        ])
    except NameError:
        pass
    for candidate in candidates:
        if candidate.is_file():
            return str(candidate)
    return shutil.which("clustalw2") or shutil.which("clustalw")


# 常量定义
class Config:
    """配置常量"""
    BASE_DIR = Path(
        os.environ.get("ANTIBODY_ANALYSIS_BASE", get_base_dir())
    ).expanduser().resolve()
    ANALYSE_DIR = BASE_DIR / "Analyse_file"
    DATABASE_DIR = BASE_DIR / "database"
    BIN_DIR = BASE_DIR / "bin"
    OPTIONAL_FILE_DIR = BASE_DIR / "optional_file"

    # scFv模式：三个linker按顺序切分VH-linker-VL
    SCFV_LINKERS = {
        'linker1': 'GGCCCAGGCGGCC',
        'linker2': 'GGTGGCGGCGGTTCAGGCGGAGGTGGCTCTGGCGGTGGCGGATCG',
        'linker3': 'GGCCAGGCCGGCC'
    }

    # Fab模式：VH/VL分别位于不同的.seq文件中，通过特征上下游序列定位
    FAB_FLANK_PATTERNS = {
        'VH': {
            'upstream': 'CTCGCTGCCCAACCAGCCATGGCC',
            'downstream': 'GCTTCCACCAAGGGCCCGTCCGTG'
        },
        'VL': {
            'upstream': 'TTCGCTACCGTGGCCCAGGCGGCC',
            'downstream': 'AGAACAGTGGCTGCTCCGTCCGTG'
        }
    }

    # VHH模式：通过SfiI位点定位VHH序列
    VHH_SITES = {
        'upstream': 'GGCCCAGGCGGCC',
        'downstream': 'GGCCAGGCCGGCC'
    }

    # 输入文件格式选项
    INPUT_FORMAT_OPTIONS = ('.seq', '.txt', '.fa', '.fas', '.fasta', '.fna')
    DEFAULT_INPUT_FORMATS = ('.seq',)

    # 物种数据库映射
    SPECIES_DB = {
        'human': {
            'vh': 'Human_IGHV.blastdb',
            'd': 'Human_IGHD.blastdb',
            'jh': 'Human_IGHJ.blastdb',
            'vl': 'Human_IGLV.blastdb',
            'jl': 'Human_IGLJ.blastdb',
            'aux': 'human_gl.aux'
        },
        'mouse': {
            'vh': 'Mouse_IGHV.blastdb',
            'd': 'Mouse_IGHD.blastdb',
            'jh': 'Mouse_IGHJ.blastdb',
            'vl': 'Mouse_IGLV.blastdb',
            'jl': 'Mouse_IGLJ.blastdb',
            'aux': 'mouse_gl.aux'
        },
        'rabbit': {
            'vh': 'Rabbit_IGHV.blastdb',
            'd': 'Rabbit_IGHD.blastdb',
            'jh': 'Rabbit_IGHJ.blastdb',
            'vl': 'Rabbit_IGLV.blastdb',
            'jl': 'Rabbit_IGLJ.blastdb',
            'aux': 'rabbit_gl.aux'
        },
        'rat': {
            'vh': 'Rat_IGHV.blastdb',
            'd': 'Rat_IGHD.blastdb',
            'jh': 'Rat_IGHJ.blastdb',
            'vl': 'Rat_IGLV.blastdb',
            'jl': 'Rat_IGLJ.blastdb',
            'aux': 'rat_gl.aux'
        },
        'alpaca': {
            'vh': 'Alpaca_IGHV.blastdb',
            'd': 'Alpaca_IGHD.blastdb',
            'jh': 'Alpaca_IGHJ.blastdb',
            'aux': 'Alpaca_gl.aux'
        },
        'human_VLL': {
            'vh': 'Human_IGLV.blastdb',
            'd': 'Human_IGHD.blastdb',
            'jh': 'Human_IGLJ.blastdb',
            'vl': 'Human_IGLV.blastdb',
            'jl': 'Human_IGLJ.blastdb',
            'aux': 'human_gl.aux'
        }
    }

    # 各分析模式可选的物种
    MODE_SPECIES = {
        'scfv': ['human', 'mouse', 'rabbit', 'rat'],
        'fab': ['human', 'mouse', 'rabbit', 'rat'],
        'vhh': ['human', 'alpaca', 'human_VLL']
    }


def run_clustalw(infile_path, base_dir=None):
    """
    使用 subprocess 调用 clustalw2.exe 进行多序列比对
    参数:
        infile_path: 输入文件的完整路径
        base_dir: 程序所在目录（clustalw2.exe 所在目录），如果为 None 则自动获取
    返回:
        (success, message)  success为布尔值，message为信息字符串
    """
    if base_dir is None:
        base_dir = get_base_dir()

    # 按程序目录、bin 目录和 PATH 查找 clustalw
    clustalw_exe = resolve_clustalw(base_dir)
    if not clustalw_exe:
        return False, "错误: 未找到 clustalw2，请将其放在程序目录或 bin 目录，或加入 PATH"

    # 检查输入文件是否存在
    if not os.path.exists(infile_path):
        return False, f"错误: 输入文件不存在: {infile_path}"

    # 检查序列条数，避免单序列卡死
    try:
        seq_count = sum(1 for _ in SeqIO.parse(infile_path, 'fasta'))
        if seq_count < 2:
            return False, f"跳过: 序列数量少于2条 ({seq_count}条) - {os.path.basename(infile_path)}"
    except Exception as e:
        return False, f"警告: 无法解析序列文件: {str(e)}"

    # 切换到输入文件所在目录执行（clustalw2 会在当前目录生成输出文件）
    work_dir = os.path.dirname(infile_path)
    filename = os.path.basename(infile_path)
    aln_path = os.path.splitext(infile_path)[0] + '.aln'
    dnd_path = os.path.splitext(infile_path)[0] + '.dnd'
    for stale_path in (aln_path, dnd_path):
        if os.path.exists(stale_path):
            try:
                os.remove(stale_path)
            except OSError as e:
                return False, f"无法清理旧比对文件 {stale_path}: {str(e)}"

    # 构建命令
    cmd = [clustalw_exe, '-infile=' + filename]

    try:
        run_kwargs = {}
        if os.name == 'nt':
            startupinfo_class = getattr(subprocess, "STARTUPINFO", None)
            if startupinfo_class is not None:
                startupinfo = startupinfo_class()
                startupinfo.dwFlags |= getattr(subprocess, "STARTF_USESHOWWINDOW", 0)
                startupinfo.wShowWindow = getattr(subprocess, "SW_HIDE", 0)
                run_kwargs['startupinfo'] = startupinfo
            run_kwargs['creationflags'] = getattr(subprocess, "CREATE_NO_WINDOW", 0)
        result = subprocess.run(
            cmd,
            cwd=work_dir,
            capture_output=True,
            text=True,
            timeout=300,  # 5分钟超时
            **run_kwargs
        )

        if result.returncode == 0:
            if not os.path.exists(aln_path) or os.path.getsize(aln_path) == 0:
                return False, f"比对返回成功但未生成有效 .aln 文件: {filename}"
            return True, f"比对成功: {filename}"
        else:
            error_msg = result.stderr if result.stderr else result.stdout
            return False, f"比对失败: {filename} - {error_msg[:200]}"

    except subprocess.TimeoutExpired:
        return False, f"比对超时 (5分钟): {filename}"
    except Exception as e:
        return False, f"执行出错: {filename} - {str(e)}"


class LogQueueHandler:
    """日志队列处理器"""

    def __init__(self, log_callback):
        self.log_callback = log_callback
        self.queue = queue.Queue()

    def put(self, message, level="INFO"):
        timestamp = datetime.now().strftime("%H:%M:%S")
        self.queue.put((timestamp, level, message))

    def process_queue(self, max_items=200):
        processed = 0
        while max_items is None or processed < max_items:
            try:
                timestamp, level, message = self.queue.get_nowait()
            except queue.Empty:
                break
            self.log_callback(timestamp, level, message)
            processed += 1
        return processed

    def clear(self):
        while True:
            try:
                self.queue.get_nowait()
            except queue.Empty:
                break


class AntibodyGUI:
    """主GUI类"""

    APP_VERSION = "V4.2.2"
    MAX_LOG_LINES = 5000
    MAX_OUTPUT_LINES = 3000
    MODE_INFO = {
        "scfv": ("scFv", "单文件 VH-linker-VL"),
        "fab": ("Fab", "VH / VL 分文件"),
        "vhh": ("VHH", "上下游位点间 VHH"),
    }
    MODE_DESCRIPTIONS = {
        "scfv": "单个序列文件中包含 VH-linker-VL，适用于 scFv 抗体序列。",
        "fab": "VH 和 VL 分别位于不同的 .seq 文件中，适用于 Fab 抗体序列。",
        "vhh": "单个序列文件中包含上下游位点之间的 VHH，默认使用 SfiI 位点，可在参数设置中自定义。",
    }
    SEPARATOR_SPECS = {
        "scfv": [
            ("linker1", "Linker 1", Config.SCFV_LINKERS["linker1"]),
            ("linker2", "Linker 2", Config.SCFV_LINKERS["linker2"]),
            ("linker3", "Linker 3", Config.SCFV_LINKERS["linker3"]),
        ],
        "fab": [
            ("vh_upstream", "VH 上游", Config.FAB_FLANK_PATTERNS["VH"]["upstream"]),
            ("vh_downstream", "VH 下游", Config.FAB_FLANK_PATTERNS["VH"]["downstream"]),
            ("vl_upstream", "VL 上游", Config.FAB_FLANK_PATTERNS["VL"]["upstream"]),
            ("vl_downstream", "VL 下游", Config.FAB_FLANK_PATTERNS["VL"]["downstream"]),
        ],
        "vhh": [
            ("upstream", "上游位点", Config.VHH_SITES["upstream"]),
            ("downstream", "下游位点", Config.VHH_SITES["downstream"]),
        ],
    }

    def __init__(self, root):
        self.root = root
        self.root.title(f"抗体测序数据处理系统 {self.APP_VERSION}")
        self.root.minsize(1000, 660)

        # 初始化变量
        self.seq_files = []
        self.analysis_mode = tk.StringVar(value="scfv")
        self.current_species = tk.StringVar(value="human")
        self.current_scheme = tk.StringVar(value="imgt")
        self.is_processing = False
        self.include_subfolders = tk.BooleanVar(value=False)
        self.skip_alignment = tk.BooleanVar(value=False)
        self.input_format_vars = {
            ext: tk.BooleanVar(value=ext in Config.DEFAULT_INPUT_FORMATS)
            for ext in Config.INPUT_FORMAT_OPTIONS
        }
        self.separator_vars = {
            mode: {
                key: tk.StringVar(value=default)
                for key, _label, default in specs
            }
            for mode, specs in self.SEPARATOR_SPECS.items()
        }
        self._active_mode = None
        self._active_separators = None
        self.workflow_state = {
            "extracted": False,
            "translated": False,
            "deduplicated": False,
            "aligned": False,
        }

        # 缓存字典，避免重复读取文件
        self.vh_nt_dict = {}
        self.vl_nt_dict = {}
        self.scfv_nt_dict = {}
        self.vhh_nt_dict = {}

        self.log_handler = LogQueueHandler(self.add_log_message)
        self.ui_queue = queue.Queue()

        self.setup_styles()
        self.setup_ui()
        self._center_window(1160, 740)
        self.setup_directories()
        self.refresh_results()
        self.process_log_queue()
        self.root.after(50, self._process_ui_queue)

    def _center_window(self, width, height):
        """将窗口居中显示，并避免超出屏幕可用区域"""
        self.root.update_idletasks()
        screen_w = self.root.winfo_screenwidth()
        screen_h = self.root.winfo_screenheight()
        width = min(width, max(screen_w - 40, 1000))
        height = min(height, max(screen_h - 80, 660))
        x = max((screen_w - width) // 2, 0)
        y = max((screen_h - height) // 2, 0)
        self.root.geometry(f"{width}x{height}+{x}+{y}")

    def setup_styles(self):
        style = ttk.Style()
        try:
            style.theme_use('clam')
        except tk.TclError:
            pass

        available_fonts = set(tkfont.families())
        available_lower = {name.lower(): name for name in available_fonts}
        for family in (
            "Microsoft YaHei UI", "Microsoft YaHei", "PingFang SC",
            "Noto Sans CJK SC", "Source Han Sans SC", "WenQuanYi Zen Hei",
            "Song Ti", "FangSong Ti", "Segoe UI", "Helvetica", "Arial"
        ):
            if family.lower() in available_lower:
                self.font_family = available_lower[family.lower()]
                break
        else:
            self.font_family = "TkDefaultFont"

        for family in ("Consolas", "Cascadia Mono", "Courier New", "DejaVu Sans Mono", "Courier"):
            if family.lower() in available_lower:
                self.mono_font = available_lower[family.lower()]
                break
        else:
            self.mono_font = "TkFixedFont"

        self.bg_color = "#eef2f7"
        self.panel_color = "#ffffff"
        self.primary_color = "#2563eb"
        self.primary_dark = "#1e40af"
        self.text_color = "#1f2937"
        self.muted_color = "#6b7280"
        self.border_color = "#d7dee8"
        self.error_color = "#dc2626"
        self.warning_color = "#b45309"
        self.log_bg = "#fbfdff"

        self.root.configure(bg=self.bg_color)
        self.root.option_add("*Font", (self.font_family, 10))

        style.configure('.', font=(self.font_family, 10), background=self.bg_color, foreground=self.text_color)
        style.configure('TFrame', background=self.bg_color)
        style.configure('TLabel', background=self.bg_color, foreground=self.text_color)
        style.configure('Muted.TLabel', background=self.bg_color, foreground=self.muted_color)
        style.configure('Section.TLabel', font=(self.font_family, 10, 'bold'), background=self.bg_color, foreground=self.text_color)
        style.configure('TButton', padding=(12, 6))
        style.configure('Accent.TButton', background=self.primary_color, foreground='#ffffff', borderwidth=0, padding=(16, 8))
        style.map('Accent.TButton',
                  background=[('active', self.primary_dark), ('disabled', '#9ca3af')],
                  foreground=[('disabled', '#f3f4f6')])
        style.configure('TLabelframe', background=self.bg_color, bordercolor=self.border_color, relief='solid', borderwidth=1)
        style.configure('TLabelframe.Label', background=self.bg_color, foreground=self.text_color, font=(self.font_family, 10, 'bold'))
        style.configure('TNotebook', background=self.bg_color, borderwidth=0)
        style.configure('TNotebook.Tab', padding=(20, 9), font=(self.font_family, 10))
        style.map('TNotebook.Tab', background=[('selected', self.panel_color)], foreground=[('selected', self.primary_color)])
        style.configure('Treeview', background=self.panel_color, fieldbackground=self.panel_color,
                        foreground=self.text_color, rowheight=26, font=(self.font_family, 9), borderwidth=0)
        style.configure('Treeview.Heading', font=(self.font_family, 9, 'bold'),
                        background="#f1f5f9", foreground=self.text_color, relief='flat', padding=(8, 6))
        style.map('Treeview.Heading', background=[('active', "#e2e8f0")])
        style.configure('TRadiobutton', background=self.bg_color, foreground=self.text_color, font=(self.font_family, 10))
        style.map('TRadiobutton', background=[('active', self.bg_color)])
        style.configure('TCheckbutton', background=self.bg_color, foreground=self.text_color)
        style.map('TCheckbutton', background=[('active', self.bg_color)])
        style.configure('Horizontal.TProgressbar', thickness=8, background=self.primary_color,
                        troughcolor="#dbe3ee", borderwidth=0)

    def setup_ui(self):
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(1, weight=1)

        header = tk.Frame(self.root, bg=self.primary_color, height=76)
        header.grid(row=0, column=0, sticky='ew')
        header.grid_propagate(False)
        header.columnconfigure(0, weight=1)

        title_box = tk.Frame(header, bg=self.primary_color)
        title_box.grid(row=0, column=0, sticky='w', padx=24, pady=(12, 10))
        tk.Label(title_box, text="抗体测序数据处理系统", bg=self.primary_color, fg="#ffffff",
                 font=(self.font_family, 19, 'bold'), anchor='w').pack(anchor='w')
        tk.Label(title_box, text="scFv / Fab / VHH 序列提取 · 去重 · 比对 · IgBlast 注释",
                 bg=self.primary_color, fg="#dbeafe", font=(self.font_family, 10), anchor='w').pack(anchor='w', pady=(2, 0))

        badge_box = tk.Frame(header, bg=self.primary_color)
        badge_box.grid(row=0, column=1, sticky='e', padx=24)
        self.header_mode_label = tk.Label(
            badge_box, text="当前模式：scFv", bg=self.primary_dark, fg="#ffffff",
            font=(self.font_family, 9, 'bold'), padx=12, pady=5
        )
        self.header_mode_label.pack(side=tk.LEFT, padx=(0, 8))
        tk.Label(badge_box, text=self.APP_VERSION, bg="#ffffff", fg=self.primary_dark,
                 font=(self.font_family, 9, 'bold'), padx=10, pady=5).pack(side=tk.LEFT)

        body = ttk.Frame(self.root, padding=(16, 12, 16, 6))
        body.grid(row=1, column=0, sticky='nsew')
        body.columnconfigure(0, weight=1)
        body.rowconfigure(0, weight=1)

        self.notebook = ttk.Notebook(body)
        self.notebook.grid(row=0, column=0, sticky='nsew')

        self.create_process_tab()
        self.create_settings_tab()
        self.create_results_tab()
        self.create_log_tab()

        footer = ttk.Frame(self.root, padding=(16, 0, 16, 10))
        footer.grid(row=2, column=0, sticky='ew')
        footer.columnconfigure(0, weight=1)

        self.progress_var = tk.DoubleVar()
        self.progress_bar = ttk.Progressbar(footer, variable=self.progress_var, maximum=100, mode='determinate')
        self.progress_bar.grid(row=0, column=0, sticky='ew', pady=(0, 6))

        status_row = ttk.Frame(footer)
        status_row.grid(row=1, column=0, sticky='ew')
        status_row.columnconfigure(0, weight=1)
        self.status_label = ttk.Label(status_row, text="就绪", style='Muted.TLabel')
        self.status_label.grid(row=0, column=0, sticky='w')
        self.file_count_label = ttk.Label(status_row, text="已选择 0 个文件", style='Muted.TLabel')
        self.file_count_label.grid(row=0, column=1, sticky='e')

        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)

    def create_process_tab(self):
        tab = ttk.Frame(self.notebook, padding=(12, 12, 12, 12))
        self.notebook.add(tab, text="数据处理")
        tab.columnconfigure(0, weight=0, minsize=430)
        tab.columnconfigure(1, weight=1)
        tab.rowconfigure(0, weight=1)

        left = ttk.Frame(tab)
        left.grid(row=0, column=0, sticky='nsew', padx=(0, 10))
        left.columnconfigure(0, weight=1)
        left.rowconfigure(1, weight=1)

        right = ttk.Frame(tab)
        right.grid(row=0, column=1, sticky='nsew', padx=(10, 0))
        right.columnconfigure(0, weight=1)
        right.rowconfigure(1, weight=1)

        mode_frame = ttk.LabelFrame(left, text="分析模式", padding=(12, 10))
        mode_frame.grid(row=0, column=0, sticky='ew', pady=(0, 10))
        mode_frame.columnconfigure(1, weight=1)

        self.mode_radios = []
        for row, (mode, (name, hint)) in enumerate(self.MODE_INFO.items()):
            radio = ttk.Radiobutton(
                mode_frame, text=name, value=mode, variable=self.analysis_mode,
                command=self._on_mode_changed
            )
            radio.grid(row=row, column=0, sticky='w', pady=3, padx=(0, 12))
            self.mode_radios.append(radio)
            ttk.Label(mode_frame, text=hint, style='Muted.TLabel').grid(row=row, column=1, sticky='w', pady=3)

        self.mode_label = ttk.Label(
            mode_frame, text=self.MODE_DESCRIPTIONS['scfv'], style='Muted.TLabel',
            wraplength=380, justify='left'
        )
        self.mode_label.grid(row=len(self.MODE_INFO), column=0, columnspan=2, sticky='ew', pady=(8, 4))

        file_frame = ttk.LabelFrame(left, text="序列文件", padding=(12, 10))
        file_frame.grid(row=1, column=0, sticky='nsew')
        file_frame.columnconfigure(0, weight=1)
        file_frame.rowconfigure(0, weight=1)

        tree_box = ttk.Frame(file_frame)
        tree_box.grid(row=0, column=0, sticky='nsew')
        tree_box.columnconfigure(0, weight=1)
        tree_box.rowconfigure(0, weight=1)

        self.file_tree = ttk.Treeview(
            tree_box, columns=("name", "folder"), show="headings",
            height=7, selectmode="extended"
        )
        self.file_tree.heading("name", text="文件名")
        self.file_tree.heading("folder", text="所在目录")
        self.file_tree.column("name", width=200, minwidth=140, stretch=True)
        self.file_tree.column("folder", width=170, minwidth=100, stretch=True)
        file_scroll = ttk.Scrollbar(tree_box, orient=tk.VERTICAL, command=self.file_tree.yview)
        self.file_tree.configure(yscrollcommand=file_scroll.set)
        self.file_tree.grid(row=0, column=0, sticky='nsew')
        file_scroll.grid(row=0, column=1, sticky='ns')
        self.file_tree.tag_configure('odd', background="#f8fafc")
        self.file_tree.tag_configure('even', background=self.panel_color)

        action_row = ttk.Frame(file_frame)
        action_row.grid(row=1, column=0, sticky='ew', pady=(10, 6))
        self.select_files_btn = ttk.Button(action_row, text="选择文件", command=self.select_seq_files)
        self.select_files_btn.pack(side=tk.LEFT, padx=(0, 6))
        self.select_folder_btn = ttk.Button(action_row, text="选择文件夹", command=self.select_seq_folder)
        self.select_folder_btn.pack(side=tk.LEFT, padx=(0, 6))
        self.remove_files_btn = ttk.Button(action_row, text="移除选中", command=self.remove_selected_files)
        self.remove_files_btn.pack(side=tk.LEFT, padx=(0, 6))
        self.clear_files_btn = ttk.Button(action_row, text="清空列表", command=self.clear_file_list)
        self.clear_files_btn.pack(side=tk.LEFT)

        option_row = ttk.Frame(file_frame)
        option_row.grid(row=2, column=0, sticky='ew')
        self.include_subfolders_check = ttk.Checkbutton(
            option_row, text="包含子文件夹", variable=self.include_subfolders
        )
        self.include_subfolders_check.pack(side=tk.LEFT)
        self.folder_info_label = ttk.Label(file_frame, text="", style='Muted.TLabel')
        self.folder_info_label.grid(row=3, column=0, sticky='w', pady=(6, 0))

        steps_frame = ttk.LabelFrame(right, text="处理步骤", padding=(12, 10))
        steps_frame.grid(row=0, column=0, sticky='ew', pady=(0, 10))
        steps_frame.columnconfigure(0, weight=1)

        steps_grid = ttk.Frame(steps_frame)
        steps_grid.grid(row=0, column=0, sticky='ew')
        for col in range(3):
            steps_grid.columnconfigure(col, weight=1)

        self.step1_btn = ttk.Button(steps_grid, text="提取序列", command=self.extract_sequences, state=tk.DISABLED)
        self.step1_btn.grid(row=0, column=0, sticky='ew', padx=4, pady=4)
        self.step2_btn = ttk.Button(steps_grid, text="翻译序列", command=self.translate_sequences, state=tk.DISABLED)
        self.step2_btn.grid(row=0, column=1, sticky='ew', padx=4, pady=4)
        self.step3_btn = ttk.Button(steps_grid, text="序列去重", command=self.remove_duplicates, state=tk.DISABLED)
        self.step3_btn.grid(row=0, column=2, sticky='ew', padx=4, pady=4)
        self.step4_btn = ttk.Button(steps_grid, text="多序列比对", command=self.run_alignment, state=tk.DISABLED)
        self.step4_btn.grid(row=1, column=0, sticky='ew', padx=4, pady=4)
        self.step5_btn = ttk.Button(steps_grid, text="IgBlast 分析", command=self.run_igblast, state=tk.DISABLED)
        self.step5_btn.grid(row=1, column=1, sticky='ew', padx=4, pady=4)
        self.step6_btn = ttk.Button(steps_grid, text="处理重复文件", command=self.process_duplicate_files, state=tk.DISABLED)
        self.step6_btn.grid(row=1, column=2, sticky='ew', padx=4, pady=4)

        run_row = ttk.Frame(steps_frame)
        run_row.grid(row=1, column=0, sticky='ew', pady=(8, 0))
        run_row.columnconfigure(0, weight=1)
        self.run_all_btn = ttk.Button(
            run_row, text="一键运行全部步骤", command=self.run_all_steps,
            style='Accent.TButton', state=tk.DISABLED
        )
        self.run_all_btn.grid(row=0, column=0, sticky='ew', padx=(0, 8))
        self.skip_alignment_check = ttk.Checkbutton(
            run_row, text="跳过比对", variable=self.skip_alignment
        )
        self.skip_alignment_check.grid(row=0, column=1, sticky='e')

        output_frame = ttk.LabelFrame(right, text="实时输出", padding=(12, 10))
        output_frame.grid(row=1, column=0, sticky='nsew')
        output_frame.columnconfigure(0, weight=1)
        output_frame.rowconfigure(0, weight=1)
        self.output_text = scrolledtext.ScrolledText(
            output_frame, width=40, height=10, bg=self.log_bg, fg=self.text_color,
            relief=tk.FLAT, borderwidth=1, font=(self.mono_font, 9), padx=8, pady=6
        )
        self.output_text.grid(row=0, column=0, sticky='nsew')

    def is_fab_mode(self, mode=None):
        """当前是否为Fab分析模式；可传入任务快照中的 mode"""
        return (mode or self.analysis_mode.get()) == "fab"

    def is_vhh_mode(self, mode=None):
        """当前是否为VHH分析模式；可传入任务快照中的 mode"""
        return (mode or self.analysis_mode.get()) == "vhh"

    def _update_species_options(self):
        """根据当前分析模式更新可选物种"""
        mode = self.analysis_mode.get()
        values = Config.MODE_SPECIES.get(mode, Config.MODE_SPECIES['scfv'])
        if hasattr(self, "species_combo"):
            self.species_combo.config(values=values)
        if self.current_species.get() not in values:
            self.current_species.set(values[0])
            self.log(f"物种已自动切换为: {values[0]}")

    def _get_separator_values(self, mode=None):
        """读取指定模式的分隔序列，未初始化控件时回退到默认值"""
        mode = mode or self.analysis_mode.get()
        values = {}
        var_map = getattr(self, "separator_vars", {}).get(mode, {})
        for key, _label, default in self.SEPARATOR_SPECS.get(mode, []):
            var = var_map.get(key)
            raw = var.get() if var is not None else default
            values[key] = ''.join(str(raw).split()).upper()
        return values

    def _render_separator_fields(self):
        """按当前分析模式重建分隔序列输入框"""
        if not hasattr(self, "separator_fields_frame"):
            return
        for child in self.separator_fields_frame.winfo_children():
            child.destroy()

        mode = self.analysis_mode.get()
        specs = self.SEPARATOR_SPECS.get(mode, [])
        per_row = 1 if len(specs) <= 2 else 2
        for col in range(4):
            self.separator_fields_frame.columnconfigure(col, weight=0)
        self.separator_fields_frame.columnconfigure(1, weight=1)
        if per_row == 2:
            self.separator_fields_frame.columnconfigure(3, weight=1)

        var_map = self.separator_vars.get(mode, {})
        for index, (key, label, _default) in enumerate(specs):
            row, slot = divmod(index, per_row)
            label_col = slot * 2
            entry_col = label_col + 1
            entry_span = 3 if per_row == 2 and slot == 0 and index == len(specs) - 1 else 1
            ttk.Label(self.separator_fields_frame, text=label).grid(
                row=row, column=label_col, sticky='w',
                padx=(0 if label_col == 0 else 18, 8), pady=4
            )
            entry = ttk.Entry(
                self.separator_fields_frame, textvariable=var_map[key], width=32
            )
            entry.grid(
                row=row, column=entry_col, columnspan=entry_span,
                sticky='ew', pady=4
            )

    def reset_separators(self):
        """将当前模式的分隔序列恢复为默认值"""
        mode = self.analysis_mode.get()
        var_map = self.separator_vars.get(mode, {})
        for key, _label, default in self.SEPARATOR_SPECS.get(mode, []):
            if key in var_map:
                var_map[key].set(default)
        self._render_separator_fields()
        name = self.MODE_INFO.get(mode, ("", ""))[0]
        self.log(f"已恢复 {name} 模式的默认分隔序列")

    def _validate_separators(self):
        """校验当前模式的分隔序列，返回规范化后的字典；失败时返回 None"""
        mode = self.analysis_mode.get()
        specs = self.SEPARATOR_SPECS.get(mode, [])
        values = self._get_separator_values(mode)
        var_map = getattr(self, "separator_vars", {}).get(mode, {})
        errors = []
        allowed = set("ACGT")

        for key, label, _default in specs:
            value = values.get(key, "")
            if not value:
                errors.append(f"{label}不能为空")
                continue
            invalid = sorted(set(value) - allowed)
            if invalid:
                errors.append(f"{label}包含非 A/C/G/T 字符：{''.join(invalid)}")
            if key in var_map:
                var_map[key].set(value)

        if mode == "scfv":
            if len(set(values.values())) < len(specs):
                errors.append("三个 Linker 序列不能重复")
        elif mode == "fab":
            seen_values = {}
            for key, label, _default in specs:
                value = values.get(key, "")
                if value in seen_values:
                    errors.append(f"{label}与{seen_values[value]}序列不能相同")
                else:
                    seen_values[value] = label
        elif mode == "vhh":
            if values.get("upstream") == values.get("downstream"):
                errors.append("VHH 上游与下游位点序列不能相同")

        if not errors:
            pairs = []
            if mode == "scfv":
                pairs = [
                    ("Linker 1", "linker1", "Linker 2", "linker2"),
                    ("Linker 1", "linker1", "Linker 3", "linker3"),
                    ("Linker 2", "linker2", "Linker 3", "linker3"),
                ]
            elif mode == "fab":
                items = [(label, key) for key, label, _default in specs]
                for i in range(len(items)):
                    for j in range(i + 1, len(items)):
                        pairs.append((items[i][0], items[i][1], items[j][0], items[j][1]))
            elif mode == "vhh":
                pairs = [("上游位点", "upstream", "下游位点", "downstream")]

            for label_a, key_a, label_b, key_b in pairs:
                value_a = values.get(key_a, "")
                value_b = values.get(key_b, "")
                if value_a and value_b and value_a != value_b and (
                        value_a in value_b or value_b in value_a):
                    errors.append(f"{label_a}与{label_b}不能互为子串")

        if errors:
            message = "\n".join(errors)
            messagebox.showerror("分隔序列错误", message)
            self.log("分隔序列校验失败：" + "；".join(errors), "ERROR")
            return None

        defaults = {key: default for key, _label, default in specs}
        if values != defaults:
            self.log("已应用自定义分隔序列")
        return values

    def _on_mode_changed(self):
        """分析模式切换后更新说明、标题栏、物种选项与分隔序列"""
        mode = self.analysis_mode.get()
        name = self.MODE_INFO.get(mode, ("scFv", ""))[0]
        self.mode_label.config(text=self.MODE_DESCRIPTIONS.get(mode, ""))
        if hasattr(self, "header_mode_label"):
            self.header_mode_label.config(text=f"当前模式：{name}")
        if hasattr(self, "separator_frame"):
            self.separator_frame.config(text=f"分隔序列 - {name}（可自定义）")
        self._render_separator_fields()
        self._update_species_options()
        self._invalidate_workflow()
        self.log(f"分析模式已切换为: {name}")

    def create_settings_tab(self):
        tab = ttk.Frame(self.notebook, padding=(16, 16, 16, 16))
        self.notebook.add(tab, text="参数设置")
        tab.columnconfigure(0, weight=1)
        tab.columnconfigure(1, weight=1)
        tab.rowconfigure(1, weight=1)

        params_frame = ttk.LabelFrame(tab, text="分析参数", padding=(14, 12))
        params_frame.grid(row=0, column=0, columnspan=2, sticky='ew', pady=(0, 12))
        params_frame.columnconfigure(1, weight=1)
        params_frame.columnconfigure(3, weight=1)

        ttk.Label(params_frame, text="物种").grid(row=0, column=0, sticky='w', padx=(0, 10), pady=4)
        self.species_combo = ttk.Combobox(
            params_frame, textvariable=self.current_species,
            values=Config.MODE_SPECIES['scfv'], state='readonly', width=22
        )
        self.species_combo.grid(row=0, column=1, sticky='w', pady=4)
        self.species_combo.bind('<<ComboboxSelected>>', self.on_species_changed)

        ttk.Label(params_frame, text="编号方案").grid(row=0, column=2, sticky='w', padx=(24, 10), pady=4)
        self.scheme_combo = ttk.Combobox(
            params_frame, textvariable=self.current_scheme,
            values=["imgt", "kabat"], state='readonly', width=22
        )
        self.scheme_combo.grid(row=0, column=3, sticky='w', pady=4)
        self.scheme_combo.bind('<<ComboboxSelected>>', self.on_scheme_changed)

        format_frame = ttk.Frame(params_frame)
        format_frame.grid(row=1, column=0, columnspan=4, sticky='ew', pady=(6, 0))
        ttk.Label(
            format_frame,
            text="输入格式（默认 .seq，勾选后同时支持 .gz）：",
            style='Muted.TLabel'
        ).grid(row=0, column=0, sticky='w', padx=(0, 12))
        self.input_format_checks = []
        for i, ext in enumerate(Config.INPUT_FORMAT_OPTIONS):
            check = ttk.Checkbutton(
                format_frame, text=ext, variable=self.input_format_vars[ext],
                command=self._on_input_format_changed
            )
            check.grid(row=0, column=i + 1, sticky='w', padx=(0, 8))
            self.input_format_checks.append(check)
        format_frame.columnconfigure(len(Config.INPUT_FORMAT_OPTIONS) + 1, weight=1)

        mode_name = self.MODE_INFO.get(self.analysis_mode.get(), ("scFv", ""))[0]
        self.separator_frame = ttk.LabelFrame(
            params_frame, text=f"分隔序列 - {mode_name}（可自定义）", padding=(12, 10)
        )
        self.separator_frame.grid(row=2, column=0, columnspan=4, sticky='ew', pady=(10, 0))
        self.separator_frame.columnconfigure(0, weight=1)

        sep_header = ttk.Frame(self.separator_frame)
        sep_header.grid(row=0, column=0, sticky='ew', pady=(0, 6))
        sep_header.columnconfigure(0, weight=1)
        ttk.Label(
            sep_header,
            text="以下序列用于定位目标序列，默认为当前模式的标准序列；可手工修改，仅允许 A/C/G/T。",
            style='Muted.TLabel', wraplength=760, justify='left'
        ).grid(row=0, column=0, sticky='w')
        self.reset_separators_btn = ttk.Button(
            sep_header, text="恢复当前模式默认序列", command=self.reset_separators
        )
        self.reset_separators_btn.grid(row=0, column=1, sticky='e', padx=(12, 0))

        self.separator_fields_frame = ttk.Frame(self.separator_frame)
        self.separator_fields_frame.grid(row=1, column=0, sticky='ew')
        self._render_separator_fields()

        dir_frame = ttk.LabelFrame(tab, text="目录设置", padding=(14, 12))
        dir_frame.grid(row=1, column=0, sticky='nsew', padx=(0, 6))
        dir_frame.columnconfigure(1, weight=1)

        dir_info = [
            ("分析目录", Config.ANALYSE_DIR),
            ("数据库目录", Config.DATABASE_DIR),
            ("程序目录", Config.BIN_DIR),
            ("可选文件目录", Config.OPTIONAL_FILE_DIR)
        ]
        for i, (label, path) in enumerate(dir_info):
            ttk.Label(dir_frame, text=label).grid(row=i, column=0, sticky='w', padx=(0, 10), pady=4)
            path_label = ttk.Label(dir_frame, text=str(path), style='Muted.TLabel')
            path_label.grid(row=i, column=1, sticky='w', pady=4)
            path_label.bind("<Enter>", lambda e, p=str(path): self.show_tooltip(e, p))

        check_frame = ttk.LabelFrame(tab, text="外部工具", padding=(14, 12))
        check_frame.grid(row=1, column=1, sticky='nsew', padx=(6, 0))
        check_frame.columnconfigure(0, weight=1)
        check_frame.rowconfigure(2, weight=1)

        ttk.Label(
            check_frame, text="检查 clustalw2、seqkit 与 igblastn 是否已正确配置。",
            style='Muted.TLabel', wraplength=420, justify='left'
        ).grid(row=0, column=0, sticky='w', pady=(0, 8))
        self.check_tools_btn = ttk.Button(check_frame, text="检查外部工具", command=self.check_tools)
        self.check_tools_btn.grid(row=1, column=0, sticky='w', pady=(0, 8))
        self.tool_status_text = scrolledtext.ScrolledText(
            check_frame, width=40, height=9, bg=self.log_bg, fg=self.text_color,
            relief=tk.FLAT, borderwidth=1, font=(self.mono_font, 9), padx=8, pady=6
        )
        self.tool_status_text.grid(row=2, column=0, sticky='nsew')

    def create_results_tab(self):
        tab = ttk.Frame(self.notebook, padding=(12, 12, 12, 12))
        self.notebook.add(tab, text="结果查看")
        tab.columnconfigure(0, weight=1)
        tab.rowconfigure(0, weight=1)

        self.results_paned = ttk.PanedWindow(tab, orient=tk.VERTICAL)
        paned = self.results_paned
        paned.grid(row=0, column=0, sticky='nsew')

        results_frame = ttk.LabelFrame(paned, text="生成的结果文件", padding=(12, 10))
        paned.add(results_frame, weight=2)
        results_frame.columnconfigure(0, weight=1)
        results_frame.rowconfigure(1, weight=1)

        toolbar = ttk.Frame(results_frame)
        toolbar.grid(row=0, column=0, sticky='ew', pady=(0, 8))
        self.results_count_label = ttk.Label(toolbar, text="共 0 个结果文件", style='Muted.TLabel')
        self.results_count_label.pack(side=tk.LEFT)
        ttk.Label(toolbar, text="双击文件可快速预览", style='Muted.TLabel').pack(side=tk.LEFT, padx=(16, 0))

        tree_box = ttk.Frame(results_frame)
        tree_box.grid(row=1, column=0, sticky='nsew')
        tree_box.columnconfigure(0, weight=1)
        tree_box.rowconfigure(0, weight=1)

        columns = ("文件名", "大小", "修改时间")
        self.results_tree = ttk.Treeview(tree_box, columns=columns, show="headings", height=8)
        self.results_tree.heading("文件名", text="文件名")
        self.results_tree.heading("大小", text="大小")
        self.results_tree.heading("修改时间", text="修改时间")
        self.results_tree.column("文件名", width=320, minwidth=180, stretch=True)
        self.results_tree.column("大小", width=90, minwidth=70)
        self.results_tree.column("修改时间", width=150, minwidth=120)
        scrollbar = ttk.Scrollbar(tree_box, orient=tk.VERTICAL, command=self.results_tree.yview)
        self.results_tree.configure(yscrollcommand=scrollbar.set)
        self.results_tree.grid(row=0, column=0, sticky='nsew')
        scrollbar.grid(row=0, column=1, sticky='ns')
        self.results_tree.tag_configure('odd', background="#f8fafc")
        self.results_tree.tag_configure('even', background=self.panel_color)
        self.results_tree.bind("<Double-1>", self._on_result_double_click)

        button_frame = ttk.Frame(results_frame)
        button_frame.grid(row=2, column=0, sticky='ew', pady=(8, 0))
        refresh_btn = ttk.Button(button_frame, text="刷新列表", command=self.refresh_results)
        refresh_btn.pack(side=tk.LEFT, padx=(0, 6))
        preview_btn = ttk.Button(button_frame, text="预览文件", command=self.preview_selected_file)
        preview_btn.pack(side=tk.LEFT)

        preview_frame = ttk.LabelFrame(paned, text="文件预览", padding=(12, 10))
        paned.add(preview_frame, weight=1)
        preview_frame.columnconfigure(0, weight=1)
        preview_frame.rowconfigure(0, weight=1)
        self.preview_text = scrolledtext.ScrolledText(
            preview_frame, width=40, height=10, bg=self.log_bg, fg=self.text_color,
            relief=tk.FLAT, borderwidth=1, font=(self.mono_font, 9), padx=8, pady=6
        )
        self.preview_text.grid(row=0, column=0, sticky='nsew')

    def _on_result_double_click(self, event):
        if self.results_tree.identify_region(event.x, event.y) == "cell":
            self.preview_selected_file()

    def create_log_tab(self):
        tab = ttk.Frame(self.notebook, padding=(12, 12, 12, 12))
        self.notebook.add(tab, text="运行日志")
        tab.columnconfigure(0, weight=1)
        tab.rowconfigure(1, weight=1)

        toolbar = ttk.Frame(tab)
        toolbar.grid(row=0, column=0, sticky='ew', pady=(0, 8))
        ttk.Label(toolbar, text="程序运行日志", style='Section.TLabel').pack(side=tk.LEFT)
        self.save_log_btn = ttk.Button(toolbar, text="保存日志", command=self.save_log)
        self.save_log_btn.pack(side=tk.RIGHT)
        self.clear_log_btn = ttk.Button(toolbar, text="清空日志", command=self.clear_log)
        self.clear_log_btn.pack(side=tk.RIGHT, padx=(0, 6))

        self.log_text = scrolledtext.ScrolledText(
            tab, width=40, bg=self.log_bg, fg=self.text_color,
            relief=tk.FLAT, borderwidth=1, font=(self.mono_font, 9), padx=8, pady=6
        )
        self.log_text.grid(row=1, column=0, sticky='nsew')
        self.log_text.tag_config("error", foreground=self.error_color)
        self.log_text.tag_config("warning", foreground=self.warning_color)
        self.log_text.tag_config("info", foreground=self.text_color)

    def show_tooltip(self, event, text):
        if getattr(self, "_tooltip", None) is not None:
            try:
                self._tooltip.destroy()
            except tk.TclError:
                pass
        tooltip = tk.Toplevel(self.root)
        tooltip.wm_overrideredirect(True)
        tooltip.wm_geometry(f"+{event.x_root + 12}+{event.y_root + 14}")
        tk.Label(
            tooltip, text=text, bg="#1f2937", fg="#ffffff",
            padx=8, pady=4, font=(self.font_family, 9), justify='left'
        ).pack()
        self._tooltip = tooltip
        tooltip.after(2500, tooltip.destroy)

    def setup_directories(self):
        try:
            Config.ANALYSE_DIR.mkdir(parents=True, exist_ok=True)
            Config.DATABASE_DIR.mkdir(parents=True, exist_ok=True)
            Config.BIN_DIR.mkdir(parents=True, exist_ok=True)
            Config.OPTIONAL_FILE_DIR.mkdir(parents=True, exist_ok=True)
            self.log("目录初始化完成")
        except Exception as e:
            self.log(f"目录初始化失败: {str(e)}", "ERROR")

    @staticmethod
    def _is_generated_output(path):
        """只归档程序生成的结果，避免误移动用户输入文件"""
        name = path.name
        exact = {
            'Error_seq.txt',
            'VH_nt_seq.txt', 'VL_nt_seq.txt', 'scfv_nt_seq.txt', 'VHH_nt_seq.txt',
            'VH_aa_seq.txt', 'VL_aa_seq.txt', 'scfv_aa_seq.txt', 'VHH_aa_seq.txt',
            'run_manifest.json',
        }
        if name in exact:
            return True
        patterns = (
            '*_unique.txt', '*_pureunique.txt', '*scfvunique.txt', '*fabunique.txt',
            'duplicates.detail_*', '*.aln', '*.dnd', '*-Seq_unset.csv',
            '*-Seq_set.xlsx', '*.tmp',
        )
        return any(fnmatch.fnmatch(name, pattern) for pattern in patterns)

    def _archive_old_outputs(self):
        """新流程开始前归档旧结果，避免新旧文件混用"""
        if not Config.ANALYSE_DIR.exists():
            return
        old_files = [
            path for path in Config.ANALYSE_DIR.iterdir()
            if path.is_file() and self._is_generated_output(path)
        ]
        if not old_files:
            return

        archive_dir = Config.ANALYSE_DIR / 'archive' / datetime.now().strftime('%Y%m%d_%H%M%S')
        archive_dir.mkdir(parents=True, exist_ok=True)
        moved = 0
        for path in old_files:
            target = archive_dir / path.name
            counter = 1
            while target.exists():
                target = archive_dir / f"{path.stem}_{counter}{path.suffix}"
                counter += 1
            try:
                os.replace(path, target)
                moved += 1
            except OSError as e:
                self.log(f"警告: 无法归档旧文件 {path.name}: {str(e)}", "WARNING")
        if moved:
            self.log(f"已将 {moved} 个旧结果归档到: {archive_dir}")
            self.refresh_results()
        else:
            try:
                archive_dir.rmdir()
            except OSError:
                pass

    def select_seq_files(self):
        patterns = []
        for ext in self._selected_input_formats():
            patterns.extend([f"*{ext}", f"*{ext}.gz"])
        files = filedialog.askopenfilenames(
            title="选择序列文件",
            filetypes=[
                ("Sequence files", tuple(patterns)),
                ("All files", "*.*"),
            ]
        )
        if files:
            self.seq_files = list(files)
            self.update_file_listbox()
            self.file_count_label.config(text=f"已选择 {len(files)} 个文件")
            self.log(f"已选择 {len(files)} 个序列文件")
            self.folder_info_label.config(text="")
            self._invalidate_workflow()

    def select_seq_folder(self):
        folder_path = filedialog.askdirectory(title="选择包含序列文件的文件夹")
        if folder_path:
            self.find_seq_files_in_folder(folder_path)

    def find_seq_files_in_folder(self, folder_path):
        try:
            seq_files = []
            folder_name = os.path.basename(folder_path)

            if self.include_subfolders.get():
                for root, _dirs, files in os.walk(folder_path):
                    for file in files:
                        if self._is_supported_input(file):
                            seq_files.append(os.path.join(root, file))
                self.log(f"扫描文件夹 '{folder_name}'（包含子文件夹）")
            else:
                for file in os.listdir(folder_path):
                    if self._is_supported_input(file):
                        seq_files.append(os.path.join(folder_path, file))
                self.log(f"扫描文件夹 '{folder_name}'")

            if seq_files:
                self.seq_files = sorted(seq_files)
                self.update_file_listbox()

                folder_count = len(set(os.path.dirname(f) for f in seq_files))
                if self.include_subfolders.get() and folder_count > 1:
                    self.folder_info_label.config(
                        text=f"从 {folder_count} 个文件夹中找到 {len(seq_files)} 个序列文件"
                    )
                else:
                    self.folder_info_label.config(
                        text=f"从文件夹 '{folder_name}' 中找到 {len(seq_files)} 个序列文件"
                    )

                self.file_count_label.config(text=f"已选择 {len(seq_files)} 个文件")
                self.log(f"从文件夹中找到 {len(seq_files)} 个序列文件")
                self._invalidate_workflow()
            else:
                self.log(f"在文件夹 '{folder_name}' 中未找到支持的序列文件", "WARNING")
                messagebox.showwarning("警告", f"在文件夹 '{folder_name}' 中未找到支持的序列文件")

        except Exception as e:
            self.log(f"扫描文件夹时出错: {str(e)}", "ERROR")
            messagebox.showerror("错误", f"扫描文件夹时出错: {str(e)}")

    def update_file_listbox(self):
        """刷新文件列表"""
        self.file_tree.delete(*self.file_tree.get_children())
        for index, file in enumerate(self.seq_files):
            folder = os.path.dirname(file)
            try:
                relative = os.path.relpath(folder)
                if len(relative) < len(folder):
                    folder = relative
            except ValueError:
                pass
            tag = 'odd' if index % 2 else 'even'
            self.file_tree.insert("", tk.END, values=(os.path.basename(file), folder), tags=(tag,))

    def remove_selected_files(self):
        """移除列表中选中的文件"""
        selection = self.file_tree.selection()
        if not selection:
            messagebox.showinfo("提示", "请先选择要移除的文件")
            return
        indexes = sorted((self.file_tree.index(item) for item in selection), reverse=True)
        for index in indexes:
            del self.seq_files[index]
        self.update_file_listbox()
        self.file_count_label.config(text=f"已选择 {len(self.seq_files)} 个文件")
        self.log(f"已移除 {len(indexes)} 个文件")
        self._invalidate_workflow()
        if not self.seq_files:
            self.folder_info_label.config(text="")

    def clear_file_list(self):
        self.seq_files = []
        self.file_tree.delete(*self.file_tree.get_children())
        self.file_count_label.config(text="已选择 0 个文件")
        self.folder_info_label.config(text="")
        self.log("已清空文件列表")
        self._invalidate_workflow()

    # ==================== 优化后的核心功能 ====================

    def safe_subprocess_run(self, cmd, shell=False, timeout=None):
        """统一处理子进程调用，避免 Windows 下弹出控制台窗口"""
        if os.name == 'nt':
            startupinfo_class = getattr(subprocess, "STARTUPINFO", None)
            startupinfo = None
            if startupinfo_class is not None:
                startupinfo = startupinfo_class()
                startupinfo.dwFlags |= getattr(subprocess, "STARTF_USESHOWWINDOW", 0)
                startupinfo.wShowWindow = getattr(subprocess, "SW_HIDE", 0)
            return subprocess.run(
                cmd,
                shell=shell,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                stdin=subprocess.DEVNULL,
                text=True,
                timeout=timeout,
                startupinfo=startupinfo,
                creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0)
            )
        return subprocess.run(
            cmd, shell=shell, capture_output=True, text=True,
            timeout=timeout
        )

    def extract_sequences(self):
        if not self.seq_files:
            messagebox.showwarning("警告", "请先选择序列文件")
            return
        separators = self._validate_separators()
        if separators is None:
            return
        mode = self.analysis_mode.get()
        self._archive_old_outputs()
        self._invalidate_workflow()
        self._write_run_manifest(
            mode, separators, self.current_species.get(),
            self.current_scheme.get(), self.skip_alignment.get())
        self._active_mode = mode
        self._active_separators = separators
        self.set_processing_state(True)
        thread = threading.Thread(
            target=self._extract_sequences_thread, args=(mode, separators))
        thread.daemon = True
        thread.start()

    def _extract_sequences_thread(self, mode=None, separators=None):
        status = "failed"
        try:
            error_path = Config.ANALYSE_DIR / 'Error_seq.txt'
            if error_path.exists():
                error_path.unlink()
                self.log("已清理旧的错误文件")
            mode = mode or getattr(self, "_active_mode", None) or self.analysis_mode.get()
            mode_name = self.MODE_INFO.get(mode, ("scFv", ""))[0]
            self.log(f"开始提取序列（{mode_name}模式）...")
            success_count, fail_count = self._extract_current_sequences(mode, separators)

            if success_count > 0:
                status = "completed_with_warnings" if fail_count else "completed"
                self.log(f"序列提取完成: 成功 {success_count} 条, 失败 {fail_count} 条")
                self._post_ui(self._mark_extracted)
                self._post_ui(lambda: messagebox.showinfo(
                    "完成", f"序列提取完成！\n成功: {success_count} 条\n失败: {fail_count} 条"))
                if error_path.exists():
                    self.log(f"错误文件保存在: {error_path.absolute()}")
            else:
                self.log("未能成功提取任何序列", "ERROR")
                self._post_ui(lambda: messagebox.showwarning("警告", "未能成功提取任何序列"))
        except Exception as e:
            self.log(f"序列提取失败: {str(e)}", "ERROR")
            self.log(traceback.format_exc(), "ERROR")
        finally:
            self._update_run_manifest(status)
            self._post_ui(lambda: self.set_processing_state(False))

    def _extract_current_sequences(self, mode=None, separators=None):
        """按当前模式提取序列"""
        mode = mode or getattr(self, "_active_mode", None) or self.analysis_mode.get()
        if separators is None:
            separators = getattr(self, "_active_separators", None)
        if mode == "fab":
            return self.extract_fab_sequences(self.seq_files, separators)
        if mode == "vhh":
            return self.extract_vhh_sequences(self.seq_files, separators)
        return self.extract_scfv_sequences(self.seq_files, separators)

    def _selected_input_formats(self):
        vars = getattr(self, "input_format_vars", None)
        if not vars:
            return tuple(Config.DEFAULT_INPUT_FORMATS)
        selected = tuple(
            ext for ext in Config.INPUT_FORMAT_OPTIONS if vars[ext].get()
        )
        return selected or tuple(Config.DEFAULT_INPUT_FORMATS)

    def _supported_input_suffixes(self):
        suffixes = []
        for ext in self._selected_input_formats():
            suffixes.append(ext)
            suffixes.append(ext + '.gz')
        return tuple(suffixes)

    def _is_supported_input(self, filename):
        lower = str(filename).lower()
        return any(lower.endswith(suffix) for suffix in self._supported_input_suffixes())

    def _on_input_format_changed(self):
        vars = getattr(self, "input_format_vars", None)
        if not vars:
            return
        if not any(var.get() for var in vars.values()):
            vars['.seq'].set(True)
            messagebox.showinfo("提示", "至少选择一种输入格式，已恢复 .seq")
        before = len(self.seq_files)
        self.seq_files = [f for f in self.seq_files if self._is_supported_input(f)]
        removed = before - len(self.seq_files)
        if removed:
            self.update_file_listbox()
            self.file_count_label.config(text=f"已选择 {len(self.seq_files)} 个文件")
            if not self.seq_files:
                self.folder_info_label.config(text="")
            self.log(f"输入格式变更，已移除 {removed} 个不再支持的文件")
        self._invalidate_workflow()
        self.log("输入格式: " + ", ".join(self._selected_input_formats()))

    @staticmethod
    def _iter_fasta_records(filename):
        """按 FASTA 读取普通或 gzip 压缩文件"""
        if str(filename).lower().endswith('.gz'):
            with gzip.open(filename, 'rt', encoding='utf-8') as handle:
                yield from SeqIO.parse(handle, 'fasta')
        else:
            yield from SeqIO.parse(filename, 'fasta')

    def _resolve_tool(self, name):
        """优先使用 bin/程序目录下的工具，其次使用 PATH；统一返回绝对路径"""
        candidates = [
            Config.BIN_DIR / name,
            Config.BIN_DIR / f"{name}.exe",
            Config.BASE_DIR / name,
            Config.BASE_DIR / f"{name}.exe",
        ]
        for candidate in candidates:
            if candidate.is_file():
                return str(candidate.resolve())
        found = shutil.which(name) or shutil.which(f"{name}.exe")
        if found:
            return str(Path(found).resolve())
        return None

    @staticmethod
    def _find_all_positions(sequence, pattern):
        """返回 pattern 在 sequence 中的所有起始位置"""
        positions = []
        start = 0
        while True:
            pos = sequence.find(pattern, start)
            if pos == -1:
                break
            positions.append(pos)
            start = pos + 1
        return positions

    def _make_unique_id(self, seq_id, used_ids):
        """为重复的 FASTA ID 生成稳定且唯一的后缀"""
        if seq_id not in used_ids:
            used_ids.add(seq_id)
            return seq_id
        original_id = seq_id
        suffix = 2
        while f"{original_id}__dup{suffix}" in used_ids:
            suffix += 1
        unique_id = f"{original_id}__dup{suffix}"
        used_ids.add(unique_id)
        self.log(f"    警告: 序列ID重复 ({original_id})，已重命名为 {unique_id}", "WARNING")
        return unique_id

    def extract_scfv_sequences(self, seq_files, separators=None):
        """提取VH/VL/scFv序列，返回字典缓存，同时原子写入文件"""
        vh_nt_path = Config.ANALYSE_DIR / 'VH_nt_seq.txt'
        vl_nt_path = Config.ANALYSE_DIR / 'VL_nt_seq.txt'
        error_path = Config.ANALYSE_DIR / 'Error_seq.txt'
        scfv_nt_path = Config.ANALYSE_DIR / 'scfv_nt_seq.txt'

        Config.ANALYSE_DIR.mkdir(parents=True, exist_ok=True)

        self.vh_nt_dict.clear()
        self.vl_nt_dict.clear()
        self.scfv_nt_dict.clear()
        self.vhh_nt_dict.clear()

        success_count = 0
        fail_count = 0
        total_files = len(seq_files)
        used_ids = set()

        if separators is None:
            separators = getattr(self, "_active_separators", None)
        if not separators:
            separators = self._get_separator_values("scfv")
        linker1 = separators['linker1']
        linker2 = separators['linker2']
        linker3 = separators['linker3']
        linker1_len = len(linker1)
        linker2_len = len(linker2)
        linker_labels = {
            'linker1': 'Linker 1',
            'linker2': 'Linker 2',
            'linker3': 'Linker 3',
        }

        vh_tmp = vh_nt_path.with_name(vh_nt_path.name + '.tmp')
        vl_tmp = vl_nt_path.with_name(vl_nt_path.name + '.tmp')
        error_tmp = error_path.with_name(error_path.name + '.tmp')
        scfv_tmp = scfv_nt_path.with_name(scfv_nt_path.name + '.tmp')
        tmp_paths = (vh_tmp, vl_tmp, error_tmp, scfv_tmp)

        try:
            with open(vh_tmp, 'w', encoding='utf-8', buffering=1024 * 1024) as vh_nt, \
                    open(vl_tmp, 'w', encoding='utf-8', buffering=1024 * 1024) as vl_nt, \
                    open(error_tmp, 'w', encoding='utf-8', buffering=1024 * 1024) as error_seq, \
                    open(scfv_tmp, 'w', encoding='utf-8', buffering=1024 * 1024) as scfv_nt:
                last_progress = -1
                for i, filename in enumerate(seq_files):
                    file_basename = os.path.basename(filename)
                    if not self._is_supported_input(filename):
                        self.log(f"跳过不支持的文件: {file_basename}", "WARNING")
                        continue
                    self.log(f"处理文件: {file_basename} ({i + 1}/{total_files})")

                    try:
                        for seq_record in self._iter_fasta_records(filename):
                            seq_id = seq_record.id
                            seq_str = str(seq_record.seq).upper()

                            missing = []
                            for key, label in linker_labels.items():
                                if separators[key] not in seq_str:
                                    missing.append(label)
                            if missing:
                                error_seq.write(f'>{seq_id}\n缺少linker序列: {missing}\n')
                                self.log(f"    错误: 缺少linker: {missing}")
                                fail_count += 1
                                continue

                            l1_positions = self._find_all_positions(seq_str, linker1)
                            l2_positions = self._find_all_positions(seq_str, linker2)
                            l3_positions = self._find_all_positions(seq_str, linker3)
                            best = None
                            for p1 in l1_positions:
                                for p2 in l2_positions:
                                    if p2 < p1 + linker1_len:
                                        continue
                                    vh_len = p2 - (p1 + linker1_len)
                                    if vh_len < 50:
                                        continue
                                    for p3 in l3_positions:
                                        if p3 < p2 + linker2_len:
                                            continue
                                        vl_len = p3 - (p2 + linker2_len)
                                        if vl_len < 50:
                                            continue
                                        span = p3 - p1
                                        if best is None or span < best[0]:
                                            best = (span, p1, p2, p3)

                            if best is None:
                                error_seq.write(
                                    f'>{seq_id}\n未找到合法的linker组合: '
                                    f'L1={l1_positions}, L2={l2_positions}, L3={l3_positions}\n')
                                self.log("    错误: 未找到合法的linker组合", "WARNING")
                                fail_count += 1
                                continue

                            _span, linker1_pos, linker2_pos, linker3_pos = best
                            vh_seq = seq_str[linker1_pos + linker1_len:linker2_pos]
                            vl_seq = seq_str[linker2_pos + linker2_len:linker3_pos]
                            scfv_seq = seq_str[linker1_pos + linker1_len:linker3_pos]

                            seq_id = self._make_unique_id(seq_id, used_ids)
                            vh_nt.write(f'>{seq_id}\n{vh_seq}\n')
                            vl_nt.write(f'>{seq_id}\n{vl_seq}\n')
                            scfv_nt.write(f'>{seq_id}\n{scfv_seq}\n')

                            self.vh_nt_dict[seq_id] = vh_seq
                            self.vl_nt_dict[seq_id] = vl_seq
                            self.scfv_nt_dict[seq_id] = scfv_seq

                            success_count += 1

                    except Exception as e:
                        error_seq.write(f'处理文件{file_basename}时出错: {str(e)}\n')
                        self.log(f"    异常错误: {str(e)}", "ERROR")
                        fail_count += 1

                    progress = int((i + 1) / total_files * 100)
                    if progress != last_progress:
                        last_progress = progress
                        self._post_ui(lambda p=progress: self.progress_var.set(p))

            os.replace(vh_tmp, vh_nt_path)
            os.replace(vl_tmp, vl_nt_path)
            os.replace(error_tmp, error_path)
            os.replace(scfv_tmp, scfv_nt_path)
        finally:
            for tmp in tmp_paths:
                if tmp.exists():
                    try:
                        tmp.unlink()
                    except OSError:
                        pass

        self._post_ui(lambda: self.progress_var.set(100))
        return success_count, fail_count

    def extract_fab_sequences(self, seq_files, separators=None):
        """提取Fab的VH/VL序列，根据特征上下游序列定位，分别原子写入文件"""
        vh_nt_path = Config.ANALYSE_DIR / 'VH_nt_seq.txt'
        vl_nt_path = Config.ANALYSE_DIR / 'VL_nt_seq.txt'
        error_path = Config.ANALYSE_DIR / 'Error_seq.txt'

        Config.ANALYSE_DIR.mkdir(parents=True, exist_ok=True)

        self.vh_nt_dict.clear()
        self.vl_nt_dict.clear()
        self.scfv_nt_dict.clear()
        self.vhh_nt_dict.clear()

        success_count = 0
        fail_count = 0
        total_files = len(seq_files)
        used_vh_ids = set()
        used_vl_ids = set()

        if separators is None:
            separators = getattr(self, "_active_separators", None)
        if not separators:
            separators = self._get_separator_values("fab")
        flank_patterns = {
            'VH': {
                'upstream': separators['vh_upstream'],
                'downstream': separators['vh_downstream'],
            },
            'VL': {
                'upstream': separators['vl_upstream'],
                'downstream': separators['vl_downstream'],
            },
        }

        vh_tmp = vh_nt_path.with_name(vh_nt_path.name + '.tmp')
        vl_tmp = vl_nt_path.with_name(vl_nt_path.name + '.tmp')
        error_tmp = error_path.with_name(error_path.name + '.tmp')
        tmp_paths = (vh_tmp, vl_tmp, error_tmp)

        try:
            with open(vh_tmp, 'w', encoding='utf-8', buffering=1024 * 1024) as vh_nt, \
                    open(vl_tmp, 'w', encoding='utf-8', buffering=1024 * 1024) as vl_nt, \
                    open(error_tmp, 'w', encoding='utf-8', buffering=1024 * 1024) as error_seq:
                last_progress = -1
                for i, filename in enumerate(seq_files):
                    file_basename = os.path.basename(filename)
                    if not self._is_supported_input(filename):
                        self.log(f"跳过不支持的文件: {file_basename}", "WARNING")
                        continue
                    self.log(f"处理文件: {file_basename} ({i + 1}/{total_files})")

                    try:
                        for seq_record in self._iter_fasta_records(filename):
                            seq_id = seq_record.id
                            seq_str = str(seq_record.seq).upper()

                            found_chain = False
                            for chain_type, flank in flank_patterns.items():
                                upstream = flank['upstream']
                                downstream = flank['downstream']
                                up_pos = -1
                                down_pos = -1
                                best_len = None
                                for candidate in self._find_all_positions(seq_str, upstream):
                                    candidate_down = seq_str.find(
                                        downstream, candidate + len(upstream))
                                    if candidate_down == -1:
                                        continue
                                    candidate_len = candidate_down - (candidate + len(upstream))
                                    if candidate_len < 50:
                                        continue
                                    if best_len is None or candidate_len < best_len:
                                        best_len = candidate_len
                                        up_pos, down_pos = candidate, candidate_down

                                if up_pos == -1 or down_pos == -1:
                                    continue

                                chain_seq = seq_str[up_pos + len(upstream):down_pos]
                                if chain_type == 'VH':
                                    chain_id = self._make_unique_id(seq_id, used_vh_ids)
                                    vh_nt.write(f'>{chain_id}\n{chain_seq}\n')
                                    self.vh_nt_dict[chain_id] = chain_seq
                                else:
                                    chain_id = self._make_unique_id(seq_id, used_vl_ids)
                                    vl_nt.write(f'>{chain_id}\n{chain_seq}\n')
                                    self.vl_nt_dict[chain_id] = chain_seq

                                found_chain = True

                            if found_chain:
                                success_count += 1
                            else:
                                error_seq.write(f'>{seq_id}\n未找到VH或VL特征上下游序列\n')
                                self.log("    错误: 未找到VH/VL特征上下游序列", "WARNING")
                                fail_count += 1

                    except Exception as e:
                        error_seq.write(f'处理文件{file_basename}时出错: {str(e)}\n')
                        self.log(f"    异常错误: {str(e)}", "ERROR")
                        fail_count += 1

                    progress = int((i + 1) / total_files * 100)
                    if progress != last_progress:
                        last_progress = progress
                        self._post_ui(lambda p=progress: self.progress_var.set(p))

            os.replace(vh_tmp, vh_nt_path)
            os.replace(vl_tmp, vl_nt_path)
            os.replace(error_tmp, error_path)
        finally:
            for tmp in tmp_paths:
                if tmp.exists():
                    try:
                        tmp.unlink()
                    except OSError:
                        pass

        self._post_ui(lambda: self.progress_var.set(100))
        return success_count, fail_count

    def extract_vhh_sequences(self, seq_files, separators=None):
        """提取VHH序列（自定义上下游位点之间），原子写入VHH_nt_seq.txt"""
        vhh_nt_path = Config.ANALYSE_DIR / 'VHH_nt_seq.txt'
        error_path = Config.ANALYSE_DIR / 'Error_seq.txt'

        Config.ANALYSE_DIR.mkdir(parents=True, exist_ok=True)

        self.vh_nt_dict.clear()
        self.vl_nt_dict.clear()
        self.scfv_nt_dict.clear()
        self.vhh_nt_dict.clear()

        success_count = 0
        fail_count = 0
        total_files = len(seq_files)
        used_ids = set()

        if separators is None:
            separators = getattr(self, "_active_separators", None)
        if not separators:
            separators = self._get_separator_values("vhh")
        upstream = separators['upstream']
        downstream = separators['downstream']

        vhh_tmp = vhh_nt_path.with_name(vhh_nt_path.name + '.tmp')
        error_tmp = error_path.with_name(error_path.name + '.tmp')
        tmp_paths = (vhh_tmp, error_tmp)

        try:
            with open(vhh_tmp, 'w', encoding='utf-8', buffering=1024 * 1024) as vhh_nt, \
                    open(error_tmp, 'w', encoding='utf-8', buffering=1024 * 1024) as error_seq:
                last_progress = -1
                for i, filename in enumerate(seq_files):
                    file_basename = os.path.basename(filename)
                    if not self._is_supported_input(filename):
                        self.log(f"跳过不支持的文件: {file_basename}", "WARNING")
                        continue
                    self.log(f"处理文件: {file_basename} ({i + 1}/{total_files})")

                    try:
                        for seq_record in self._iter_fasta_records(filename):
                            seq_id = seq_record.id
                            seq_str = str(seq_record.seq).upper()

                            up_pos = -1
                            down_pos = -1
                            best_len = None
                            for candidate in self._find_all_positions(seq_str, upstream):
                                candidate_down = seq_str.find(
                                    downstream, candidate + len(upstream))
                                if candidate_down == -1:
                                    continue
                                candidate_len = candidate_down - (candidate + len(upstream))
                                if candidate_len < 50:
                                    continue
                                if best_len is None or candidate_len < best_len:
                                    best_len = candidate_len
                                    up_pos, down_pos = candidate, candidate_down

                            if up_pos == -1 or down_pos == -1:
                                error_seq.write(f'>{seq_id}\n缺少上游或下游位点序列\n')
                                self.log("    错误: 缺少上游或下游位点序列", "WARNING")
                                fail_count += 1
                                continue

                            vhh_seq = seq_str[up_pos + len(upstream):down_pos]
                            seq_id = self._make_unique_id(seq_id, used_ids)
                            vhh_nt.write(f'>{seq_id}\n{vhh_seq}\n')
                            self.vhh_nt_dict[seq_id] = vhh_seq
                            success_count += 1

                    except Exception as e:
                        error_seq.write(f'处理文件{file_basename}时出错: {str(e)}\n')
                        self.log(f"    异常错误: {str(e)}", "ERROR")
                        fail_count += 1

                    progress = int((i + 1) / total_files * 100)
                    if progress != last_progress:
                        last_progress = progress
                        self._post_ui(lambda p=progress: self.progress_var.set(p))

            os.replace(vhh_tmp, vhh_nt_path)
            os.replace(error_tmp, error_path)
        finally:
            for tmp in tmp_paths:
                if tmp.exists():
                    try:
                        tmp.unlink()
                    except OSError:
                        pass

        self._post_ui(lambda: self.progress_var.set(100))
        return success_count, fail_count

    def translate_sequences(self):
        mode = self.analysis_mode.get()
        state = getattr(self, "workflow_state", None)
        if state is not None:
            state["translated"] = False
            state["deduplicated"] = False
            state["aligned"] = False
        self._refresh_step_buttons()
        self.set_processing_state(True)
        thread = threading.Thread(target=self._translate_sequences_thread, args=(mode,))
        thread.daemon = True
        thread.start()

    def _translate_sequences_thread(self, mode=None):
        try:
            self.log("开始翻译核酸序列为氨基酸序列...")
            if not self._translate_current_sequences(mode):
                self.log("序列翻译未完成，请检查缺失的核酸文件", "ERROR")
                return
            self.log("序列翻译完成")
            self._post_ui(self._mark_translated)
        except Exception as e:
            self.log(f"翻译序列时出错: {str(e)}", "ERROR")
            self.log(traceback.format_exc(), "ERROR")
        finally:
            self._post_ui(lambda: self.set_processing_state(False))

    def _get_translation_files(self, mode=None):
        """当前模式下需要翻译的文件对"""
        if self.is_vhh_mode(mode):
            return [('VHH_nt_seq.txt', 'VHH_aa_seq.txt')]
        if self.is_fab_mode(mode):
            return [
                ('VH_nt_seq.txt', 'VH_aa_seq.txt'),
                ('VL_nt_seq.txt', 'VL_aa_seq.txt')
            ]
        return [
            ('VH_nt_seq.txt', 'VH_aa_seq.txt'),
            ('VL_nt_seq.txt', 'VL_aa_seq.txt'),
            ('scfv_nt_seq.txt', 'scfv_aa_seq.txt')
        ]

    def _translate_current_sequences(self, mode=None):
        """按当前模式翻译序列；逐链处理，空链不会覆盖已有结果"""
        mode = mode or self.analysis_mode.get()
        if self.is_vhh_mode(mode):
            if self.vhh_nt_dict:
                return self._translate_dict_to_file(
                    self.vhh_nt_dict, Config.ANALYSE_DIR / 'VHH_aa_seq.txt')
            return self._translate_nt_file('VHH_nt_seq.txt', 'VHH_aa_seq.txt')

        chain_files = [
            (self.vh_nt_dict, 'VH_nt_seq.txt', 'VH_aa_seq.txt'),
            (self.vl_nt_dict, 'VL_nt_seq.txt', 'VL_aa_seq.txt'),
        ]
        if not self.is_fab_mode(mode):
            chain_files.append((self.scfv_nt_dict, 'scfv_nt_seq.txt', 'scfv_aa_seq.txt'))

        use_cache = any(seq_dict for seq_dict, _in, _out in chain_files)
        success = True
        for seq_dict, input_file, output_file in chain_files:
            if use_cache:
                if seq_dict:
                    ok = self._translate_dict_to_file(seq_dict, Config.ANALYSE_DIR / output_file)
                else:
                    self.log(f"警告: 当前提取结果缺少 {output_file} 对应的核酸序列", "WARNING")
                    ok = False
            else:
                ok = self._translate_nt_file(input_file, output_file)
            success = success and ok
        return success

    def _translate_nt_file(self, input_file, output_file):
        """从核酸文件读取并翻译为氨基酸文件，空输入不覆盖旧结果"""
        input_path = Config.ANALYSE_DIR / input_file
        output_path = Config.ANALYSE_DIR / output_file
        if not input_path.exists():
            self.log(f"警告: 文件 {input_file} 不存在", "WARNING")
            return False

        tmp_path = output_path.with_name(output_path.name + ".tmp")
        count = 0
        non_multiple = 0
        internal_stops = 0
        ambiguous = 0
        try:
            with open(tmp_path, 'w', encoding='utf-8', buffering=1024 * 1024) as out_f:
                for rec in SeqIO.parse(input_path, 'fasta'):
                    aa_seq = str(rec.seq.translate())
                    out_f.write(f'>{rec.id}\n{aa_seq}\n')
                    count += 1
                    if len(rec.seq) % 3:
                        non_multiple += 1
                    if '*' in aa_seq[:-1]:
                        internal_stops += 1
                    if 'X' in aa_seq:
                        ambiguous += 1
            if count == 0:
                self.log(f"警告: {input_file} 为空，未生成 {output_file}", "WARNING")
                return False
            os.replace(tmp_path, output_path)
            self.log(f"已翻译: {input_file} -> {output_file} ({count}条)")
            if non_multiple or internal_stops or ambiguous:
                self.log(
                    f"翻译质控 {output_file}: 非3倍数={non_multiple}, "
                    f"内部终止={internal_stops}, 含X={ambiguous}", "WARNING")
            return True
        finally:
            if tmp_path.exists():
                try:
                    tmp_path.unlink()
                except OSError:
                    pass

    def _translate_dict_to_file(self, seq_dict, output_path):
        """将字典中的核酸序列翻译并原子写入文件"""
        if not seq_dict:
            self.log(f"警告: 没有可翻译的序列，未生成 {output_path.name}", "WARNING")
            return False

        tmp_path = output_path.with_name(output_path.name + ".tmp")
        non_multiple = 0
        internal_stops = 0
        ambiguous = 0
        try:
            with open(tmp_path, 'w', encoding='utf-8', buffering=1024 * 1024) as f:
                for seq_id, nt_seq in seq_dict.items():
                    aa_seq = str(Seq(nt_seq).translate())
                    f.write(f'>{seq_id}\n{aa_seq}\n')
                    if len(nt_seq) % 3:
                        non_multiple += 1
                    if '*' in aa_seq[:-1]:
                        internal_stops += 1
                    if 'X' in aa_seq:
                        ambiguous += 1
            os.replace(tmp_path, output_path)
            self.log(f"已翻译: {output_path.name} ({len(seq_dict)}条)")
            if non_multiple or internal_stops or ambiguous:
                self.log(
                    f"翻译质控 {output_path.name}: 非3倍数={non_multiple}, "
                    f"内部终止={internal_stops}, 含X={ambiguous}", "WARNING")
            return True
        finally:
            if tmp_path.exists():
                try:
                    tmp_path.unlink()
                except OSError:
                    pass

    def remove_duplicates(self):
        mode = self.analysis_mode.get()
        state = getattr(self, "workflow_state", None)
        if state is not None:
            state["deduplicated"] = False
            state["aligned"] = False
        self._refresh_step_buttons()
        self.set_processing_state(True)
        thread = threading.Thread(target=self._remove_duplicates_thread, args=(mode,))
        thread.daemon = True
        thread.start()

    def _remove_duplicates_thread(self, mode=None):
        try:
            self.log("开始序列去重...")
            if self._deduplicate_current_sequences(mode):
                self.log("序列去重完成")
                self._post_ui(self._mark_deduplicated)
            else:
                self.log("序列去重失败", "ERROR")
        except Exception as e:
            self.log(f"去重时出错: {str(e)}", "ERROR")
        finally:
            self._post_ui(lambda: self.set_processing_state(False))

    def _deduplicate_current_sequences(self, mode=None):
        """按当前模式对氨基酸序列去重，并提取对应的去重核酸序列"""
        mode = mode or self.analysis_mode.get()
        if self.is_vhh_mode(mode):
            if not self.run_seqkit_command(
                    "VHH_aa_seq.txt", "VHH_aa_unique.txt", "duplicates.detail_VHH_aa.txt"):
                return False
            return self.extract_vhh_unique_nt()

        if self.is_fab_mode(mode):
            vh_ok = self.run_seqkit_command(
                "VH_aa_seq.txt", "VH_aa_unique.txt", "duplicates.detail_VH_aa.txt")
            vl_ok = self.run_seqkit_command(
                "VL_aa_seq.txt", "VL_aa_unique.txt", "duplicates.detail_VL_aa.txt")
            if not (vh_ok and vl_ok):
                return False
            return self.extract_fab_unique_sequences()

        if not self.run_seqkit_command(
                "scfv_aa_seq.txt", "scfv_aa_unique.txt", "duplicates.detail_scfv_aa.txt"):
            return False
        unique_ok = self.extract_unique_sequences()
        vh_ok = self.run_seqkit_command(
            "VH_aa_scfvunique.txt", "VH_aa_pureunique.txt",
            "duplicates.detail_VH_aa_scfvunique.txt")
        vl_ok = self.run_seqkit_command(
            "VL_aa_scfvunique.txt", "VL_aa_pureunique.txt",
            "duplicates.detail_VL_aa_scfvunique.txt")
        return unique_ok and vh_ok and vl_ok

    def run_seqkit_command(self, input_file, output_file, dup_detail_file):
        input_path = Config.ANALYSE_DIR / input_file
        output_path = Config.ANALYSE_DIR / output_file
        dup_path = Config.ANALYSE_DIR / dup_detail_file
        if not input_path.exists():
            self.log(f"错误: seqkit输入文件不存在: {input_file}", "ERROR")
            return False

        out_tmp = output_path.with_name(output_path.name + ".tmp")
        dup_tmp = dup_path.with_name(dup_path.name + ".tmp")
        for tmp in (out_tmp, dup_tmp):
            if tmp.exists():
                try:
                    tmp.unlink()
                except OSError:
                    pass

        seqkit_exe = self._resolve_tool("seqkit")
        if not seqkit_exe:
            self.log("错误: 未找到seqkit，请放入bin目录或加入PATH", "ERROR")
            return False
        cmd = [
            seqkit_exe, "rmdup", "-s",
            str(input_path),
            "-o", str(out_tmp),
            "-D", str(dup_tmp)
        ]
        try:
            result = self.safe_subprocess_run(cmd, shell=False, timeout=600)
            if result.returncode != 0:
                self.log(f"seqkit命令执行失败: {result.stderr}", "ERROR")
                return False
            if not out_tmp.exists() or out_tmp.stat().st_size == 0:
                self.log(f"错误: seqkit未生成有效输出: {output_file}", "ERROR")
                return False
            try:
                count = sum(1 for _ in SeqIO.parse(out_tmp, "fasta"))
            except Exception as e:
                self.log(f"错误: seqkit输出无法解析: {str(e)}", "ERROR")
                return False
            if count == 0:
                self.log(f"错误: seqkit输出为空: {output_file}", "ERROR")
                return False

            os.replace(out_tmp, output_path)
            if dup_tmp.exists():
                os.replace(dup_tmp, dup_path)
            else:
                dup_path.write_text('', encoding='utf-8')
            self.log(f"seqkit去重完成: {input_file} -> {output_file} ({count}条)")
            return True
        except FileNotFoundError:
            self.log("错误: 未找到seqkit，请确保它在系统路径中", "ERROR")
            return False
        except Exception as e:
            self.log(f"seqkit命令执行出错: {str(e)}", "ERROR")
            return False
        finally:
            for tmp in (out_tmp, dup_tmp):
                if tmp.exists():
                    try:
                        tmp.unlink()
                    except OSError:
                        pass

    def extract_unique_sequences(self):
        """从去重后的scFv氨基酸文件中提取唯一ID，然后从缓存字典中筛选"""
        unique_aa_path = Config.ANALYSE_DIR / "scfv_aa_unique.txt"
        if not unique_aa_path.exists():
            self.log("错误: 未找到去重后的序列文件", "ERROR")
            return False

        ordered_ids = []
        seen_ids = set()
        try:
            for seq_record in SeqIO.parse(unique_aa_path, "fasta"):
                if seq_record.id not in seen_ids:
                    seen_ids.add(seq_record.id)
                    ordered_ids.append(seq_record.id)
        except Exception as e:
            self.log(f"读取唯一序列ID时出错: {str(e)}", "ERROR")
            return False
        if not ordered_ids:
            self.log("错误: 去重后的scFv序列为空", "ERROR")
            return False
        self.log(f"找到 {len(ordered_ids)} 个唯一序列ID")

        # 任一缓存为空时重新加载，避免部分链缺失
        if not self.vh_nt_dict or not self.vl_nt_dict or not self.scfv_nt_dict:
            self._load_nt_dicts_from_files('scfv')

        nt_mappings = [
            ("VH_nt_scfvunique.txt", self.vh_nt_dict),
            ("VL_nt_scfvunique.txt", self.vl_nt_dict),
            ("scfv_nt_unique.txt", self.scfv_nt_dict),
        ]
        missing = []
        for seq_id in ordered_ids:
            for out_file, nt_dict in nt_mappings:
                if seq_id not in nt_dict:
                    missing.append(f"{out_file}:{seq_id}")
        if missing:
            shown = ", ".join(missing[:10])
            more = f" 等 {len(missing)} 项" if len(missing) > 10 else ""
            self.log(f"错误: 唯一ID在核酸缓存中缺失: {shown}{more}", "ERROR")
            return False

        aa_mappings = [
            ("VH_aa_seq.txt", "VH_aa_scfvunique.txt"),
            ("VL_aa_seq.txt", "VL_aa_scfvunique.txt"),
        ]
        for input_file, _output_file in aa_mappings:
            if not (Config.ANALYSE_DIR / input_file).exists():
                self.log(f"错误: 缺少氨基酸文件 {input_file}", "ERROR")
                return False

        staged = []
        try:
            for out_file, nt_dict in nt_mappings:
                out_path = Config.ANALYSE_DIR / out_file
                tmp_path = out_path.with_name(out_path.name + ".tmp")
                with open(tmp_path, 'w', encoding='utf-8', buffering=1024 * 1024) as f:
                    for seq_id in ordered_ids:
                        f.write(f'>{seq_id}\n{nt_dict[seq_id]}\n')
                staged.append((tmp_path, out_path, out_file, len(ordered_ids)))

            for input_file, output_file in aa_mappings:
                in_path = Config.ANALYSE_DIR / input_file
                out_path = Config.ANALYSE_DIR / output_file
                tmp_path = out_path.with_name(out_path.name + ".tmp")
                count = 0
                written_ids = set()
                with open(tmp_path, 'w', encoding='utf-8', buffering=1024 * 1024) as out_f:
                    for rec in SeqIO.parse(in_path, 'fasta'):
                        if rec.id in seen_ids and rec.id not in written_ids:
                            out_f.write(f'>{rec.id}\n{rec.seq}\n')
                            written_ids.add(rec.id)
                            count += 1
                if count != len(ordered_ids):
                    self.log(
                        f"错误: {output_file} 记录数 {count} 与唯一ID数 {len(ordered_ids)} 不一致",
                        "ERROR"
                    )
                    return False
                staged.append((tmp_path, out_path, output_file, count))

            for tmp_path, out_path, label, count in staged:
                os.replace(tmp_path, out_path)
                self.log(f"已生成: {label} ({count}条)")
            return True
        finally:
            for tmp_path, _out_path, _label, _count in staged:
                if tmp_path.exists():
                    try:
                        tmp_path.unlink()
                    except OSError:
                        pass

    def _generate_nt_from_aa(self, aa_file, nt_dict, out_file):
        """根据去重后的氨基酸ID，从核酸字典中生成去重核酸文件"""
        aa_path = Config.ANALYSE_DIR / aa_file
        if not aa_path.exists():
            self.log(f"警告: 未找到去重后的氨基酸文件 {aa_file}", "WARNING")
            return False

        ordered_ids = []
        seen_ids = set()
        for seq_record in SeqIO.parse(aa_path, "fasta"):
            if seq_record.id not in seen_ids:
                seen_ids.add(seq_record.id)
                ordered_ids.append(seq_record.id)
        if not ordered_ids:
            self.log(f"错误: {aa_file} 中没有序列", "ERROR")
            return False

        missing = [seq_id for seq_id in ordered_ids if seq_id not in nt_dict]
        if missing:
            shown = ", ".join(missing[:10])
            more = f" 等 {len(missing)} 个ID" if len(missing) > 10 else ""
            self.log(f"错误: 以下ID不在核酸字典中: {shown}{more}", "ERROR")
            return False

        out_path = Config.ANALYSE_DIR / out_file
        tmp_path = out_path.with_name(out_path.name + ".tmp")
        try:
            with open(tmp_path, 'w', encoding='utf-8', buffering=1024 * 1024) as f:
                for seq_id in ordered_ids:
                    f.write(f'>{seq_id}\n{nt_dict[seq_id]}\n')
            os.replace(tmp_path, out_path)
            self.log(f"已生成: {out_file} ({len(ordered_ids)}条)")
            return True
        finally:
            if tmp_path.exists():
                try:
                    tmp_path.unlink()
                except OSError:
                    pass

    def extract_fab_unique_sequences(self):
        """根据去重后的VH/VL氨基酸文件提取对应核苷酸序列"""
        if not self.vh_nt_dict or not self.vl_nt_dict:
            self._load_nt_dicts_from_files('fab')

        vh_ok = self._generate_nt_from_aa(
            "VH_aa_unique.txt", self.vh_nt_dict, "VH_nt_fabunique.txt")
        vl_ok = self._generate_nt_from_aa(
            "VL_aa_unique.txt", self.vl_nt_dict, "VL_nt_fabunique.txt")
        return vh_ok and vl_ok

    def extract_vhh_unique_nt(self):
        """根据去重后的VHH氨基酸文件提取对应核苷酸序列"""
        if not self.vhh_nt_dict:
            self._load_nt_dicts_from_files('vhh')
        return self._generate_nt_from_aa(
            "VHH_aa_unique.txt", self.vhh_nt_dict, "VHH_nt_unique.txt")

    def _load_nt_dicts_from_files(self, mode=None):
        """从文件加载当前模式需要的核酸字典"""
        self.vh_nt_dict.clear()
        self.vl_nt_dict.clear()
        self.scfv_nt_dict.clear()
        self.vhh_nt_dict.clear()
        if mode == "vhh":
            file_dicts = [("VHH_nt_seq.txt", self.vhh_nt_dict)]
        elif mode == "fab":
            file_dicts = [
                ("VH_nt_seq.txt", self.vh_nt_dict),
                ("VL_nt_seq.txt", self.vl_nt_dict),
            ]
        elif mode == "scfv":
            file_dicts = [
                ("VH_nt_seq.txt", self.vh_nt_dict),
                ("VL_nt_seq.txt", self.vl_nt_dict),
                ("scfv_nt_seq.txt", self.scfv_nt_dict),
            ]
        else:
            file_dicts = [
                ("VH_nt_seq.txt", self.vh_nt_dict),
                ("VL_nt_seq.txt", self.vl_nt_dict),
                ("scfv_nt_seq.txt", self.scfv_nt_dict),
                ("VHH_nt_seq.txt", self.vhh_nt_dict),
            ]
        for filename, target_dict in file_dicts:
            path = Config.ANALYSE_DIR / filename
            if not path.exists():
                continue
            for record in SeqIO.parse(path, "fasta"):
                seq_id = record.id
                if seq_id in target_dict:
                    self.log(
                        f"警告: {filename} 中存在重复ID {seq_id}，已保留首次出现的记录",
                        "WARNING"
                    )
                    continue
                target_dict[seq_id] = str(record.seq)

    def _get_duplicate_file_pairs(self, mode=None):
        """当前模式下需要处理的（重复明细文件, 去重文件）列表"""
        if self.is_vhh_mode(mode):
            return [
                ("duplicates.detail_VHH_aa.txt", "VHH_aa_unique.txt")
            ]
        if self.is_fab_mode(mode):
            return [
                ("duplicates.detail_VH_aa.txt", "VH_aa_unique.txt"),
                ("duplicates.detail_VL_aa.txt", "VL_aa_unique.txt")
            ]
        return [
            ("duplicates.detail_scfv_aa.txt", "scfv_aa_unique.txt"),
            ("duplicates.detail_VH_aa_scfvunique.txt", "VH_aa_pureunique.txt"),
            ("duplicates.detail_VL_aa_scfvunique.txt", "VL_aa_pureunique.txt")
        ]

    def process_duplicate_files(self):
        mode = self.analysis_mode.get()
        self.set_processing_state(True)
        thread = threading.Thread(target=self._process_duplicate_files_thread, args=(mode,))
        thread.daemon = True
        thread.start()

    def _process_duplicate_files_thread(self, mode=None):
        try:
            self.log("开始处理重复文件...")
            duplicate_files = self._get_duplicate_file_pairs(mode)
            processed = 0
            for i, (dup_file, unique_file) in enumerate(duplicate_files):
                dup_path = Config.ANALYSE_DIR / dup_file
                unique_path = Config.ANALYSE_DIR / unique_file
                if dup_path.exists() and unique_path.exists():
                    if self._process_single_duplicate_file(dup_path, unique_path):
                        processed += 1
                        self.log(f"已处理: {dup_file}")
                else:
                    self.log(f"跳过: {dup_file} 或 {unique_file} 不存在", "WARNING")
                progress = int((i + 1) / len(duplicate_files) * 100)
                self._post_ui(lambda p=progress: self.progress_var.set(p))
            self.log(f"重复文件处理完成: 成功处理 {processed} 个文件")
            self._post_ui(self.refresh_results)
        except Exception as e:
            self.log(f"处理重复文件时出错: {str(e)}", "ERROR")
            self.log(traceback.format_exc(), "ERROR")
        finally:
            self._post_ui(lambda: self.set_processing_state(False))

    def _process_single_duplicate_file(self, dup_path, unique_path):
        """将seqkit重复明细文件与唯一序列合并后生成Excel"""
        try:
            rows = []
            existing_ids = set()
            with open(dup_path, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    if '\t' in line:
                        group, id_text = line.split('\t', 1)
                    else:
                        group, id_text = "1", line
                    group = group.strip() or "1"
                    for seq_id in (part.strip() for part in id_text.split(',')):
                        if seq_id and seq_id not in existing_ids:
                            rows.append([group, seq_id])
                            existing_ids.add(seq_id)

            if unique_path.exists():
                for seq_record in SeqIO.parse(unique_path, 'fasta'):
                    seq_id = str(seq_record.id)
                    if seq_id not in existing_ids:
                        rows.append(["1", seq_id])
                        existing_ids.add(seq_id)

            df = pd.DataFrame(rows, columns=["Group", "Sequence_ID"])
            excel_path = dup_path.with_suffix('.xlsx')
            tmp_path = excel_path.with_name(excel_path.name + ".tmp")
            try:
                df.to_excel(tmp_path, index=False)
                os.replace(tmp_path, excel_path)
            finally:
                if tmp_path.exists():
                    try:
                        tmp_path.unlink()
                    except OSError:
                        pass
            self.log(f"已生成Excel文件: {excel_path.name} ({len(df)}行)")
            try:
                dup_path.unlink()
            except OSError:
                pass
            return True
        except Exception as e:
            self.log(f"处理文件 {dup_path.name} 时出错: {str(e)}", "ERROR")
            return False

    def _get_alignment_files(self, mode=None):
        """当前模式下需要进行多序列比对的氨基酸文件"""
        if self.is_vhh_mode(mode):
            return [
                "VHH_aa_unique.txt"
            ]
        if self.is_fab_mode(mode):
            return [
                "VH_aa_unique.txt",
                "VL_aa_unique.txt"
            ]
        return [
            "scfv_aa_unique.txt",
            "VH_aa_pureunique.txt",
            "VL_aa_pureunique.txt"
        ]

    def run_alignment(self):
        mode = self.analysis_mode.get()
        state = getattr(self, "workflow_state", None)
        if state is not None:
            state["aligned"] = False
        self._refresh_step_buttons()
        self.set_processing_state(True)
        thread = threading.Thread(target=self._run_alignment_thread, args=(mode,))
        thread.daemon = True
        thread.start()

    def _run_alignment_thread(self, mode=None):
        failures = []
        try:
            self.log("开始多序列比对...")
            files_to_align = self._get_alignment_files(mode)
            # 获取程序所在目录（clustalw2.exe 所在位置）
            base_dir = get_base_dir()

            for i, filename in enumerate(files_to_align):
                file_path = Config.ANALYSE_DIR / filename
                if file_path.exists():
                    self.log(f"正在比对: {filename}")
                    success, message = run_clustalw(str(file_path), base_dir)
                    self.log(message)
                    if success:
                        self.log(f"已完成比对: {filename}")
                    else:
                        if "序列数量少于2条" not in message:
                            failures.append(filename)
                        self.log(f"比对问题: {message}", "WARNING")
                else:
                    failures.append(filename)
                    self.log(f"警告: 比对文件 {filename} 不存在，跳过", "WARNING")

                progress = int((i + 1) / len(files_to_align) * 100)
                self._post_ui(lambda p=progress: self.progress_var.set(p))

            if failures:
                self.log(f"多序列比对结束，但有 {len(failures)} 个文件未完成", "WARNING")
            else:
                self.log("多序列比对完成")
            self._post_ui(self._mark_aligned)
        except Exception as e:
            self.log(f"比对时出错: {str(e)}", "ERROR")
            self.log(traceback.format_exc(), "ERROR")
        finally:
            self._post_ui(lambda: self.set_processing_state(False))

    def _get_chain_nt_unique_file(self, chain_type, mode=None):
        """当前模式下用于IgBlast的去重核酸文件名"""
        if self.is_vhh_mode(mode):
            return "VHH_nt_unique.txt"
        suffix = "fabunique" if self.is_fab_mode(mode) else "scfvunique"
        return f"{chain_type}_nt_{suffix}.txt"

    def _get_igblast_chains(self, mode=None):
        """当前模式下需要运行IgBlast的链列表"""
        return ['VHH'] if self.is_vhh_mode(mode) else ['VH', 'VL']

    def run_igblast(self):
        mode = self.analysis_mode.get()
        species = self.current_species.get()
        scheme = self.current_scheme.get()
        if not self._preflight_igblast(mode, species):
            return
        for chain in self._get_igblast_chains(mode):
            for suffix in ('-Seq_set.xlsx', '-Seq_unset.csv'):
                stale = Config.ANALYSE_DIR / f"{chain}{suffix}"
                if stale.exists():
                    try:
                        stale.unlink()
                    except OSError as e:
                        self.log(f"警告: 无法清理旧IgBlast文件 {stale.name}: {str(e)}", "WARNING")
        self.set_processing_state(True)
        thread = threading.Thread(
            target=self._run_igblast_thread, args=(species, scheme, mode))
        thread.daemon = True
        thread.start()

    def _run_igblast_thread(self, species=None, scheme=None, mode=None):
        failures = []
        try:
            species = species or self.current_species.get()
            scheme = scheme or self.current_scheme.get()
            mode = mode or self.analysis_mode.get()
            self.log(f"开始IgBlast分析 - 物种: {species}, 编号方案: {scheme}")
            chains = self._get_igblast_chains(mode)
            for i, chain_type in enumerate(chains):
                input_file = Config.ANALYSE_DIR / self._get_chain_nt_unique_file(chain_type, mode)
                if not input_file.exists():
                    failures.append(chain_type)
                    self.log(f"警告: 输入文件 {input_file} 不存在，跳过 {chain_type} 分析", "WARNING")
                    continue
                cmd = self.build_igblast_args(chain_type, species, scheme, mode)
                # 不显示具体命令，保持日志简洁
                # self.log(f"执行命令: {self.build_igblast_command(chain_type, species, scheme)}")
                try:
                    result = self.safe_subprocess_run(cmd, shell=False, timeout=900)
                    if result.returncode == 0:
                        if self.process_igblast_output(chain_type):
                            self.log(f"igblast分析完成: {chain_type}")
                        else:
                            failures.append(chain_type)
                            self.log(f"igblast输出处理失败: {chain_type}", "ERROR")
                    else:
                        failures.append(chain_type)
                        self.log(f"igblast分析失败: {chain_type}", "ERROR")
                        self.log(f"错误输出: {result.stderr}", "ERROR")
                except Exception as e:
                    failures.append(chain_type)
                    self.log(f"执行igblast命令时出错: {str(e)}", "ERROR")
                progress = int((i + 1) / len(chains) * 100)
                self._post_ui(lambda p=progress: self.progress_var.set(p))
            if failures:
                self.log(f"IgBlast分析结束，但有 {len(failures)} 个链未完成: {', '.join(failures)}", "WARNING")
            else:
                self.log("IgBlast分析完成")
            self._post_ui(self.refresh_results)
        except Exception as e:
            self.log(f"IgBlast分析时出错: {str(e)}", "ERROR")
        finally:
            self._post_ui(lambda: self.set_processing_state(False))

    def build_igblast_args(self, chain_type, species, numbering_scheme, mode=None):
        """构造 igblastn 参数列表，路径含空格时也不会被 shell 拆分"""
        if species not in Config.SPECIES_DB:
            raise ValueError(f"不支持的物种: {species}")
        db_info = Config.SPECIES_DB[species]

        if chain_type == 'VHH':
            organism = 'human' if species == 'human_VLL' else species
            cmd_parts = [
                str(Config.BIN_DIR / "igblastn"),
                "-germline_db_V", str(Config.DATABASE_DIR / db_info['vh']),
                "-germline_db_D", str(Config.DATABASE_DIR / db_info['d']),
                "-germline_db_J", str(Config.DATABASE_DIR / db_info['jh']),
                "-auxiliary_data", str(Config.OPTIONAL_FILE_DIR / db_info['aux']),
                "-show_translation",
                "-domain_system", numbering_scheme,
                "-outfmt", "19",
                "-extend_align3end",
                "-extend_align5end",
                "-num_alignments_V", "1",
                "-num_alignments_D", "1",
                "-num_alignments_J", "1",
                "-organism", organism,
                "-query", str(Config.ANALYSE_DIR / self._get_chain_nt_unique_file(chain_type, mode)),
                "-out", str(Config.ANALYSE_DIR / f"{chain_type}-Seq_unset.csv")
            ]
            return cmd_parts

        if chain_type == 'VH':
            v_db = db_info['vh']
            d_db = db_info['d']
            j_db = db_info['jh']
        else:
            if 'vl' not in db_info or 'jl' not in db_info:
                raise ValueError(f"物种 {species} 未配置VL数据库，无法分析VL")
            v_db = db_info['vl']
            d_db = db_info['d']
            j_db = db_info['jl']
        cmd_parts = [
            str(Config.BIN_DIR / "igblastn"),
            "-germline_db_V", str(Config.DATABASE_DIR / v_db),
            "-germline_db_D", str(Config.DATABASE_DIR / d_db),
            "-germline_db_J", str(Config.DATABASE_DIR / j_db),
            "-auxiliary_data", str(Config.OPTIONAL_FILE_DIR / db_info['aux']),
            "-show_translation",
            "-domain_system", numbering_scheme,
            "-outfmt", "19",
            "-organism", species,
            "-query", str(Config.ANALYSE_DIR / self._get_chain_nt_unique_file(chain_type, mode)),
            "-out", str(Config.ANALYSE_DIR / f"{chain_type}-Seq_unset.csv")
        ]
        return cmd_parts

    def build_igblast_command(self, chain_type, species, numbering_scheme, mode=None):
        """返回供日志/兼容使用的命令字符串；实际执行请使用 build_igblast_args"""
        return shlex.join(
            self.build_igblast_args(chain_type, species, numbering_scheme, mode))

    def process_igblast_output(self, chain_type):
        input_path = Config.ANALYSE_DIR / f"{chain_type}-Seq_unset.csv"
        output_path = Config.ANALYSE_DIR / f"{chain_type}-Seq_set.xlsx"
        if not input_path.exists():
            self.log(f"警告: igblast输出文件 {input_path} 不存在", "WARNING")
            return False
        if input_path.stat().st_size == 0:
            self.log(f"警告: {input_path.name} 为空，未生成Excel", "WARNING")
            return False

        tmp_path = output_path.with_name(output_path.name + ".tmp")
        try:
            with open(input_path, 'r', encoding='utf-8') as f:
                first_line = f.readline()
            if not first_line.strip():
                self.log(f"警告: {input_path.name} 为空，未生成Excel", "WARNING")
                return False
            sep = '\t' if '\t' in first_line else ','
            df = pd.read_csv(
                input_path, sep=sep, dtype=str, keep_default_na=False,
                encoding='utf-8'
            )
            if df.empty and len(df.columns) == 0:
                self.log(f"警告: {input_path.name} 没有可写入的数据", "WARNING")
                return False
            df.to_excel(tmp_path, index=False)
            os.replace(tmp_path, output_path)
            self.log(f"已生成Excel文件: {output_path.name} ({len(df)}行)")
            try:
                input_path.unlink()
            except OSError as e:
                self.log(f"警告: 无法删除临时文件 {input_path.name}: {str(e)}", "WARNING")
            return True
        except pd.errors.EmptyDataError:
            self.log(f"警告: {input_path.name} 为空，未生成Excel", "WARNING")
            return False
        except Exception as e:
            self.log(f"处理igblast输出文件时出错: {str(e)}", "ERROR")
            return False
        finally:
            if tmp_path.exists():
                try:
                    tmp_path.unlink()
                except OSError:
                    pass

    def run_all_steps(self):
        if not self.seq_files:
            messagebox.showwarning("警告", "请先选择序列文件")
            return
        mode = self.analysis_mode.get()
        separators = self._validate_separators()
        if separators is None:
            return
        species = self.current_species.get()
        scheme = self.current_scheme.get()
        skip_alignment = self.skip_alignment.get()
        skip_note = "\n（已勾选跳过比对步骤）" if skip_alignment else ""
        if messagebox.askyesno("确认",
                               "将依次执行以下步骤:\n1. 序列提取\n2. 序列翻译\n3. 序列去重\n4. 处理重复文件\n5. 多序列比对\n6. IgBlast分析" + skip_note + "\n\n是否继续?"):
            if not self._preflight_for_run_all(mode, species, scheme, skip_alignment):
                return
            self._archive_old_outputs()
            self._invalidate_workflow()
            self._write_run_manifest(mode, separators, species, scheme, skip_alignment)
            self._active_mode = mode
            self._active_separators = separators
            self.set_processing_state(True)
            thread = threading.Thread(
                target=self._run_all_steps_thread,
                args=(mode, separators, species, scheme, skip_alignment)
            )
            thread.daemon = True
            thread.start()

    def _run_all_steps_thread(self, mode=None, separators=None, species=None,
                              scheme=None, skip_alignment=None):
        manifest_status = "failed"
        had_failure = False
        try:
            mode = mode or getattr(self, "_active_mode", None) or self.analysis_mode.get()
            if separators is None:
                separators = getattr(self, "_active_separators", None)
            if species is None:
                species = self.current_species.get()
            if scheme is None:
                scheme = self.current_scheme.get()
            if skip_alignment is None:
                skip_alignment = self.skip_alignment.get()

            self.log("步骤1: 提取序列...")
            success, fail = self._extract_current_sequences(mode, separators)
            if success == 0:
                self.log("步骤1失败，终止后续步骤", "ERROR")
                return
            if fail:
                had_failure = True
                self.log(f"步骤1有 {fail} 条序列提取失败，详见 Error_seq.txt", "WARNING")
            self._post_ui(self._mark_extracted)
            self._post_ui(lambda: self.progress_var.set(15))

            self.log("步骤2: 翻译序列...")
            if not self._translate_current_sequences(mode):
                self.log("翻译失败，终止后续步骤", "ERROR")
                return
            self._post_ui(self._mark_translated)
            self._post_ui(lambda: self.progress_var.set(30))

            self.log("步骤3: 序列去重...")
            if not self._deduplicate_current_sequences(mode):
                self.log("去重失败，终止后续步骤", "ERROR")
                return
            self._post_ui(self._mark_deduplicated)
            self._post_ui(lambda: self.progress_var.set(50))

            self.log("步骤4: 处理重复文件...")
            dup_files = self._get_duplicate_file_pairs(mode)
            for dup, uniq in dup_files:
                dup_path = Config.ANALYSE_DIR / dup
                uniq_path = Config.ANALYSE_DIR / uniq
                if dup_path.exists() and uniq_path.exists():
                    if not self._process_single_duplicate_file(dup_path, uniq_path):
                        had_failure = True
                else:
                    had_failure = True
            self._post_ui(lambda: self.progress_var.set(65))

            if skip_alignment:
                self.log("已勾选跳过比对步骤，跳过步骤5")
            else:
                self.log("步骤5: 多序列比对...")
                # 获取程序所在目录（clustalw2.exe 所在位置）
                base_dir = get_base_dir()
                for fname in self._get_alignment_files(mode):
                    fpath = Config.ANALYSE_DIR / fname
                    if fpath.exists():
                        self.log(f"正在比对: {fname}")
                        success, message = run_clustalw(str(fpath), base_dir)
                        self.log(message)
                        if not success and "序列数量少于2条" not in message:
                            had_failure = True
                            self.log(f"注意: {message}", "WARNING")
                    else:
                        had_failure = True
                        self.log(f"警告: 比对文件 {fname} 不存在，跳过", "WARNING")
            self._post_ui(self._mark_aligned)
            self._post_ui(lambda: self.progress_var.set(80))

            self.log("步骤6: IgBlast分析...")
            for chain in self._get_igblast_chains(mode):
                infile = Config.ANALYSE_DIR / self._get_chain_nt_unique_file(chain, mode)
                if infile.exists():
                    cmd = self.build_igblast_args(chain, species, scheme, mode)
                    res = self.safe_subprocess_run(cmd, shell=False, timeout=900)
                    if res.returncode == 0:
                        if self.process_igblast_output(chain):
                            self.log(f"IgBlast {chain} 完成")
                        else:
                            had_failure = True
                            self.log(f"IgBlast输出处理失败: {chain}", "ERROR")
                    else:
                        had_failure = True
                        self.log(f"IgBlast {chain} 失败", "ERROR")
                else:
                    had_failure = True
                    self.log(f"警告: IgBlast输入文件不存在: {infile.name}", "WARNING")
            self._post_ui(lambda: self.progress_var.set(100))
            if had_failure:
                manifest_status = "completed_with_warnings"
                self.log("所有步骤执行结束，但存在失败项，请检查日志", "WARNING")
            else:
                manifest_status = "completed"
                self.log("所有步骤完成！")
            self._post_ui(self.refresh_results)
        except Exception as e:
            self.log(f"运行所有步骤时出错: {str(e)}", "ERROR")
            self.log(traceback.format_exc(), "ERROR")
        finally:
            self._update_run_manifest(manifest_status)
            self._post_ui(lambda: self.set_processing_state(False))

    def on_species_changed(self, event=None):
        self.log(f"物种已设置为: {self.current_species.get()}")

    def on_scheme_changed(self, event=None):
        self.log(f"编号方案已设置为: {self.current_scheme.get()}")

    @staticmethod
    def _database_exists(filename):
        """BLAST数据库以 .nhr/.nin/.nsq 等索引文件存在，按前缀判断"""
        base_path = Config.DATABASE_DIR / filename
        if base_path.is_file():
            return True
        return any(Config.DATABASE_DIR.glob(f"{filename}.*"))

    def _check_missing_databases(self, mode, species):
        missing = []
        db_info = Config.SPECIES_DB.get(species)
        if not db_info:
            return [f"未知物种数据库映射: {species}"]
        keys = ['vh', 'd', 'jh']
        if mode != 'vhh':
            keys.extend(['vl', 'jl'])
        for key in keys:
            filename = db_info.get(key)
            if not filename:
                missing.append(f"{species}.{key} 未配置")
                continue
            if not self._database_exists(filename):
                missing.append(str(Config.DATABASE_DIR / filename))
        aux = db_info.get('aux')
        if aux and not (Config.OPTIONAL_FILE_DIR / aux).is_file():
            missing.append(str(Config.OPTIONAL_FILE_DIR / aux))
        return missing

    def _preflight_igblast(self, mode, species):
        missing = []
        if not self._resolve_tool("igblastn"):
            missing.append("igblastn（bin 目录或 PATH）")
        missing.extend(self._check_missing_databases(mode, species))
        if not missing:
            self.log("预检通过: IgBlast 工具和数据库配置完整")
            return True
        message = "检测到以下 IgBlast 工具或数据库缺失:\n" + "\n".join(
            f"- {item}" for item in missing[:15])
        if len(missing) > 15:
            message += f"\n- 其余 {len(missing) - 15} 项..."
        self.log("IgBlast预检发现问题: " + "；".join(missing[:10]), "WARNING")
        return messagebox.askyesno("预检警告", message + "\n\n是否仍继续?")

    def _preflight_for_run_all(self, mode, species, scheme, skip_alignment):
        missing = []
        if not self._resolve_tool("seqkit"):
            missing.append("seqkit（bin 目录或 PATH）")
        if not skip_alignment and not resolve_clustalw():
            missing.append("clustalw2（程序目录、bin 目录或 PATH）")
        if not self._resolve_tool("igblastn"):
            missing.append("igblastn（bin 目录或 PATH）")
        missing.extend(self._check_missing_databases(mode, species))
        if not missing:
            self.log("预检通过: 外部工具和数据库配置完整")
            return True
        message = "检测到以下工具或数据库缺失:\n" + "\n".join(
            f"- {item}" for item in missing[:15])
        if len(missing) > 15:
            message += f"\n- 其余 {len(missing) - 15} 项..."
        self.log("预检发现问题: " + "；".join(missing[:10]), "WARNING")
        return messagebox.askyesno("预检警告", message + "\n\n是否仍继续?")

    def _write_run_manifest(self, mode, separators, species, scheme, skip_alignment):
        manifest = {
            "program_version": self.APP_VERSION,
            "started_at": datetime.now().isoformat(timespec='seconds'),
            "status": "running",
            "mode": mode,
            "species": species,
            "scheme": scheme,
            "skip_alignment": bool(skip_alignment),
            "separators": dict(separators),
            "input_files": [os.path.abspath(path) for path in self.seq_files],
        }
        path = Config.ANALYSE_DIR / 'run_manifest.json'
        tmp_path = path.with_name(path.name + '.tmp')
        try:
            Config.ANALYSE_DIR.mkdir(parents=True, exist_ok=True)
            with open(tmp_path, 'w', encoding='utf-8') as f:
                json.dump(manifest, f, ensure_ascii=False, indent=2)
            os.replace(tmp_path, path)
            self._manifest_path = path
        except Exception as e:
            self.log(f"警告: 无法写入运行清单: {str(e)}", "WARNING")
        finally:
            if tmp_path.exists():
                try:
                    tmp_path.unlink()
                except OSError:
                    pass

    def _update_run_manifest(self, status, message=""):
        path = getattr(self, "_manifest_path", None)
        if path is None or not path.exists():
            return
        try:
            data = json.loads(path.read_text(encoding='utf-8'))
            data["status"] = status
            data["finished_at"] = datetime.now().isoformat(timespec='seconds')
            if message:
                data["message"] = message
            tmp_path = path.with_name(path.name + '.tmp')
            with open(tmp_path, 'w', encoding='utf-8') as f:
                json.dump(data, f, ensure_ascii=False, indent=2)
            os.replace(tmp_path, path)
        except Exception as e:
            self.log(f"警告: 无法更新运行清单: {str(e)}", "WARNING")

    def check_tools(self):
        self.tool_status_text.delete(1.0, tk.END)

        # 检查 clustalw2
        clustalw_path = resolve_clustalw()
        if clustalw_path:
            self.tool_status_text.insert(tk.END, f"[正常] clustalw2: {clustalw_path}\n")
        else:
            self.tool_status_text.insert(tk.END, "[缺失] clustalw2: 请放在程序目录/bin，或加入 PATH\n")

        # 检查 seqkit（正确命令是 seqkit version）
        seqkit_path = self._resolve_tool("seqkit")
        if seqkit_path:
            try:
                result = self.safe_subprocess_run([seqkit_path, "version"], shell=False)
                output = ((result.stdout or "") + (result.stderr or "")).strip()
                if result.returncode == 0:
                    detail = f" - {output.splitlines()[0]}" if output else ""
                    self.tool_status_text.insert(tk.END, f"[正常] seqkit: {seqkit_path}{detail}\n")
                else:
                    self.tool_status_text.insert(
                        tk.END, f"[异常] seqkit: 返回码 {result.returncode} ({seqkit_path})\n")
                    if output:
                        self.tool_status_text.insert(tk.END, f"  {output[:300]}\n")
            except Exception as e:
                self.tool_status_text.insert(tk.END, f"[异常] seqkit: 检查失败 - {str(e)}\n")
        else:
            self.tool_status_text.insert(tk.END, "[缺失] seqkit: 未找到\n")

        # 检查 igblastn
        igblast_path = self._resolve_tool("igblastn")
        if igblast_path:
            try:
                result = self.safe_subprocess_run([igblast_path, "-version"], shell=False)
                if result.returncode == 0 or "version" in result.stderr.lower():
                    self.tool_status_text.insert(tk.END, f"[正常] igblastn: {igblast_path}\n")
                else:
                    self.tool_status_text.insert(tk.END, f"[异常] igblastn: 不可用 ({igblast_path})\n")
            except Exception as e:
                self.tool_status_text.insert(tk.END, f"[异常] igblastn: 检查失败 - {str(e)}\n")
        else:
            self.tool_status_text.insert(tk.END, f"[缺失] igblastn: 请检查 {Config.BIN_DIR}\n")

    def refresh_results(self):
        for item in self.results_tree.get_children():
            self.results_tree.delete(item)
        count = 0
        if Config.ANALYSE_DIR.exists():
            for file_path in sorted(Config.ANALYSE_DIR.glob("*"), key=lambda p: p.name.lower()):
                if file_path.is_file():
                    stat = file_path.stat()
                    size = self.format_file_size(stat.st_size)
                    mtime = datetime.fromtimestamp(stat.st_mtime).strftime("%Y-%m-%d %H:%M:%S")
                    tag = 'odd' if count % 2 else 'even'
                    self.results_tree.insert("", tk.END, values=(file_path.name, size, mtime), tags=(tag,))
                    count += 1
        if hasattr(self, "results_count_label"):
            self.results_count_label.config(text=f"共 {count} 个结果文件")

    def preview_selected_file(self):
        selection = self.results_tree.selection()
        if not selection:
            messagebox.showinfo("提示", "请先选择一个文件")
            return
        item = self.results_tree.item(selection[0])
        file_name = item['values'][0]
        file_path = Config.ANALYSE_DIR / file_name
        binary_suffixes = {'.xlsx', '.xls', '.png', '.jpg', '.jpeg', '.gif', '.pdf', '.zip'}
        if file_path.suffix.lower() in binary_suffixes:
            self.preview_text.delete(1.0, tk.END)
            self.preview_text.insert(1.0, "该文件为二进制格式，暂不支持文本预览。")
            return
        try:
            with open(file_path, 'r', encoding='utf-8', errors='replace') as f:
                content = f.read(5000)
            self.preview_text.delete(1.0, tk.END)
            self.preview_text.insert(1.0, content)
            if len(content) == 5000:
                self.preview_text.insert(tk.END, "\n\n...（仅显示前 5000 个字符）")
        except Exception as e:
            self.preview_text.delete(1.0, tk.END)
            self.preview_text.insert(1.0, f"无法读取文件: {str(e)}")

    def format_file_size(self, size):
        for unit in ['B', 'KB', 'MB', 'GB']:
            if size < 1024.0:
                return f"{size:.1f} {unit}"
            size /= 1024.0
        return f"{size:.1f} TB"

    def log(self, message, level="INFO"):
        self.log_handler.put(message, level)

    def add_log_message(self, timestamp, level, message):
        self.output_text.insert(tk.END, f"[{timestamp}] {message}\n")
        self.output_text.see(tk.END)
        self.log_text.insert(tk.END, f"[{timestamp}] ")
        if level == "ERROR":
            self.log_text.insert(tk.END, f"ERROR: {message}\n", "error")
        elif level == "WARNING":
            self.log_text.insert(tk.END, f"WARNING: {message}\n", "warning")
        else:
            self.log_text.insert(tk.END, f"{message}\n")
        self.log_text.see(tk.END)
        self.status_label.config(text=message[:80])
        self._trim_text_widget(self.output_text, self.MAX_OUTPUT_LINES)
        self._trim_text_widget(self.log_text, self.MAX_LOG_LINES)

    @staticmethod
    def _trim_text_widget(widget, max_lines):
        end_line = int(widget.index('end-1c').split('.')[0])
        if end_line > max_lines:
            widget.delete('1.0', f'{end_line - max_lines + 1}.0')

    def _post_ui(self, callback, *args):
        """把 UI 更新投递到主线程队列，避免工作线程直接操作 Tk"""
        ui_queue = getattr(self, "ui_queue", None)
        if ui_queue is None:
            callback(*args)
            return
        ui_queue.put((callback, args))

    def _process_ui_queue(self):
        processed = 0
        ui_queue = getattr(self, "ui_queue", None)
        if ui_queue is None:
            return
        while processed < 200:
            try:
                callback, args = ui_queue.get_nowait()
            except queue.Empty:
                break
            try:
                callback(*args)
            except Exception as e:
                try:
                    self.log(f"UI回调执行出错: {str(e)}", "ERROR")
                except Exception:
                    pass
            processed += 1
        delay = 10 if processed >= 200 else 50
        self.root.after(delay, self._process_ui_queue)

    def process_log_queue(self):
        processed = self.log_handler.process_queue(200)
        delay = 10 if processed >= 200 else 100
        self.root.after(delay, self.process_log_queue)

    def clear_log(self):
        self.log_handler.clear()
        self.log_text.delete(1.0, tk.END)
        self.output_text.delete(1.0, tk.END)
        self.log("日志已清空")

    def save_log(self):
        file_path = filedialog.asksaveasfilename(
            title="保存日志",
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
        )
        if file_path:
            try:
                with open(file_path, 'w', encoding='utf-8') as f:
                    f.write(self.log_text.get(1.0, tk.END))
                self.log(f"日志已保存到: {file_path}")
            except Exception as e:
                messagebox.showerror("错误", f"保存日志失败: {str(e)}")

    def _refresh_step_buttons(self):
        """按工作流状态恢复步骤按钮，避免无条件启用所有步骤"""
        if not hasattr(self, "step1_btn"):
            return
        if self.is_processing or not self.seq_files:
            for name in ("step1_btn", "step2_btn", "step3_btn",
                         "step4_btn", "step5_btn", "step6_btn", "run_all_btn"):
                button = getattr(self, name, None)
                if button is not None:
                    button.config(state=tk.DISABLED)
            return

        self.step1_btn.config(state=tk.NORMAL)
        self.run_all_btn.config(state=tk.NORMAL)
        self.step2_btn.config(
            state=tk.NORMAL if self.workflow_state["extracted"] else tk.DISABLED)
        self.step3_btn.config(
            state=tk.NORMAL if self.workflow_state["translated"] else tk.DISABLED)
        dedup_state = (
            tk.NORMAL if self.workflow_state["deduplicated"] else tk.DISABLED
        )
        for name in ("step4_btn", "step5_btn", "step6_btn"):
            getattr(self, name).config(state=dedup_state)

    def _invalidate_workflow(self):
        """文件列表或分析模式变化后，清空缓存并失效下游步骤"""
        state = getattr(self, "workflow_state", None)
        if state is not None:
            for key in state:
                state[key] = False
        for cache in (self.vh_nt_dict, self.vl_nt_dict,
                      self.scfv_nt_dict, self.vhh_nt_dict):
            cache.clear()
        self._active_mode = None
        self._active_separators = None
        self._refresh_step_buttons()

    def _mark_extracted(self):
        state = getattr(self, "workflow_state", None)
        if state is not None:
            state["extracted"] = True
            state["translated"] = False
            state["deduplicated"] = False
            state["aligned"] = False
        self._refresh_step_buttons()
        self.refresh_results()

    def _mark_translated(self):
        state = getattr(self, "workflow_state", None)
        if state is not None:
            state["translated"] = True
            state["deduplicated"] = False
            state["aligned"] = False
        self._refresh_step_buttons()
        self.refresh_results()

    def _mark_deduplicated(self):
        state = getattr(self, "workflow_state", None)
        if state is not None:
            state["deduplicated"] = True
            state["aligned"] = False
        self._refresh_step_buttons()
        self.refresh_results()

    def _mark_aligned(self):
        state = getattr(self, "workflow_state", None)
        if state is not None:
            state["aligned"] = True
        self._refresh_step_buttons()

    def _set_inputs_state(self, processing):
        state = tk.DISABLED if processing else tk.NORMAL
        for name in ("include_subfolders_check", "skip_alignment_check",
                     "reset_separators_btn"):
            widget = getattr(self, name, None)
            if widget is not None:
                try:
                    widget.config(state=state)
                except tk.TclError:
                    pass
        for name in ("species_combo", "scheme_combo"):
            combo = getattr(self, name, None)
            if combo is not None:
                try:
                    combo.config(state='disabled' if processing else 'readonly')
                except tk.TclError:
                    pass
        if hasattr(self, "separator_fields_frame"):
            for child in self.separator_fields_frame.winfo_children():
                if isinstance(child, ttk.Entry):
                    try:
                        child.config(state=state)
                    except tk.TclError:
                        pass
        for check in getattr(self, "input_format_checks", []):
            try:
                check.config(state=state)
            except tk.TclError:
                pass
        for radio in getattr(self, "mode_radios", []):
            try:
                radio.config(state=state)
            except tk.TclError:
                pass

    def set_processing_state(self, processing):
        self.is_processing = processing
        state = tk.DISABLED if processing else tk.NORMAL
        buttons = [
            self.select_files_btn,
            self.select_folder_btn,
            self.remove_files_btn,
            self.clear_files_btn,
            self.step1_btn,
            self.step2_btn,
            self.step3_btn,
            self.step4_btn,
            self.step5_btn,
            self.step6_btn,
            self.run_all_btn,
            self.check_tools_btn
        ]
        for btn in buttons:
            try:
                btn.config(state=state)
            except tk.TclError:
                pass
        self._set_inputs_state(processing)
        if not processing:
            self.progress_var.set(0)
        self._refresh_step_buttons()

    def on_closing(self):
        if self.is_processing:
            if not messagebox.askokcancel(
                    "任务运行中", "任务仍在运行，退出可能留下不完整文件。确定退出吗？"):
                return
        elif not messagebox.askokcancel("退出", "确定要退出程序吗？"):
            return
        self.root.destroy()


def main():
    root = tk.Tk()
    AntibodyGUI(root)
    root.mainloop()


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"程序运行出错: {str(e)}")
        traceback.print_exc()
