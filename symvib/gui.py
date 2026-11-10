"""振動モード解析の GUI (Tkinter + matplotlib)。

起動: python vibration_gui.py  (配布版は SymVib.exe)

画面構成:
    上段   入力ファイル / 解析の種類 (分子・結晶) / 許容差 / 解析・保存ボタン
    左     解析結果の要約とモード一覧 (クリックで右の 3D 表示が切り替わる)
           下にテキストレポート (CLI と同じ内容)
    右     選んだモードの 3D 表示 (原子・結合・変位矢印)
           左ドラッグ: 回転 / ホイール: ズーム / 右ドラッグ: 平行移動
           左下に結晶軸 (a, b, c; 分子は x, y, z) の向きを VESTA 風に表示
           縮退成分の切替、矢印の倍率、矢印・アニメーション・単位胞の表示切替

解析は別スレッドで実行し、画面が固まらないようにする。

3D 表示の設計:
    構造が変わったときだけ原子・結合の描画物を作り直し、モードの切替・アニメーションの
    各コマでは座標だけを更新する。アニメーションは Tk のタイマー (after) で回し、
    世代番号で古いタイマーを確実に止める (matplotlib の FuncAnimation は、止めた後も
    ウィンドウの再描画で再開することがあり、古いモードが動き続ける原因になっていた)。
"""
from __future__ import annotations

import math
import os
import queue
import sys
import threading
import traceback
import warnings
import tkinter as tk
import tkinter.font as tkfont
from tkinter import filedialog, messagebox, ttk
from tkinter.scrolledtext import ScrolledText

import numpy as np

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg      # noqa: E402
from mpl_toolkits.mplot3d.art3d import Line3DCollection              # noqa: E402

# 3D 表示のタイトルに日本語を出せるよう、日本語フォントを優先 (なければ既定のフォント)
matplotlib.rcParams["font.sans-serif"] = ["Yu Gothic", "Meiryo", "MS Gothic", "BIZ UDGothic",
                                          "Hiragino Sans", "Noto Sans CJK JP", "IPAGothic",
                                          "IPAexGothic", "DejaVu Sans"]
from matplotlib.figure import Figure                                 # noqa: E402

from . import __version__                                            # noqa: E402
from .common import (AnalysisError, ACTIVITIES, bond_pairs,          # noqa: E402
                     covalent_radius, display_z, element_symbol)
from .i18n import LANGUAGES, get_language, set_language, t            # noqa: E402

_SETTINGS = os.path.join(os.path.expanduser("~"), ".symvib.json")
AUTHOR = "Yosuke Ishii"
CONTACT = "ishii.yosuke@nitech.ac.jp"


def _load_settings():
    """前回の設定 (表示言語など)。読めなければ空。"""
    import json
    try:
        with open(_SETTINGS, encoding="utf-8") as f:
            data = json.load(f)
    except (OSError, ValueError):
        return {}
    return data if isinstance(data, dict) else {}         # 壊れた設定ファイルは無視


_ERROR_LOG = os.path.join(os.path.expanduser("~"), "symvib_error.log")


def _write_error_log(text):
    """予期しないエラーの詳細をログファイルに追記し、そのパスを返す
    (配布用 exe はコンソールがないので、画面には要点だけを出す)。"""
    import datetime
    try:
        with open(_ERROR_LOG, "a", encoding="utf-8") as f:
            f.write(f"--- {datetime.datetime.now():%Y-%m-%d %H:%M:%S} "
                    f"(symvib {__version__})\n{text}\n")
    except OSError:
        pass
    return _ERROR_LOG


def add_bond_collection(ax, segments):
    """結合の線をまとめて 3D 軸に追加して返す。結合が 1 本もないときは None
    (空の Line3DCollection を add_collection3d すると matplotlib が例外を出す;
    VESTA の結合表にない組だけの結晶 (例: AsI3) で起きた)。"""
    if len(segments) == 0:
        return None
    col = Line3DCollection(segments, colors="#9a9a9a", linewidths=2)
    ax.add_collection3d(col)
    return col


def _save_settings(**kw):
    import json
    data = _load_settings()
    data.update(kw)
    try:
        with open(_SETTINGS, "w", encoding="utf-8") as f:
            json.dump(data, f)
    except OSError:
        pass

# 原子の表示色 (CPK 配色に近いもの)。ないものは灰色
_COLORS = {"H": "#f2f2f2", "C": "#505050", "N": "#3050f8", "O": "#ff0d0d", "F": "#90e050",
           "Cl": "#1ff01f", "Br": "#a62929", "I": "#940094", "S": "#e6c619", "P": "#ff8000",
           "B": "#ffb5b5", "Si": "#2f5fe0", "Na": "#ab5cf2", "K": "#8f40d4", "Ca": "#3dff00",
           "Mg": "#8aff00", "Al": "#bfa6a6", "Ti": "#bfc2c7", "Fe": "#e06633", "Co": "#f090a0",
           "Ni": "#50d050", "Cu": "#c88033", "Zn": "#7d80b0", "Sr": "#00ff00", "Ba": "#00c900",
           "Xe": "#429eb0", "Pt": "#d0d0e0", "U": "#008fff", "Li": "#cc80ff", "Ga": "#c28f8f",
           "V": "#a6a6ab", "Cr": "#8a99c7", "Mn": "#9c7ac7", "Zr": "#94e0e0", "Nb": "#73c2c9",
           "Mo": "#54b5b5", "W": "#2194d6", "Pb": "#575961", "Bi": "#9e4fb5", "La": "#70d4ff"}
_ACTIVITY_COLORS = {"Raman": "#c0392b", "IR": "#1f5fbf", "Raman+IR": "#8e44ad",
                    "Silent": "#7f7f7f"}
_ARROW_COLOR = "#e67e22"
_FRAMES = 24              # アニメーション 1 周期のコマ数
_FRAME_MS = 50            # コマの間隔 [ms]


# 日本語の字形 (日本の漢字の形) を持つフォント。前の方ほど優先する。
# 等幅 (レポートの桁がそろう) のものを先に並べる。Windows・macOS・Linux の順
JAPANESE_FONTS = ("BIZ UDゴシック", "BIZ UDGothic", "ＭＳ ゴシック", "MS Gothic",
                  "Osaka−等幅", "Osaka-Mono", "Noto Sans Mono CJK JP", "IPAGothic",
                  "IPAゴシック", "VL Gothic",
                  "Yu Gothic UI", "Yu Gothic", "游ゴシック", "Meiryo", "メイリオ",
                  "BIZ UDPゴシック", "BIZ UDPGothic", "Hiragino Sans", "ヒラギノ角ゴシック",
                  "Noto Sans CJK JP", "IPAexGothic", "IPAexゴシック")
DEFAULT_FONT_SIZE = 11
FONT_SIZES = tuple(range(8, 25))


def _font_families(root):
    """選べるフォント名 (日本語のフォントを先頭に。Windows の縦書き用 @ フォントは除く)。"""
    families = {f for f in tkfont.families(root) if f and not f.startswith("@")}
    japanese = [f for f in JAPANESE_FONTS if f in families]
    return japanese + sorted(families - set(japanese), key=str.lower)


def _report_font(root, settings=None):
    """モード一覧・レポートのフォント (名前, 大きさ)。保存した設定があればそれを使い、なければ
    日本語の字形を持つ等幅フォントを優先する (無いと中国語の字形のフォントで代用されることがある)。"""
    settings = settings or {}
    families = set(tkfont.families(root))
    size = settings.get("font_size")
    size = size if isinstance(size, int) and size in FONT_SIZES else DEFAULT_FONT_SIZE
    family = settings.get("font_family")
    if isinstance(family, str) and family in families:
        return (family, size)
    for name in JAPANESE_FONTS + ("Consolas",):
        if name in families:
            return (name, size)
    return ("TkFixedFont", size)


def _best_view(axes):
    """既定の視点 (elev, azim): 3 本の軸 (a, b, c / x, y, z) がどれも見やすい向き。

    各軸を画面に投影した長さの最小値 × 投影どうしの開き角 (|sin|) の最小値 を
    最大にする向きを、仰角 10–40°・方位角 5° 刻みで探す (結晶系によらず同じ規則)。
    軸の和の側 (軸が手前を向く側) から見る向きに限り、同点ならより正面に近い向き。"""
    axes = np.asarray(axes, float)
    axes = axes / np.linalg.norm(axes, axis=1, keepdims=True)
    front = axes.sum(axis=0)
    front /= max(np.linalg.norm(front), 1e-9)
    best, best_score = (20.0, -60.0), (-1.0, -1.0)
    for elev in range(10, 41, 5):
        e = np.radians(elev)
        for azim in range(-180, 180, 5):
            a = np.radians(azim)
            eye = np.array([np.cos(e) * np.cos(a), np.cos(e) * np.sin(a), np.sin(e)])
            if eye @ front <= 0:
                continue
            proj = axes - np.outer(axes @ eye, eye)               # 画面上の軸
            lens = np.linalg.norm(proj, axis=1)
            unit = proj / np.maximum(lens[:, None], 1e-9)
            spread = min(np.linalg.norm(np.cross(unit[i], unit[j]))   # |sin(開き角)|
                         for i, j in ((0, 1), (0, 2), (1, 2)))
            score = (round(lens.min() * spread, 3), round(float(eye @ front), 3))
            if score > best_score:
                best, best_score = (float(elev), float(azim)), score
    return best


def _atom_color(z):
    dz = display_z(z)
    return _COLORS.get(element_symbol(dz), "#a0a0a0") if dz > 0 else "#a0a0a0"


class AutoScrollbar(ttk.Scrollbar):
    """中身が収まるときは隠れ、はみ出したときだけ現れるスクロールバー (grid で配置する)。"""

    def set(self, first, last):
        if float(first) <= 0.0 and float(last) >= 1.0:
            self.grid_remove()
        else:
            self.grid()
        super().set(first, last)


class VibrationApp(ttk.Frame):
    """メインウィンドウ。"""

    def __init__(self, master):
        super().__init__(master, padding=6)
        self.master = master
        self.result = None          # 直近の AnalysisResult
        self.source_name = ""       # 入力ファイル名 (保存先の既定名に使う)
        self.scene = None           # 表示用の原子 (座標・元の原子番号・結合・色)
        self._artists = None        # 原子・結合・矢印の描画物
        self._anim_job = None       # アニメーションのタイマー ID
        self._anim_gen = 0          # アニメーションの世代 (古いタイマーを無効にする)
        self._anim_frame = 0
        self._jobs = queue.Queue()  # 解析スレッド → GUI への結果受け渡し
        self._i18n = []             # 言語切替時に文字列を付け直す関数
        self._warnings = []         # 読み込み時の警告 (レポートの末尾に表示)
        self._job_id = 0            # 解析の通し番号 (古い解析の結果を捨てる)
        self._running = False       # 解析中か
        self._lock_widgets = []     # 解析中は押せなくするボタン類
        self._build_widgets()
        self._set_kind("molecule")
        self.pack(fill="both", expand=True)

    def _update_row_height(self):
        """一覧の行の高さを、フォントの実際の高さに余白を足して決める。"""
        line = self.text_font.metrics("linespace")
        style = ttk.Style(self)
        style.configure("Modes.Treeview", font=self.text_font,
                        rowheight=line + max(6, line // 3))
        style.configure("Modes.Treeview.Heading", font=self.head_font)

    def apply_font(self, family, size, save=True):
        """モード一覧とレポートのフォントを変える (選んだものは次回の起動時にも使う)。"""
        for f in (self.text_font, self.head_font):
            f.configure(family=family, size=size)
        self._update_row_height()
        self._fit_columns()
        if save:
            _save_settings(font_family=family, font_size=size)

    def choose_font(self):
        """フォントの種類と大きさを選ぶ小さな画面。"""
        win = tk.Toplevel(self)
        win.title(t("gui.font_title"))
        win.transient(self.winfo_toplevel())
        win.resizable(False, False)
        families = _font_families(self)
        fam_var = tk.StringVar(value=self.text_font.actual("family"))
        size_var = tk.IntVar(value=abs(int(self.text_font.actual("size"))) or DEFAULT_FONT_SIZE)
        preview_font = tkfont.Font(win, font=self.text_font)
        body = ttk.Frame(win, padding=10)
        body.pack(fill="both", expand=True)
        ttk.Label(body, text=t("gui.font_family")).grid(row=0, column=0, sticky="w")
        fam_box = ttk.Combobox(body, textvariable=fam_var, values=families, state="readonly",
                               width=34, height=20)
        fam_box.grid(row=0, column=1, sticky="ew", padx=4, pady=2)
        ttk.Label(body, text=t("gui.font_size")).grid(row=1, column=0, sticky="w")
        size_box = ttk.Combobox(body, textvariable=size_var, values=FONT_SIZES, state="readonly",
                                width=5)
        size_box.grid(row=1, column=1, sticky="w", padx=4, pady=2)
        sample = tk.Label(body, text=t("gui.font_sample"), font=preview_font, anchor="w",
                          justify="left", relief="groove", padx=6, pady=6)
        sample.grid(row=2, column=0, columnspan=2, sticky="ew", pady=(8, 4))
        ttk.Label(body, text=t("gui.font_note"), foreground="#555555", wraplength=420,
                  justify="left").grid(row=3, column=0, columnspan=2, sticky="w")

        def preview(event=None):
            preview_font.configure(family=fam_var.get(), size=int(size_var.get()))
        fam_box.bind("<<ComboboxSelected>>", preview)
        size_box.bind("<<ComboboxSelected>>", preview)

        def ok():
            self.apply_font(fam_var.get(), int(size_var.get()))
            win.destroy()

        def default():
            family, size = _report_font(self)
            fam_var.set(family)
            size_var.set(size)
            preview()
        buttons = ttk.Frame(body)
        buttons.grid(row=4, column=0, columnspan=2, sticky="e", pady=(8, 0))
        ttk.Button(buttons, text=t("gui.font_default"), command=default).pack(side="left", padx=4)
        ttk.Button(buttons, text=t("gui.ok"), command=ok).pack(side="left", padx=4)
        ttk.Button(buttons, text=t("gui.cancel"), command=win.destroy).pack(side="left")
        win.bind("<Escape>", lambda e: win.destroy())
        win.grab_set()
        return win

    def _fit_columns(self):
        """モード一覧の各列の幅を、見出しと中身の文字列の幅 (現在のフォント) に合わせる。
        合計が一覧の幅を超えるときは横のスクロールバーが現れる。"""
        pad = 2 * self.text_font.measure("0")
        rows = [self.tree.item(i, "values") for i in self.tree.get_children()]
        for k, c in enumerate(self._cols):
            texts = [str(self.tree.heading(c, "text"))]
            width = self.head_font.measure(texts[0])
            for values in rows:
                if k < len(values):
                    width = max(width, self.text_font.measure(str(values[k])))
            self.tree.column(c, width=width + pad, minwidth=pad)

    # ── 表示言語 ───────────────────────────────────────────────────────
    def _tr(self, widget, key, option="text"):
        """widget の文字列を言語キーで設定し、言語切替時に付け直せるよう登録する。"""
        def apply():
            widget.configure(**{option: t(key)})
        apply()
        self._i18n.append(apply)
        return widget

    def _on_language(self, event=None):
        lang = {v: k for k, v in LANGUAGES.items()}[self.lang_var.get()]
        set_language(lang)
        _save_settings(language=lang)
        self._apply_language()

    def _apply_language(self):
        """画面の文字列・一覧・レポートを現在の言語で付け直す (解析結果はそのまま)。"""
        for apply in self._i18n:
            apply()
        self.lang_var.set(LANGUAGES[get_language()])
        self.master.title(t("gui.title"))
        self._set_kind(self.kind_var.get())
        if self.result is not None:
            self._fill_result_texts(keep_selection=True)
            if self.scene is not None:
                self.update_mode()                   # 3D 表示のタイトルも付け直す
        else:
            self.summary_var.set(t("gui.start"))
        if self._running:
            self.status_var.set(t("gui.running"))

    # ── 画面の組み立て ─────────────────────────────────────────────────
    def _build_widgets(self):
        top = ttk.Frame(self)
        top.pack(fill="x")
        self._tr(ttk.Label(top), "gui.input").pack(side="left")
        self.path_var = tk.StringVar()
        path_entry = ttk.Entry(top, textvariable=self.path_var, width=48)
        path_entry.pack(side="left", padx=4)
        open_btn = self._tr(ttk.Button(top, command=self.open_file), "gui.open")
        open_btn.pack(side="left")
        demo = self._tr(ttk.Menubutton(top), "gui.demo")
        self._lock_widgets += [path_entry, open_btn, demo]
        menu = tk.Menu(demo, tearoff=False)
        menu.add_command(command=lambda: self.run_demo("molecule"))
        menu.add_command(command=lambda: self.run_demo("crystal"))
        for idx, key in enumerate(("gui.demo_water", "gui.demo_graphite")):
            def relabel(idx=idx, key=key):
                menu.entryconfigure(idx, label=t(key))
            relabel()
            self._i18n.append(relabel)
        demo["menu"] = menu
        demo.pack(side="left", padx=4)

        ttk.Separator(top, orient="vertical").pack(side="left", fill="y", padx=6)
        self.kind_var = tk.StringVar(value="molecule")
        for value, key in (("molecule", "gui.molecule"), ("crystal", "gui.crystal")):
            rb = self._tr(ttk.Radiobutton(top, value=value, variable=self.kind_var,
                                          command=lambda: self._set_kind(self.kind_var.get())),
                          key)
            rb.pack(side="left")
            self._lock_widgets.append(rb)
        self.tol_label = ttk.Label(top)
        self.tol_label.pack(side="left", padx=(8, 2))
        self.tol_var = tk.StringVar()
        tol_entry = ttk.Entry(top, textvariable=self.tol_var, width=8)
        tol_entry.pack(side="left")
        self._lock_widgets.append(tol_entry)

        self.run_btn = self._tr(ttk.Button(top, command=self.run_analysis), "gui.run")
        self.run_btn.pack(side="left", padx=(10, 4))
        self.save_btn = self._tr(ttk.Button(top, command=self.save_modes, state="disabled"),
                                 "gui.save_modes")
        self.save_btn.pack(side="left")
        self.report_btn = self._tr(ttk.Button(top, command=self.save_report, state="disabled"),
                                   "gui.save_report")
        self.report_btn.pack(side="left", padx=4)

        # 右端: 言語の選択と進捗表示
        self.lang_var = tk.StringVar(value=LANGUAGES[get_language()])
        lang_box = ttk.Combobox(top, textvariable=self.lang_var, values=list(LANGUAGES.values()),
                                width=9, state="readonly")
        lang_box.bind("<<ComboboxSelected>>", self._on_language)
        lang_box.pack(side="right")
        self._tr(ttk.Label(top), "gui.language").pack(side="right", padx=(10, 2))
        self._tr(ttk.Button(top, command=self.choose_font), "gui.font").pack(side="right", padx=(10, 0))
        self.progress = ttk.Progressbar(top, mode="determinate", length=90)
        self.progress.pack(side="right")

        panes = ttk.PanedWindow(self, orient="horizontal")
        panes.pack(fill="both", expand=True, pady=(6, 0))

        # 左: 要約 + モード一覧 + レポート
        left = ttk.Frame(panes)
        panes.add(left, weight=1)
        self.summary_var = tk.StringVar(value=t("gui.start"))
        ttk.Label(left, textvariable=self.summary_var, wraplength=520,
                  justify="left").pack(fill="x", pady=(0, 4))
        # モード一覧とレポートは同じフォント・同じ大きさで表示する。行の高さは
        # フォントの実際の高さ (画面の拡大率を反映したピクセル数) から決める
        # (固定のピクセル値だと、高 DPI の画面で文字の上下が行からはみ出して切れる)
        font = _report_font(self, _load_settings())
        self.text_font = tkfont.Font(self, font=font)
        self.head_font = tkfont.Font(self, font=font)
        self.head_font.configure(weight="bold")
        self._update_row_height()
        self._cols = ("no", "label", "deg", "activity", "direction", "tensor")
        # 一覧 (+凡例) とレポートの境目は上下にドラッグして動かせる
        split = ttk.PanedWindow(left, orient="vertical")
        split.pack(fill="both", expand=True)
        upper = ttk.Frame(split)
        split.add(upper, weight=1)
        tree_frame = ttk.Frame(upper)
        tree_frame.pack(fill="both", expand=True)
        tree_frame.rowconfigure(0, weight=1)
        tree_frame.columnconfigure(0, weight=1)
        self.tree = ttk.Treeview(tree_frame, columns=self._cols, show="headings", height=8,
                                 selectmode="browse", style="Modes.Treeview")
        for c in self._cols:
            self.tree.column(c, anchor="w" if c == "tensor" else "center", stretch=False)

        def relabel_heads():
            for c, h in zip(self._cols, t("gui.cols").split("|")):
                self.tree.heading(c, text=h)
            self._fit_columns()
        relabel_heads()
        self._i18n.append(relabel_heads)
        for act, color in _ACTIVITY_COLORS.items():
            self.tree.tag_configure(act, foreground=color)
        self.tree.tag_configure("acoustic", foreground="#999999")
        # 一覧が縦・横に収まらないときだけスクロールバーを出す
        ysb = AutoScrollbar(tree_frame, orient="vertical", command=self.tree.yview)
        xsb = AutoScrollbar(tree_frame, orient="horizontal", command=self.tree.xview)
        self.tree.configure(yscrollcommand=ysb.set, xscrollcommand=xsb.set)
        self.tree.grid(row=0, column=0, sticky="nsew")
        ysb.grid(row=0, column=1, sticky="ns")
        xsb.grid(row=1, column=0, sticky="ew")
        self.tree.bind("<<TreeviewSelect>>", lambda e: self.show_selected_mode())

        legend = ttk.Frame(upper)
        legend.pack(side="bottom", fill="x", pady=2, before=tree_frame)   # 狭いときは一覧の側が縮む
        for act in ACTIVITIES:
            tk.Label(legend, text=f"■ {act}", fg=_ACTIVITY_COLORS[act]).pack(side="left", padx=4)
        self._tr(ttk.Label(legend, foreground="#555555"), "gui.legend").pack(side="left")

        self.report = ScrolledText(split, height=8, wrap="none", font=self.text_font)
        split.add(self.report, weight=1)
        self.report.configure(state="disabled")

        # 右: 3D 表示
        right = ttk.Frame(panes)
        panes.add(right, weight=1)
        self.fig = Figure(figsize=(5.5, 5.0), dpi=100)
        self.fig.subplots_adjust(left=0, right=1, bottom=0, top=0.95)
        self.ax = self.fig.add_subplot(111, projection="3d")
        # 左下の小さな軸表示 (結晶軸 a, b, c の向き)。本体の視点に合わせて回す
        self.axes_ax = self.fig.add_axes([0.0, 0.0, 0.22, 0.22], projection="3d")
        self.axes_ax.mouse_init(rotate_btn=[], pan_btn=[], zoom_btn=[])
        self.canvas = FigureCanvasTkAgg(self.fig, master=right)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        # 左ドラッグ: 回転, 右ドラッグ: 平行移動, ホイール: ズーム
        self.ax.mouse_init(rotate_btn=1, pan_btn=3, zoom_btn=[])
        self.canvas.mpl_connect("scroll_event", self._on_scroll)
        self.canvas.mpl_connect("draw_event", self._sync_axes_indicator)

        hint = ttk.Frame(right)
        hint.pack(fill="x")
        self._tr(ttk.Label(hint, foreground="#555555"), "gui.hint").pack(side="left")
        self._tr(ttk.Button(hint, command=self.save_image), "gui.save_image").pack(side="right")
        self._tr(ttk.Button(hint, command=self.reset_view),
                 "gui.reset_view").pack(side="right", padx=4)

        ctrl = ttk.Frame(right)
        ctrl.pack(fill="x", pady=4)
        self._tr(ttk.Label(ctrl), "gui.component").pack(side="left")
        self.comp_var = tk.StringVar()
        self.comp_box = ttk.Combobox(ctrl, textvariable=self.comp_var, width=4,
                                     state="readonly")
        self.comp_box.pack(side="left", padx=(2, 10))
        self.comp_box.bind("<<ComboboxSelected>>", lambda e: self.update_mode())
        self._tr(ttk.Label(ctrl), "gui.arrow_scale").pack(side="left")
        self.scale_var = tk.DoubleVar(value=1.0)
        ttk.Scale(ctrl, from_=0.2, to=3.0, variable=self.scale_var, length=110,
                  command=lambda v: self.update_mode()).pack(side="left", padx=(2, 10))
        opts = ttk.Frame(right)                   # 表示の切り替え (2 段目)
        opts.pack(fill="x", pady=(0, 4))
        self.arrow_var = tk.BooleanVar(value=True)
        self._tr(ttk.Checkbutton(opts, variable=self.arrow_var, command=self.update_mode),
                 "gui.arrows").pack(side="left")
        self.anim_var = tk.BooleanVar(value=False)
        self._tr(ttk.Checkbutton(opts, variable=self.anim_var, command=self.update_mode),
                 "gui.animation").pack(side="left", padx=10)
        self.cell_var = tk.BooleanVar(value=True)
        self.cell_check = self._tr(ttk.Checkbutton(opts, variable=self.cell_var,
                                                   command=self._toggle_cell), "gui.cell")
        self.cell_check.pack(side="left")

        # 最下部のバー: 左に状態表示、右にバージョン・作者・連絡先 (常に表示)
        bar = ttk.Frame(self, relief="sunken")
        bar.pack(fill="x", pady=(4, 0))
        self.status_var = tk.StringVar(value="")
        ttk.Label(bar, textvariable=self.status_var, anchor="w").pack(side="left", fill="x",
                                                                    expand=True, padx=(2, 0))
        ttk.Label(bar, text=f"symvib {__version__}   |   {AUTHOR} ({CONTACT})",
                  foreground="#555555").pack(side="right", padx=(8, 4))
        self._draw_empty()

    def _set_kind(self, kind):
        """解析の種類を切り替え、許容差の意味と既定値を合わせる。"""
        self.kind_var.set(kind)
        if kind == "molecule":
            self.tol_label.configure(text=t("gui.tol_mol"))
            if not self.tol_var.get() or self.tol_var.get() == "0.001":
                self.tol_var.set("0.1")
        else:
            self.tol_label.configure(text=t("gui.tol_cry"))
            if not self.tol_var.get() or self.tol_var.get() == "0.1":
                self.tol_var.set("0.001")

    # ── 解析の実行 ─────────────────────────────────────────────────────
    def open_file(self):
        path = filedialog.askopenfilename(
            title=t("gui.open_title"),
            filetypes=[(t("gui.ft_struct"), "*.xyz *.cif"),
                       (t("gui.ft_xyz"), "*.xyz"), ("CIF", "*.cif"), (t("gui.ft_all"), "*.*")])
        if not path:
            return
        self.path_var.set(path)
        ext = os.path.splitext(path)[1].lower()
        self._set_kind("crystal" if ext == ".cif" else "molecule")
        self.run_analysis()

    def run_demo(self, kind):
        self.path_var.set("")
        self._set_kind(kind)
        self.run_analysis()

    def run_analysis(self):
        path = self.path_var.get().strip()
        kind = self.kind_var.get()
        if self._running:
            return
        try:
            tol = float(self.tol_var.get())
            if not (math.isfinite(tol) and tol > 0):
                raise ValueError
        except ValueError:
            messagebox.showerror(t("gui.title"), t("gui.err_tol"))
            return
        if path and not os.path.exists(path):
            messagebox.showerror(t("gui.title"), t("gui.err_nofile", path=path))
            return
        name = (os.path.splitext(os.path.basename(path))[0] if path
                else ("water" if kind == "molecule" else "graphite"))
        self._job_id += 1
        self._busy(True, t("gui.running"))
        threading.Thread(target=self._worker, args=(self._job_id, path, kind, tol, name),
                         daemon=True).start()
        self.after(100, self._poll)

    def _worker(self, job, path, kind, tol, name):
        """別スレッドで解析し、(番号, 結果, エラー, 警告, 名前) をキューに入れる。"""
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                notes = []
                if kind == "molecule":
                    from .molecule import analyze_molecule, load_molecule, water_demo
                    cart, Z, sym = load_molecule(path, notes) if path else water_demo()
                    result = analyze_molecule(cart, Z, sym, tol=tol, notes=notes)
                else:
                    from .crystal import analyze_crystal, graphite_2h, load_crystal
                    lat, pos, Z, sym = load_crystal(path, notes) if path else graphite_2h()
                    result = analyze_crystal(lat, pos, Z, sym, symprec=tol, notes=notes)
            msgs = sorted({str(w.message) for w in caught
                           if not issubclass(w.category, DeprecationWarning)})
            self._jobs.put((job, result, None, msgs, name))
        except (AnalysisError, OSError) as e:
            self._jobs.put((job, None, str(e), [], name))
        except Exception as e:                               # 想定外: 要点と詳細 (ログ)
            log = _write_error_log(traceback.format_exc())
            self._jobs.put((job, None, t("gui.unexpected", err=f"{type(e).__name__}: {e}",
                                         log=log), [], name))

    def _poll(self):
        try:
            job, result, error, msgs, name = self._jobs.get_nowait()
        except queue.Empty:
            self.after(100, self._poll)
            return
        if job != self._job_id:                              # 古い解析の結果は捨てる
            self.after(100, self._poll)
            return
        self._busy(False)
        if error:
            self.status_var.set(t("gui.failed"))
            messagebox.showerror(t("gui.title"), t("gui.failed_dialog", err=error))
            return
        self.source_name = name              # 成功したときだけ保存先の既定名を変える
        self.show_result(result, msgs)

    def _busy(self, busy, message=""):
        """解析中は、解析・入力・保存の操作をできなくする。"""
        self._running = busy
        state = "disabled" if busy else "normal"
        for w in [self.run_btn] + self._lock_widgets:
            w.configure(state=state)
        has_result = self.result is not None
        self.save_btn.configure(state="disabled" if busy or not (has_result and self.result.modes)
                                else "normal")
        self.report_btn.configure(state="disabled" if busy or not has_result else "normal")
        if busy:
            self.progress.configure(mode="indeterminate")
            self.progress.start(12)
            self.status_var.set(message)
        else:
            self.progress.stop()
            self.progress.configure(mode="determinate", value=0)

    # ── 結果の表示 ─────────────────────────────────────────────────────
    def show_result(self, result, warnings_text=()):
        self.result = result
        self._warnings = list(warnings_text)
        self.cell_check.state(["disabled"] if result.kind == "molecule" else ["!disabled"])
        has_modes = bool(result.modes)
        self.save_btn.configure(state="normal" if has_modes else "disabled")
        self.report_btn.configure(state="normal")
        self._fill_result_texts()
        self._build_scene()
        if has_modes:
            first = next((i for i, m in enumerate(result.modes) if not m.acoustic), 0)
            self.tree.selection_set(str(first))      # → show_selected_mode が呼ばれる
            self.tree.see(str(first))
        else:
            self.comp_box.configure(values=[])
            self.comp_var.set("")
            self.comp_box.state(["disabled"])
            self.update_mode()

    def _fill_result_texts(self, keep_selection=False):
        """要約・モード一覧・レポート・状態表示を現在の言語で書く。"""
        result = self.result
        selection = self.tree.selection() if keep_selection else ()
        if result.kind == "molecule":
            head = t("gui.sum_molecule", pg=result.point_group,
                     lin=t("rep.mol.linear_tag") if result.is_linear else "")
        else:
            head = t("gui.sum_crystal", sg=result.space_group, pg=result.point_group)
        lines = [head, result.decomposition or t("gui.no_modes")] + \
                [t("rep.note", v=n).strip() for n in result.notes]
        self.summary_var.set("\n".join(lines))
        self.tree.delete(*self.tree.get_children())
        for i, m in enumerate(result.modes, 1):
            detail = ", ".join(m.tensor) if m.tensor else "—"
            if result.kind == "crystal":
                detail = f"{detail}   [{', '.join(m.contributions)}]"
            tag = "acoustic" if m.acoustic else m.activity
            direction = t("dir.acoustic") if m.acoustic else m.direction_label
            self.tree.insert("", "end", iid=str(i - 1), tags=(tag,),
                             values=(i, m.label, m.degeneracy, m.activity, direction, detail))
        self._fit_columns()
        if selection:
            self.tree.selection_set(selection)
            self.tree.see(selection[0])
        report = result.report
        if self._warnings:
            report += "\n\n" + t("gui.warnings") + "\n" + \
                "\n".join(f"  {w}" for w in self._warnings)
        self.report.configure(state="normal")
        self.report.delete("1.0", "end")
        self.report.insert("1.0", report)
        self.report.configure(state="disabled")
        n = sum(1 for m in result.modes if not m.acoustic)
        self.status_var.set(t("gui.done", n=len(result.numbers), m=n))

    def show_selected_mode(self):
        mode = self._selected_mode()
        if mode is None:
            return
        comps = [chr(97 + k) for k in range(len(mode.displacements))]
        if list(self.comp_box.cget("values")) != comps or self.comp_var.get() not in comps:
            self.comp_box.configure(values=comps)
            self.comp_var.set(comps[0])
        self.comp_box.state(["!disabled"] if len(comps) > 1 else ["disabled"])
        self.update_mode()

    def _selected_mode(self):
        sel = self.tree.selection()
        if not sel or self.result is None or not self.result.modes:
            return None
        return self.result.modes[int(sel[0])]

    # ── 3D 表示: 場面 (原子・結合) の構築 ───────────────────────────────
    def _build_scene(self):
        """構造が変わったときに一度だけ呼ぶ。表示する原子と描画物を作り直す。"""
        res = self.result
        if res.kind == "crystal":
            from .crystal import cell_atoms
            # 入力の単位胞の中 (境界上を含む) の原子。像は元の原子と同じ変位をする
            cart, _, parents = cell_atoms(res, boundary=True)
            cell = res.cell if res.cell is not None else res.lattice
        else:
            cart, parents, cell = res.cart, np.arange(len(res.cart)), None
        numbers = [res.numbers[p] for p in parents]
        self.scene = dict(
            cart=np.asarray(cart, float), parents=np.asarray(parents, int), cell=cell,
            bonds=bond_pairs(cart, numbers, molecule=res.kind == "molecule"),
            colors=[_atom_color(z) for z in numbers],
            sizes=[60 + 140 * min(covalent_radius(z), 2.0) for z in numbers],
            arrow_len=self._arrow_length(np.asarray(cart, float)))

        ax = self.ax
        self._stop_animation()
        ax.clear()
        ax.set_axis_off()                      # 座標軸・目盛は出さない (向きは左下の軸表示で)
        sc = self.scene
        cell_lines = []
        if cell is not None:
            a, b, c = cell
            corners = [np.zeros(3), a, b, c, a + b, a + c, b + c, a + b + c]
            for i, j in [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (2, 4), (2, 6), (3, 5),
                         (3, 6), (4, 7), (5, 7), (6, 7)]:
                p = np.array([corners[i], corners[j]])
                cell_lines += ax.plot(p[:, 0], p[:, 1], p[:, 2], color="#4a90d9", lw=1,
                                      ls="--")
        # 結合は 1 つのまとまりとして描く (原子数が多くても軽く動くように)
        bonds = add_bond_collection(ax, self._bond_segments(sc["cart"]))
        atoms = ax.scatter(sc["cart"][:, 0], sc["cart"][:, 1], sc["cart"][:, 2],
                           c=sc["colors"], s=sc["sizes"], edgecolors="#333333",
                           linewidths=0.6, depthshade=False)
        self._artists = dict(cell=cell_lines, bonds=bonds, atoms=atoms, quiver=None)
        self._toggle_cell(redraw=False)
        self.reset_view(redraw=False)
        self._draw_axes_indicator(cell)              # 視点を決めてから軸表示を描く

    def _bond_segments(self, pos):
        bonds = self.scene["bonds"]
        if not bonds:
            return np.zeros((0, 2, 3))
        idx = np.array(bonds)
        return np.stack([pos[idx[:, 0]], pos[idx[:, 1]]], axis=1)

    # ── 3D 表示: モードの表示・アニメーション ────────────────────────────
    def _current_arrows(self):
        """選択中のモード・成分の矢印 (表示原子ごと, Å)。モードがなければ None。"""
        mode = self._selected_mode()
        if mode is None or self.scene is None:
            return None, None
        k = max(0, min(ord((self.comp_var.get() or "a")[0]) - 97, len(mode.displacements) - 1))
        disp = np.asarray(mode.displacements[k], float)
        mx = np.max(np.linalg.norm(disp, axis=1))
        disp = disp / mx if mx > 1e-12 else disp
        arrows = disp[self.scene["parents"]] * self.scene["arrow_len"] * self.scale_var.get()
        idx = int(self.tree.selection()[0]) + 1
        suffix = f" ({chr(97 + k)})" if len(mode.displacements) > 1 else ""
        kind = t("gui.acoustic") if mode.acoustic else mode.activity
        guess = t("gui.plot_guess") if mode.refined else ""
        return arrows, t("gui.plot_mode", i=idx, label=mode.label, comp=suffix, kind=kind) + guess

    def update_mode(self):
        """モード・成分・倍率・表示設定が変わったときに呼ぶ (視点は保つ)。"""
        if self.scene is None:
            return
        self._stop_animation()
        arrows, title = self._current_arrows()
        self.ax.set_title(title or t("gui.plot_structure"), fontsize=10)
        if arrows is not None and self.anim_var.get():
            self._anim_gen += 1
            self._anim_frame = 0
            self._animate(self._anim_gen, arrows)
        else:
            self._place_atoms(self.scene["cart"], arrows)
            self.canvas.draw_idle()

    def _animate(self, gen, arrows):
        """アニメーションの 1 コマ。世代が変わっていたら (新しい表示に替わったら) 止まる。"""
        if gen != self._anim_gen:
            return
        phase = math.sin(2 * math.pi * self._anim_frame / _FRAMES)
        self._place_atoms(self.scene["cart"] + 0.5 * phase * arrows, arrows)
        self.canvas.draw_idle()
        self._anim_frame = (self._anim_frame + 1) % _FRAMES
        self._anim_job = self.after(_FRAME_MS, self._animate, gen, arrows)

    def _stop_animation(self):
        self._anim_gen += 1
        if self._anim_job is not None:
            self.after_cancel(self._anim_job)
            self._anim_job = None

    def _place_atoms(self, pos, arrows):
        """原子・結合の座標を更新し、矢印を描き直す (描画物は作り直さない)。"""
        art = self._artists
        art["atoms"]._offsets3d = (pos[:, 0], pos[:, 1], pos[:, 2])
        if art["bonds"] is not None:
            art["bonds"].set_segments(self._bond_segments(pos))
        if art["quiver"] is not None:
            art["quiver"].remove()
            art["quiver"] = None
        if arrows is not None and self.arrow_var.get():
            art["quiver"] = self.ax.quiver(pos[:, 0], pos[:, 1], pos[:, 2],
                                           arrows[:, 0], arrows[:, 1], arrows[:, 2],
                                           color=_ARROW_COLOR, linewidth=1.5,
                                           arrow_length_ratio=0.25)

    def _toggle_cell(self, redraw=True):
        if self._artists:
            show = self.cell_var.get()
            for line in self._artists["cell"]:
                line.set_visible(show)
        if redraw:
            self.canvas.draw_idle()

    def _draw_axes_indicator(self, cell):
        """左下に軸の向きを描く。結晶は a (赤), b (緑), c (青)、分子は x, y, z。"""
        ia = self.axes_ax
        ia.clear()
        ia.set_axis_off()
        if cell is not None:
            vecs, names = np.asarray(cell, float), ("a", "b", "c")
        else:
            vecs, names = np.eye(3), ("x", "y", "z")
        vecs = vecs / np.linalg.norm(vecs, axis=1, keepdims=True)
        for v, name, color in zip(vecs, names, ("#e02020", "#20b020", "#2040e0")):
            ia.quiver(0, 0, 0, *v, color=color, linewidth=2.5, arrow_length_ratio=0.3)
            # 軸名は結晶学の慣例どおりイタリック体 (a, b, c / x, y, z)
            ia.text(*(v * 1.3), name, color="black", fontsize=13, style="italic",
                    family=["Times New Roman", "DejaVu Serif"], ha="center", va="center")
        ia.scatter([0], [0], [0], color="#cccccc", s=25, edgecolors="#666666")
        for setter in (ia.set_xlim3d, ia.set_ylim3d, ia.set_zlim3d):
            setter(-1.1, 1.1)
        ia.set_box_aspect((1, 1, 1))
        ia.view_init(elev=self.ax.elev, azim=self.ax.azim, roll=self.ax.roll)

    def _sync_axes_indicator(self, event=None):
        """本体の視点が変わっていたら、左下の軸表示も同じ向きにする。"""
        ia, ax = self.axes_ax, self.ax
        if (ia.elev, ia.azim, ia.roll) != (ax.elev, ax.azim, ax.roll):
            ia.view_init(elev=ax.elev, azim=ax.azim, roll=ax.roll)
            # 描画中の draw_idle は無視されるので、今の描画が終わってから描き直す
            self.after_idle(self.canvas.draw_idle)

    def _draw_empty(self):
        self._stop_animation()
        self.axes_ax.clear()
        self.axes_ax.set_axis_off()
        self.ax.clear()
        self.ax.set_axis_off()
        self.ax.text2D(0.5, 0.5, t("gui.plot_empty"), ha="center", va="center",
                       transform=self.ax.transAxes, color="#888888")
        self.canvas.draw_idle()

    # ── 3D 表示: 視点操作 ───────────────────────────────────────────────
    def reset_view(self, redraw=True):
        """全体が収まる立方体の範囲と既定の向きに戻す。"""
        if self.scene is None:
            return
        pts = self.scene["cart"]
        pad = self.scene["arrow_len"] * 1.2
        if self.scene["cell"] is not None:
            a, b, c = self.scene["cell"]
            pts = np.vstack([pts, [np.zeros(3), a, b, c, a + b, a + c, b + c, a + b + c]])
        lo, hi = pts.min(0) - pad, pts.max(0) + pad
        center, half = (lo + hi) / 2, max((hi - lo).max() / 2, 1.0)
        self.ax.set_xlim3d(center[0] - half, center[0] + half)
        self.ax.set_ylim3d(center[1] - half, center[1] + half)
        self.ax.set_zlim3d(center[2] - half, center[2] + half)
        self.ax.set_box_aspect((1, 1, 1), zoom=1.35)   # 既定より大きめに表示
        axes = self.scene["cell"] if self.scene["cell"] is not None else np.eye(3)
        elev, azim = _best_view(axes)
        self.ax.view_init(elev=elev, azim=azim)
        if redraw:
            self.update_mode()

    def _on_scroll(self, event):
        """ホイールでズーム (表示範囲を中心から拡大・縮小)。"""
        if event.inaxes is not self.ax or self.scene is None:
            return
        factor = 0.85 if event.button == "up" else 1 / 0.85
        for get, set_ in ((self.ax.get_xlim3d, self.ax.set_xlim3d),
                          (self.ax.get_ylim3d, self.ax.set_ylim3d),
                          (self.ax.get_zlim3d, self.ax.set_zlim3d)):
            lo, hi = get()
            c, h = (lo + hi) / 2, (hi - lo) / 2 * factor
            set_(c - h, c + h)
        self.canvas.draw_idle()

    @staticmethod
    def _arrow_length(cart):
        """最大変位の矢印の長さ [Å]: 最短原子間距離の 0.6 倍 (構造の大きさに合わせる)。"""
        if len(cart) < 2:
            return 1.0
        d = np.linalg.norm(cart[:, None] - cart[None], axis=2)
        d[d < 1e-6] = np.inf
        return 0.6 * float(np.min(d))

    # ── 保存 ───────────────────────────────────────────────────────────
    def _ask_formats(self):
        """書き出す形式を選ぶ小さなダイアログ。("vesta",) / ("xsf",) / 両方、取消なら None。"""
        choices = {"vesta": ("vesta",), "xsf": ("xsf",), "both": ("xsf", "vesta")}
        var = tk.StringVar(value=_load_settings().get("save_format", "vesta"))
        if var.get() not in choices:
            var.set("vesta")
        answer = {}
        dlg = tk.Toplevel(self)
        dlg.title(t("gui.fmt_title"))
        dlg.transient(self.winfo_toplevel())
        dlg.resizable(False, False)
        body = ttk.Frame(dlg, padding=12)
        body.pack(fill="both", expand=True)
        ttk.Label(body, text=t("gui.fmt_prompt")).pack(anchor="w", pady=(0, 6))
        for key in ("vesta", "xsf", "both"):
            ttk.Radiobutton(body, text=t(f"gui.fmt_{key}"), variable=var,
                            value=key).pack(anchor="w")

        def ok():
            answer["v"] = var.get()
            dlg.destroy()

        buttons = ttk.Frame(body)
        buttons.pack(anchor="e", pady=(10, 0))
        ttk.Button(buttons, text=t("gui.ok"), command=ok).pack(side="left", padx=4)
        ttk.Button(buttons, text=t("gui.cancel"), command=dlg.destroy).pack(side="left")
        dlg.bind("<Return>", lambda e: ok())
        dlg.bind("<Escape>", lambda e: dlg.destroy())
        dlg.grab_set()
        dlg.wait_window()
        if "v" not in answer:
            return None
        _save_settings(save_format=answer["v"])
        return choices[answer["v"]]

    def save_modes(self):
        if not self.result or not self.result.modes:
            return
        formats = self._ask_formats()
        if not formats:
            return
        parent = filedialog.askdirectory(title=t("gui.save_dir_title"))
        if not parent:
            return
        outdir = os.path.join(parent, f"modes_{self.source_name}")
        if self.result.kind == "molecule":
            from .molecule import write_mode_files
        else:
            from .crystal import write_mode_files
        try:
            paths = write_mode_files(self.result, outdir, formats=formats)
        except OSError as e:
            messagebox.showerror(t("gui.title"), t("gui.save_failed", err=e))
            return
        self.status_var.set(t("gui.saved_status", n=len(paths), d=outdir))
        messagebox.showinfo(t("gui.title"), t("gui.saved_dialog", n=len(paths), d=outdir))

    def save_report(self):
        if not self.result:
            return
        path = filedialog.asksaveasfilename(
            title=t("gui.report_title"), defaultextension=".txt",
            initialfile=f"report_{self.source_name}.txt",
            filetypes=[(t("gui.ft_text"), "*.txt"), (t("gui.ft_all"), "*.*")])
        if not path:
            return
        try:
            with open(path, "w", encoding="utf-8") as f:
                f.write(self.report.get("1.0", "end"))
        except OSError as e:
            messagebox.showerror(t("gui.title"), t("gui.save_error", err=e))
            return
        self.status_var.set(t("gui.report_saved", path=path))

    def save_image(self):
        path = filedialog.asksaveasfilename(
            title=t("gui.image_title"), defaultextension=".png",
            initialfile=f"{self.source_name}.png",
            filetypes=[("PNG", "*.png"), ("SVG", "*.svg"), ("PDF", "*.pdf")])
        if path:
            try:
                self.fig.savefig(path, dpi=200)
            except (OSError, ValueError) as e:              # ValueError: 未対応の拡張子
                messagebox.showerror(t("gui.title"), t("gui.save_error", err=e))
                return
            self.status_var.set(t("gui.image_saved", path=path))


def main(argv=None):
    """GUI を起動する。引数にファイルを渡すと開いた状態で起動する。"""
    argv = sys.argv[1:] if argv is None else argv
    for mod in ("spgrep", "spglib", "pymatgen"):
        warnings.filterwarnings("ignore", category=DeprecationWarning, module=mod)
    if sys.platform == "win32":                         # 高 DPI 画面でぼやけないように
        try:
            import ctypes
            ctypes.windll.shcore.SetProcessDpiAwareness(1)
        except Exception:
            pass
    lang = _load_settings().get("language", "en")             # 既定は英語、前回の選択を記憶
    set_language(lang if lang in LANGUAGES else "en")
    root = tk.Tk()

    def report_error(exc, val, tb):
        # 画面操作中の想定外のエラー: exe では標準エラーがないので、ダイアログとログで知らせる
        log = _write_error_log("".join(traceback.format_exception(exc, val, tb)))
        messagebox.showerror(t("gui.title"), t("gui.unexpected", err=f"{exc.__name__}: {val}",
                                                               log=log))
    root.report_callback_exception = report_error
    root.title(t("gui.title"))
    root.geometry("1280x800")
    root.minsize(900, 600)
    style = ttk.Style(root)
    if "vista" in style.theme_names():
        style.theme_use("vista")
    app = VibrationApp(root)

    def on_close():
        app._stop_animation()
        root.quit()
        root.destroy()

    root.protocol("WM_DELETE_WINDOW", on_close)
    if argv:
        app.path_var.set(argv[0])
        ext = os.path.splitext(argv[0])[1].lower()
        app._set_kind("crystal" if ext == ".cif" else "molecule")
        root.after(200, app.run_analysis)
    root.mainloop()


if __name__ == "__main__":
    main()
