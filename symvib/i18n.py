"""表示言語 (英語 / 日本語) の切り替え。

使い方:
    from .i18n import t, Msg, set_language
    set_language("ja")              # 既定は英語 "en"
    t("rep.mol.molecule")           # → 現在の言語の文字列
    Msg("mol.symmetrized", pg="Td", shift=0.01)
        # → 表示するときに現在の言語で整形される文字列オブジェクト。
        #   解析結果の注記やエラーに使い、あとで言語を切り替えても正しく表示できる。

文字列はすべて MESSAGES に (英語, 日本語) の組で置く。
"""
from __future__ import annotations

import unicodedata

LANGUAGES = {"en": "English", "ja": "日本語"}
_current = "en"


def set_language(lang: str):
    global _current
    if lang not in LANGUAGES:
        raise ValueError(f"unknown language: {lang}")
    _current = lang


def get_language() -> str:
    return _current


def t(key: str, **kw) -> str:
    """キーに対応する現在の言語の文字列 (kw で {name} を埋める)。"""
    en, ja = MESSAGES[key]
    text = ja if _current == "ja" else en
    return text.format(**kw) if kw else text


class Msg:
    """表示時に現在の言語で整形される文字列。引数に Msg を入れ子にしてもよい。"""

    def __init__(self, key: str, **kw):
        self.key, self.kw = key, kw

    def __str__(self):
        return t(self.key, **self.kw)

    def __format__(self, spec):
        return format(str(self), spec)

    def __contains__(self, item):
        return item in str(self)

    def __repr__(self):
        return f"Msg({self.key!r})"


def width(text: str) -> int:
    """等幅表示での幅 (全角文字は 2)。"""
    return sum(2 if unicodedata.east_asian_width(c) in "WF" else 1 for c in str(text))


def pad(text, n: int, align: str = "<") -> str:
    """全角文字を考慮して幅 n にそろえる (表の桁合わせ用)。"""
    text = str(text)
    gap = max(n - width(text), 0)
    if align == ">":
        return " " * gap + text
    if align == "^":
        return " " * (gap // 2) + text + " " * (gap - gap // 2)
    return text + " " * gap


MESSAGES = {
    # ── 共通 ──
    "dir.exact": ("symmetry-exact", "対称性確定"),
    "dir.guess": ("guess (network)", "推測(網目)"),
    "dir.acoustic": ("acoustic", "音響"),
    "char.iso": ("isotropic", "等方的"),
    "char.z": ("z (out-of-plane)", "z(面外)"),
    "char.xy": ("xy (in-plane)", "xy(面内)"),
    "char.mixed": ("mixed", "混合"),
    "rep.note": ("  * {v}", "  ※ {v}"),
    "rep.activity_line": ("    {sym} {act}: {n} modes ({dof} DOF)",
                          "    {sym} {act}: {n} モード ({dof} 自由度)"),
    "err.partial_site": ("Structures with partially occupied sites ({sites}) cannot be analysed. "
                         "Please provide an ordered structure.",
                         "部分占有サイト ({sites}) を含む構造は解析できません。"
                         "秩序化した構造を入力してください。"),

    # ── CIF の読み込み ──
    "cif.symops_from_sg": ("No symmetry operations were listed; used the standard-setting "
                           "operations of space group {sg}.",
                           "対称操作の記載がないため、空間群 {sg} の標準設定の操作を使いました。"),
    "cif.p1": ("No symmetry operations or space group were listed; read as P1.",
               "対称操作・空間群の記載がないため P1 として読みました。"),
    "cif.fraction": ("{n} coordinate value(s) written as rounded fractions (e.g. 0.333, 0.6667) "
                     "were read as the exact fractions (1/3, 2/3, ...).",
                     "丸めた分数の小数 (0.333, 0.6667 など) で書かれた座標 {n} 個を、"
                     "正確な分数 (1/3, 2/3 など) として読みました。"),
    "cif.multi": ("The CIF contains {n} structures; used the first one, '{name}'.",
                  "CIF に {n} 個の構造があるため、最初の '{name}' を使いました。"),
    "cif.unknown": ("Sites whose element could not be identified ({sites}) were treated as "
                    "pseudo-atoms of mass 1 amu.",
                    "元素を特定できないサイト ({sites}) は質量 1 amu の仮想原子として扱いました。"),
    "cif.merged": ("Sites sharing the same position were merged into one pseudo-atom: {sites} "
                   "(mass averaged by occupancy).",
                   "同じ座標を占める複数のサイトを 1 つの仮想原子にまとめました: {sites} "
                   "(質量は占有率で平均)"),
    "cif.partial": ("Partially occupied sites were treated as fully occupied: {sites}",
                    "部分占有のサイトを完全占有として扱いました: {sites}"),
    "err.cif_read": ("Cannot read the CIF file: {path}\n  {err}",
                     "CIF ファイルを読めません: {path}\n  {err}"),
    "err.cif_no_atoms": ("The CIF has no cell parameters and fractional atomic coordinates "
                         "(_atom_site_fract_x).",
                         "CIF に格子定数と原子の分率座標 (_atom_site_fract_x) がありません。"),
    "err.cif_cell": ("Cannot read the cell parameters (_cell_length_*, _cell_angle_*) in the CIF.",
                     "CIF の格子定数 (_cell_length_*, _cell_angle_*) が読めません。"),
    "err.cif_coords": ("Cannot read the fractional coordinates of atom {label}.",
                       "原子 {label} の分率座標が読めません。"),

    # ── 結晶 ──
    "err.reduce": ("Could not map the irreducible representations onto the primitive-cell "
                   "operations (try a different symmetry tolerance, symprec).",
                   "既約表現を primitive cell の対称操作に対応づけられません "
                   "(対称性の許容差 symprec を変えて再実行してください)。"),
    "err.image_atom": ("A symmetry operation maps an atom onto no atom "
                       "(check the symmetry tolerance, symprec).",
                       "対称操作で写った先の原子が見つかりません "
                       "(対称性の許容差 symprec を見直してください)。"),
    "err.no_sg": ("Could not determine the space group (please check the structure).",
                  "空間群を判定できません (構造を確認してください)。"),
    "cry.idealized": ("The lattice and coordinates were idealized to the detected symmetry "
                      "(max. atom shift {shift:.2e} Å, max. lattice strain {strain:.1e}).",
                      "格子と座標を検出した対称性に合わせて理想化しました "
                      "(原子の最大移動 {shift:.2e} Å, 格子の最大ひずみ {strain:.1e})。"),
    "cry.centered": ("Centred lattice → analysed in the primitive cell ({n0} → {n} atoms)",
                     "中心格子 → primitive cell で解析 ({n0} → {n} 原子)"),
    "cry.inversion": ("Inversion centre present → mutual exclusion rule (Raman ↔ IR)",
                      "反転対称あり → 相互排他律 (Raman ↔ IR)"),
    "rep.cry.raman_note": ("  (x, y, z = axes of the standard setting of the space group; "
                           "letters = independent parameters, · = zero by symmetry.\n"
                           "   For degenerate modes one matrix is given per mutually orthogonal "
                           "component; any combination within the degenerate space is allowed.)",
                           "  (x, y, z は空間群の標準設定の軸。文字は独立なパラメータ、· は対称性で"
                           "ゼロの成分。\n   縮退モードは互いに直交する成分ごとに 1 つの行列を示す。"
                           "縮退空間内の任意の組み合わせが許される。)"),
    "rep.cry.tensor_form_k": ("    Raman tensor (component {k}): {v}",
                              "    ラマンテンソル (成分 {k}): {v}"),
    "rep.cry.title": (" Factor Group Analysis — Γ point", " 因子群解析 (Factor Group Analysis) — Γ point"),
    "rep.cry.structure": ("■ Crystal structure", "■ 結晶構造"),
    "rep.cry.sg": ("  Space group: {v}", "  空間群: {v}"),
    "rep.cry.pg": ("  Point group: {v} (= factor group)", "  点群:   {v} (= 因子群)"),
    "rep.cry.natoms": ("  Atoms: {n}", "  原子数: {n}"),
    "rep.cry.nops": (",  symmetry operations: {n}", ",  対称操作数: {n}"),
    "rep.cry.decomp": ("■ Irreducible decomposition (3N = {n})", "■ 既約分解 (3N = {n})"),
    "rep.cry.acoustic": ("  Acoustic: {n} ({dof} DOF)", "  音響: {n} ({dof} 自由度)"),
    "rep.cry.optical": ("  Optical: {n} ({dof} DOF)", "  光学: {n} ({dof} 自由度)"),
    "rep.cry.cols": ("Irrep|d|Activity|Displacement|Type|Direction|Sites",
                     "既約表現|d|活性|変位|種別|方向|関与サイト"),
    "rep.cry.legend": (" d = degeneracy, ♪ = acoustic mode\n"
                       " Direction: symmetry-exact  = unique from group theory alone (exact) /\n"
                       "            guess (network) = the irrep occurs more than once (m>1), so the\n"
                       "                              direction is a 1/d² spring-network guess "
                       "(true eigenvector needs DFT)",
                       " d = 縮退度, ♪ = 音響モード\n"
                       " 方向: 対称性確定 = 群論だけで一意 (厳密) /\n"
                       "       推測(網目) = 同一既約表現が複数回出現 (m>1) のため、\n"
                       "                    1/d² バネ網目モデルによる推測方向 (真の向きは要 DFT)"),
    "rep.cry.raman_title": (" ■ Polarized Raman selection rules (Raman tensors)",
                            " ■ 偏光ラマン選択則 (Raman テンソル)"),
    "rep.cry.raman_mode": ("  {label} (degeneracy {d}):", "  {label} (縮退度 {d}):"),
    "rep.cry.tensor_comps": ("    Tensor components: {v}", "    テンソル成分: {v}"),
    "rep.cry.tensor_form": ("    Raman tensor: {v}", "    ラマンテンソル: {v}"),
    "rep.cry.porto": ("    Porto geometry   (✓ = observed, ✗ = extinct):",
                      "    Porto 配置     (✓=観測, ✗=消光):"),

    # ── 分子 ──
    "err.read_struct": ("Cannot read the structure file: {path}\n  as a molecule: {e1}\n"
                        "  as a crystal: {e2}",
                        "構造ファイルを読めません: {path}\n  分子として: {e1}\n  結晶として: {e2}"),
    "err.xyz_rows": ("{n} atoms are declared but there are only {m} coordinate lines",
                     "原子数 {n} に対して座標が {m} 行しかありません"),
    "err.xyz_format": ("Invalid xyz file: {path}\n  {err}",
                       "xyz ファイルの形式が正しくありません: {path}\n  {err}"),
    "err.linear_internal": ("The symmetry operations of a linear molecule are not finite "
                            "(internal error).",
                            "直線分子の対称操作は有限個に定まりません (内部エラー)。"),
    "err.min_atoms": ("At least 2 atoms are required.", "原子が 2 個以上必要です。"),
    "err.order_mismatch": ("Point group {pg} (order {order}) does not match the number of "
                           "symmetry operations found ({n}). Try a different tolerance.",
                           "点群 {pg} (位数 {order}) と見つかった対称操作の数 ({n}) が"
                           "一致しません。許容差を変えて再実行してください。"),
    "axis.c2v_planar": ("z = C2 axis, x = normal of the molecular plane",
                        "z = C2 軸, x = 分子面の法線"),
    "axis.c2v_other": ("z = C2 axis, x = normal of the mirror plane containing more atoms",
                       "z = C2 軸, x = 原子を多く含む鏡面の法線"),
    "axis.d2_planar": ("z = C2 axis through the most atoms, x = normal of the molecular plane",
                       "z = 最も多くの原子を通る C2 軸, x = 分子面の法線"),
    "axis.d2_other": ("z = C2 axis through the most atoms, x = C2 axis through the next most atoms",
                      "z = 最も多くの原子を通る C2 軸, x = 次に多くの原子を通る C2 軸"),
    "mol.axes": ("Axis convention: {frame} (B1/B2/B3 labels depend on this choice)",
                 "軸の取り方: {frame} (B1/B2/B3 のラベルはこの軸に依存)"),
    "mol.symmetrized": ("The structure was symmetrized exactly to point group {pg} "
                        "(max. atom shift {shift:.4f} Å).",
                        "構造を点群 {pg} に厳密に対称化しました (原子の最大移動 {shift:.4f} Å)。"),
    "mol.linear": ("Linear molecule: infinite group {g} analysed exactly by continuous-group "
                   "projection. Σ = along the axis, Π = bending (2-fold degenerate)",
                   "直線分子: 無限群 {g} を連続群射影で厳密に解析。"
                   " Σ = 軸方向, Π = 変角 (2 重縮退)"),
    "mol.no_table": ("No character table for point group {pg}, so the vibrations cannot be "
                     "decomposed (number of vibrations: {n}). Available: {avail}",
                     "点群 {pg} の指標表は未収録のため、既約分解できません "
                     "(振動の本数は {n})。収録済み: {avail}"),
    "mol.tensor_frame": ("Raman tensor components (xx, xz, …) refer to the standard orientation "
                         "of the point group (z = principal axis; for C2v/D2/D2h the axes stated "
                         "above), not to the input coordinates. Displacements are in the input "
                         "coordinates.",
                         "ラマンテンソル成分 (xx, xz, …) は点群の標準の向き (z = 主軸; C2v/D2/D2h "
                         "は上記の軸) で表しており、入力座標の向きとは別です。変位は入力座標系です。"),
    "mol.arbitrary": ("{labels}: the spring network (central forces) has zero stiffness in "
                      "this subspace, so the guessed direction is arbitrary (e.g. out-of-plane "
                      "modes of planar molecules).",
                      "{labels} は網目モデル (中心力バネ) の剛性がゼロの部分空間のため、"
                      "推測方向は任意です (平面分子の面外モードなど)。"),
    "rep.mol.title": (" Molecular Vibration Analysis — no DFT required",
                      " 分子振動モード解析 (Molecular Vibration Analysis) — DFT 不要"),
    "rep.mol.molecule": ("■ Molecule", "■ 分子"),
    "rep.mol.atoms": ("  Atoms: {n}   Point group: {pg}{lin}", "  原子数: {n}   点群: {pg}{lin}"),
    "rep.mol.linear_tag": ("  (linear)", "  (直線分子)"),
    "rep.mol.dof": ("  DOF: 3N={n3}  (translation 3 + rotation {r} + vibration {v})",
                    "  自由度: 3N={n3}  (並進3 + 回転{r} + 振動{v})"),
    "rep.mol.internal": ("■ Internal vibrations ({dof} DOF)", "■ 内部振動モード ({dof} 自由度)"),
    "rep.mol.cols": ("No.|Irrep|Deg.|Activity|Direction|Tensor components",
                     "振動|既約表現|縮退|活性|方向|テンソル成分"),
    "rep.mol.legend": (" Direction: symmetry-exact  = unique from group theory (exact) /\n"
                       "            guess (network) = the irrep occurs more than once (m>1);\n"
                       "                              spring-network guess "
                       "(true direction needs quantum chemistry)",
                       " 方向: 対称性確定 = 群論で一意 (厳密) /\n"
                       "       推測(網目) = 同一既約表現が複数回 (m>1) のためバネ網目モデル推測\n"
                       "                    (真の向きは量子化学計算が必要)"),

    # ── QE 凍結フォノン ──
    "err.qe_delta": ("Displacement amplitudes must be finite and positive, one per basis vector.",
                     "変位振幅は基底ベクトルごとに有限の正の値を指定してください。"),
    "err.qe_pseudo": ("Structures containing merged partial-occupancy sites or unknown elements "
                      "cannot be used for QE frozen-phonon calculations "
                      "(please provide an ordered structure).",
                      "部分占有をまとめた仮想原子や未知の元素を含む構造は、QE の凍結フォノン計算に"
                      "使えません (秩序化した構造を入力してください)。"),
    "err.qe_pp_missing": ("The config 'pseudopotentials' has no entry for: {el}",
                          "config の pseudopotentials に {el} がありません。"),
    "err.qe_no_output": ("pw.x output not found: {f}", "pw.x の出力がありません: {f}"),
    "err.qe_forces": ("Cannot read forces: {f} ({k}/{n} atoms)", "力が読めません: {f} ({k}/{n} 原子)"),
    "err.qe_no_manifest": ("manifest.json not found: {p} (run 'generate' first)",
                           "manifest.json がありません: {p} (先に generate を実行)"),
    "qe.gen_title": ("■ Frozen-phonon input generation: {sg} ({pg})",
                     "■ 凍結フォノン入力生成: {sg} ({pg})"),
    "qe.centered": ("  * Centred lattice → primitive cell ({n} atoms)",
                    "  ※ 中心格子 → primitive cell ({n} 原子)"),
    "qe.target.selected": ("selected", "指定"),
    "qe.target.raman": ("Raman-active", "ラマン活性"),
    "qe.target.optical_all": ("all optical", "全光学"),
    "qe.target": ("  Targets: {target} irreps ({n})", "  対象: {target} 既約表現  {n} 種"),
    "qe.block": ("    {label} (d={d}, m={m}, {act}) → {n} displacements",
                 "    {label} (d={d}, m={m}, {act}) → {n} 変位"),
    "qe.inputs": ("  pw.x inputs: {n} (±δ, δ={delta} Å) + 1 reference",
                  "  pw.x 入力: {n} 本 (±δ, δ={delta} Å) + 参照 1 本"),
    "qe.no_config": ("  (no config given → manifest only, no input files written)",
                     "  (config 未指定 → 入力ファイルは書かず manifest のみ)"),
    "qe.outdir": ("  Output: {d}/", "  出力先: {d}/"),
    "qe.an_title": (" Force-based Γ-point block frozen phonon — frequencies",
                    " 力ベース Γ点ブロック凍結フォノン — 振動数"),
    "qe.an_head": ("■ {sg} ({pg}),  {n} atoms,  δ={delta} Å", "■ {sg} ({pg}),  {n} 原子,  δ={delta} Å"),
    "qe.an_count": ("  {n} optical modes ({dof} DOF; {ac} acoustic DOF excluded)",
                    "  光学モード {n} 本 (自由度 {dof}; 音響 {ac} 自由度は除外)"),
    "qe.an_cols": ("#|Irrep|Freq. [cm⁻¹]|Deg.|Activity|Tensor components",
                   "#|既約表現|振動数[cm⁻¹]|縮退|活性|テンソル成分"),
    "qe.deg": ("{d}-fold", "{d}重"),
    "qe.soft": ("  ⚠ {n} imaginary frequencies (soft modes / instability): the structure may be "
                "unrelaxed or unstable towards a phase transition.",
                "  ⚠ 虚振動数 (ソフトモード/不安定) が {n} 本。構造が未緩和か、相転移不安定性の可能性。"),
    "qe.xsf": ("  Eigenvector XSF files: {d}/  ({n})", "  固有ベクトル XSF: {d}/  ({n} 本)"),
    "qe.validate": ("full diagonalization ↔ block method, max. difference",
                    "全対角化↔ブロック法 最大差"),

    # ── CLI ──
    "cli.mol_demo": ("No input → demo (water molecule)", "入力なし → デモ (水分子)"),
    "cli.cry_demo": ("No CIF given → demo (graphite 2H)", "CIF 未指定 → デモ (黒鉛 2H)"),
    "cli.output": ("■ Output: {d}/  ({n} files: .xsf and .vesta)",
                   "■ 出力: {d}/  ({n} ファイル: .xsf と .vesta)"),
    "cli.vesta_hint": ("  .vesta opens directly in VESTA; for .xsf use Properties → Vectors → "
                       "Show vectors",
                       "  .vesta は VESTA でそのまま開けます。.xsf は Properties → Vectors → "
                       "Show vectors で矢印を表示"),
    "cli.components": ("  Degenerate modes are written per component (_a, _b, _c, ...)",
                       "  縮退モードは成分ごとに _a, _b, _c... と出力"),
    "err.not_group": ("The symmetry operations found do not form a group (try a different "
                      "symmetry tolerance).",
                      "見つかった対称操作が群になっていません (対称性の許容差を変えてみてください)。"),
    "err.class_assign": ("Could not match the symmetry operations to the {pg} character table "
                         "(try a different symmetry tolerance).",
                         "対称操作を {pg} の指標表のクラスに対応させられませんでした "
                         "(対称性の許容差を変えてみてください)。"),
    "err.projection": ("Internal check failed: the symmetry projection is inconsistent. The result "
                       "would be unreliable, so the analysis was stopped.",
                       "内部検査で不整合を検出しました (対称性の射影が不正)。結果が信頼できない"
                       "ため解析を中止しました。"),
    "err.symprec": ("The symmetry tolerance must be a positive number.",
                    "対称性の許容差は正の数にしてください。"),
    "err.no_sg_detail": ("Space-group determination failed (check for overlapping atoms or the "
                         "symmetry tolerance).\n  {err}",
                         "空間群を判定できませんでした (原子の重なりや許容差を確認してください)。"
                         "\n  {err}"),
    "cry.std_axes": ("The input cell is not in the standard setting of {sg}; irreducible-"
                     "representation labels and Raman tensor components (x, y, z) refer to the "
                     "standard setting (x ∥ a, z ∥ c of the standard cell).",
                     "入力のセルは {sg} の標準の設定ではありません。既約表現のラベルとラマン"
                     "テンソル成分 (x, y, z) は標準の設定 (標準セルの x ∥ a, z ∥ c) で表しています。"),
    "err.tol": ("The symmetry tolerance must be a positive number.",
                "対称性の許容差は正の数にしてください。"),
    "err.read_crystal": ("Cannot read the file as a crystal structure: {path}\n  {err}",
                         "結晶構造として読めません: {path}\n  {err}"),
    "err.qe_cfg_key": ("The QE settings JSON lacks required keys: {keys}",
                       "QE 設定 JSON に必須の項目がありません: {keys}"),
    "gui.plot_mode": ("Mode {i}: {label}{comp}  {kind}", "モード {i}: {label}{comp}  {kind}"),
    "gui.plot_guess": ("  [direction: network guess]", "  [方向: 網目モデルによる推測]"),
    "gui.plot_structure": ("Structure", "構造"),
    "gui.plot_empty": ("No structure", "構造がありません"),
    "gui.acoustic": ("acoustic", "音響"),
    "gui.unexpected": ("An unexpected error occurred:\n{err}\n\nDetails were saved to {log}",
                       "予期しないエラーが発生しました:\n{err}\n\n詳細を {log} に保存しました"),
    "gui.save_error": ("Could not save the file:\n{err}", "保存できませんでした:\n{err}"),
    "cli.error": ("Error: {e}", "エラー: {e}"),

    # ── GUI ──
    "gui.title": ("SymVib — Vibrational Mode Analysis", "SymVib — 振動モード解析"),
    "gui.input": ("Input:", "入力:"),
    "gui.open": ("Open…", "開く…"),
    "gui.demo": ("Demo", "デモ"),
    "gui.demo_water": ("Water H₂O (molecule)", "水 H₂O (分子)"),
    "gui.demo_graphite": ("Graphite 2H (crystal)", "黒鉛 2H (結晶)"),
    "gui.molecule": ("Molecule", "分子"),
    "gui.crystal": ("Crystal (Γ point)", "結晶 (Γ点)"),
    "gui.tol_mol": ("Symmetry tolerance [Å]:", "対称性の許容差 [Å]:"),
    "gui.tol_cry": ("symprec [Å]:", "symprec [Å]:"),
    "gui.run": ("▶ Analyze", "▶ 解析"),
    "gui.save_modes": ("Save XSF/VESTA…", "XSF/VESTA 保存…"),
    "gui.save_report": ("Save report…", "レポート保存…"),
    "gui.fmt_title": ("Save mode files", "モードファイルの保存"),
    "gui.fmt_prompt": ("File format to write:", "書き出す形式:"),
    "gui.fmt_vesta": ("VESTA only (.vesta: gray arrows and bonds preset)",
                      "VESTA のみ (.vesta: 灰色の矢印・結合の設定込み)"),
    "gui.fmt_xsf": ("XSF only (.xsf)", "XSF のみ (.xsf)"),
    "gui.fmt_both": ("Both", "両方"),
    "gui.ok": ("OK", "OK"),
    "gui.cancel": ("Cancel", "キャンセル"),
    "gui.language": ("Language:", "言語:"),
    "gui.start": ("Open a structure file or choose a demo.", "ファイルを開くか、デモを選んでください。"),
    "gui.cols": ("#|Irrep|Deg.|Activity|Direction|Raman tensor / sites",
                 "#|既約表現|縮退|活性|方向|Raman テンソル / 関与サイト"),
    "gui.legend": ("  Direction: symmetry-exact = from group theory / guess = estimated (m>1)",
                   "  方向: 対称性確定 = 群論で厳密 / 推測(網目) = m>1 のため推測"),
    "gui.hint": ("Left-drag: rotate   Wheel: zoom   Right-drag: pan",
                 "左ドラッグ: 回転   ホイール: ズーム   右ドラッグ: 平行移動"),
    "gui.save_image": ("Save image…", "画像保存…"),
    "gui.reset_view": ("Reset view", "視点リセット"),
    "gui.component": ("Component:", "縮退成分:"),
    "gui.arrow_scale": ("Arrow scale:", "矢印の倍率:"),
    "gui.arrows": ("Arrows", "矢印"),
    "gui.animation": ("Animation", "アニメーション"),
    "gui.cell": ("Unit cell", "単位胞"),
    "gui.err_tol": ("Please enter a positive number for the tolerance.",
                    "許容差には正の数を入力してください。"),
    "gui.err_nofile": ("File not found:\n{path}", "ファイルが見つかりません:\n{path}"),
    "gui.running": ("Analyzing…", "解析中…"),
    "gui.failed": ("Analysis failed", "解析に失敗しました"),
    "gui.failed_dialog": ("The analysis could not be completed.\n\n{err}",
                          "解析できませんでした。\n\n{err}"),
    "gui.sum_molecule": ("Molecule   point group {pg}{lin}", "分子  点群 {pg}{lin}"),
    "gui.sum_crystal": ("Crystal   space group {sg}   point group {pg}",
                        "結晶  空間群 {sg}  点群 {pg}"),
    "gui.no_modes": ("(no vibrational modes)", "(振動モードなし)"),
    "gui.warnings": ("■ Warnings while reading the file", "■ 読み込み時の警告"),
    "gui.done": ("Done: {n} atoms, {m} vibrational modes", "解析完了: {n} 原子, 振動モード {m} 組"),
    "gui.open_title": ("Open structure file", "構造ファイルを開く"),
    "gui.ft_struct": ("Structure files", "構造ファイル"),
    "gui.ft_xyz": ("xyz (molecule)", "xyz (分子)"),
    "gui.ft_all": ("All files", "すべて"),
    "gui.ft_text": ("Text", "テキスト"),
    "gui.save_dir_title": ("Choose a folder for the XSF / VESTA files",
                           "XSF / VESTA ファイルの保存先フォルダを選択"),
    "gui.save_failed": ("Could not save:\n{err}", "保存できませんでした:\n{err}"),
    "gui.saved_status": ("Saved {n} files: {d}", "{n} ファイル保存: {d}"),
    "gui.saved_dialog": ("Saved {n} files.\n{d}\n\n"
                         "• .vesta: opens directly in VESTA (thin grey arrows and bonds included)\n"
                         "• .xsf: open in VESTA and use Properties → Vectors → Show vectors",
                         "{n} ファイルを保存しました。\n{d}\n\n"
                         "・.vesta: VESTA でそのまま開けます (灰色の細い矢印・結合の設定込み)\n"
                         "・.xsf : VESTA で開き、Properties → Vectors → Show vectors で矢印を表示"),
    "gui.report_title": ("Save report", "レポートを保存"),
    "gui.report_saved": ("Report saved: {path}", "レポートを保存: {path}"),
    "gui.image_title": ("Save figure", "図を保存"),
    "gui.image_saved": ("Figure saved: {path}", "図を保存: {path}"),
}
