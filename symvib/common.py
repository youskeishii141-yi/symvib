"""分子・結晶の両解析で共有するデータ型と小さな道具。

- Mode / AnalysisResult : 解析結果の構造化データ (CLI 表示と GUI の両方が使う)
- atomic_mass / element_symbol : 元素情報
- network_hessian / refine_in_subspace : m>1 モードの向きを推測する 1/d² バネ網目モデル
- write_xsf / safe_label : VESTA 用 XSF 出力
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from functools import lru_cache

import numpy as np

from .i18n import Msg, t

# 活性の表示順と記号 (表の集計・GUI の色分けに使う)
ACTIVITIES = ("Raman", "IR", "Raman+IR", "Silent")
ACTIVITY_SYMBOL = {"Raman": "🔴", "IR": "🔵", "Raman+IR": "🟣", "Silent": "⚪"}


def setup_console():
    """CLI の標準出力を UTF-8 にし、依存ライブラリ内部の非推奨警告を隠す。

    日本語 Windows では出力をパイプ/リダイレクトすると cp932 になり、Σ や ✓ を
    含むレポートで UnicodeEncodeError になるため。警告は spgrep が内部で使う
    spglib の旧 API に関するもので、利用者には対処できない。"""
    import sys
    import warnings
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(encoding="utf-8", errors="replace")
    for mod in ("spgrep", "spglib", "pymatgen"):
        warnings.filterwarnings("ignore", category=DeprecationWarning, module=mod)


class AnalysisError(Exception):
    """入力構造や対称性の問題で解析を続けられないときに送出する。

    メッセージは利用者向け (日本語) で、GUI ではそのままダイアログに出す。"""


@dataclass(eq=False)          # 配列を持つので == は同一オブジェクトかどうかで判定する
class Mode:
    """1 つの振動モード (縮退している場合は縮退多重項全体)。"""
    label: str                      # 既約表現 (Mulliken 記号, 例 "A1g", "E'")
    degeneracy: int                 # 縮退度 d
    activity: str                   # "Raman" / "IR" / "Raman+IR" / "Silent"
    tensor: list[str]               # Raman テンソル成分 (例 ["xx", "yy"])
    displacements: list[np.ndarray]  # 縮退成分ごとの Cartesian 変位 (N, 3)
    refined: bool = False           # True: m>1 のため網目モデルによる推測方向
    acoustic: bool = False          # True: 音響モード (結晶のみ)
    character: str = ""             # 変位の性格 "z" / "xy" / "mixed" / "iso" (結晶のみ)
    contributions: list[str] = field(default_factory=list)  # 関与するサイト (例 "Al(c)") — 結晶のみ
    # Raman テンソルの形 (結晶のみ)。縮退成分 (互いに直交する組) ごとの 3×3 行列の文字列表現
    # (例 "[[a,·,·],[·,a,·],[·,·,b]]")。文字は独立なパラメータ、· は対称性でゼロの成分。
    raman_tensors: list[str] = field(default_factory=list)

    @property
    def direction_label(self) -> str:
        """表示用: 方向の確かさ (現在の言語)。"""
        if self.acoustic:
            return "—"
        return t("dir.guess") if self.refined else t("dir.exact")

    @property
    def character_label(self) -> str:
        """表示用: 変位の性格 (現在の言語)。"""
        return t(f"char.{self.character}") if self.character else ""


@dataclass
class AnalysisResult:
    """解析結果一式。kind は "molecule" または "crystal"。"""
    kind: str
    numbers: list[int]              # 原子番号
    symbols: list[str]              # 元素記号
    cart: np.ndarray                # 表示・出力に使う Cartesian 座標 [Å] (N, 3)
    modes: list[Mode]
    point_group: str = ""
    space_group: str = ""           # 結晶のみ (例 "P6_3/mmc (No. 194)")
    lattice: np.ndarray | None = None   # 結晶のみ: 解析に使った primitive cell (行ベクトル, Å)
    cell: np.ndarray | None = None      # 結晶のみ: 入力ファイルの単位胞 (表示・VESTA 出力用)
    is_linear: bool = False         # 分子のみ
    decomposition: str = ""         # 例 "Γ = A1g + 2E2g + ..."
    notes: list = field(default_factory=list)   # 利用者への補足・警告 (str または i18n.Msg)
    meta: dict = field(default_factory=dict)    # レポート作成用の補助情報 (Wyckoff 記号など)

    @property
    def optical_modes(self) -> list[Mode]:
        return [m for m in self.modes if not m.acoustic]

    @property
    def report(self) -> str:
        """CLI と同じテキストレポート (呼んだ時点の表示言語で作る)。"""
        if self.kind == "molecule":
            from .molecule import format_report
            return format_report(self)
        from .crystal import format_report
        return format_report(self, **self.meta)


# ════════════════════════════════════════════════════════════════════════
# 元素情報
# ════════════════════════════════════════════════════════════════════════
@lru_cache(maxsize=None)
def atomic_mass(z: int) -> float:
    """原子番号 z の原子量 [amu] (仮想原子なら登録した質量)。"""
    if int(z) in _PSEUDO:
        return _PSEUDO[int(z)]["mass"]
    from pymatgen.core import Element
    return float(Element.from_Z(int(z)).atomic_mass)


@lru_cache(maxsize=None)
def element_symbol(z: int) -> str:
    """原子番号 z の元素記号 (仮想原子ならサイト名, 例 "Ti1V2")。"""
    if int(z) in _PSEUDO:
        return _PSEUDO[int(z)]["name"]
    from pymatgen.core import Element
    return Element.from_Z(int(z)).symbol


# ── 仮想原子 (部分占有サイト・未知の元素記号) ──────────────────────────────
# 解析では原子の「種類」を整数で区別する。実在元素は原子番号 (1〜118)、
# それ以外は PSEUDO_BASE 以上の番号を割り当て、名前・質量・表示用の元素を登録する。
PSEUDO_BASE = 1000
_PSEUDO: dict[int, dict] = {}


def register_pseudo_species(name, mass, display_z, key=None):
    """仮想原子を登録して番号を返す (同じ key なら同じ番号)。

    name      : 表示名 (例 "Ti1V2")
    mass      : 質量 [amu] (占有率で平均したもの)
    display_z : XSF/VESTA に書く・色を決めるための代表元素の原子番号
    key       : 同一種とみなす判定キー (既定は name)"""
    key = key if key is not None else name
    for z, info in _PSEUDO.items():
        if info["key"] == key:
            return z
    z = PSEUDO_BASE + len(_PSEUDO)
    _PSEUDO[z] = dict(name=name, mass=float(mass), display_z=int(display_z), key=key)
    return z


def is_pseudo(z) -> bool:
    return int(z) in _PSEUDO


def display_z(z) -> int:
    """ファイル出力・色分け用の原子番号 (仮想原子は代表元素)。"""
    return _PSEUDO[int(z)]["display_z"] if int(z) in _PSEUDO else int(z)


def site_element(site) -> tuple[int, str]:
    """pymatgen の site から (原子番号, 元素記号) を返す。

    CIF に酸化数があると site.specie は "O2-" のような Species になるので、
    str() ではなく .symbol で純粋な元素記号を取り出す。部分占有サイトは
    振動解析の対象外なので明示的なエラーにする。"""
    if not site.is_ordered:
        raise AnalysisError(Msg("err.partial_site", sites=site.species_string))
    sp = site.specie
    return int(sp.Z), sp.symbol


# ════════════════════════════════════════════════════════════════════════
# 1/d² バネ網目モデル (m>1 モードの向きの推測)
# ════════════════════════════════════════════════════════════════════════
def network_hessian(cart, lattice=None, power=2.0, cutoff=None, k0=1.0, isotropic=False):
    """1/d^power 減衰の中心力バネによる力定数行列 Φ (3N×3N, 質量重み前)。

    lattice を与えると周期系として ±1 セルの像も数える (Γ 点用)。
    isotropic=True なら結合方向によらない等方バネ (どの向きの相対変位にも抵抗する)。
    中心力では剛性ゼロになる部分空間 (平面分子の面外など) の向きを決める補助に使う。
    これは「結合が短いほど強い」という幾何学的な当て推量であり、電子状態が
    決める真の力定数ではない。m>1 モードの向きの目安専用。"""
    cart = np.asarray(cart, float)
    na = len(cart)
    phi = np.zeros((3 * na, 3 * na))
    if lattice is None:
        shifts = [np.zeros(3)]
    else:
        # カットオフ内の像をすべて数えるのに必要なセル数を、格子面の間隔 d_i = V/|a_j × a_k| から決める
        # (±1 セルに固定すると、扁平な primitive cell (体心・面心・菱面体) で一部の近接原子を取りこぼし、
        #  網目の力定数が対称操作と可換でなくなる)
        lattice = np.asarray(lattice, float)
        if cutoff is None:
            reach = (1, 1, 1)
        else:
            vol = abs(np.linalg.det(lattice))
            spacing = [vol / np.linalg.norm(np.cross(lattice[(i + 1) % 3], lattice[(i + 2) % 3]))
                       for i in range(3)]
            span = np.ptp(cart @ np.linalg.inv(lattice), axis=0) if na > 1 else np.zeros(3)
            reach = [int(np.ceil((cutoff + s * d) / d)) for s, d in zip(span, spacing)]
        shifts = [np.array([i, j, k]) @ lattice
                  for i in range(-reach[0], reach[0] + 1)
                  for j in range(-reach[1], reach[1] + 1)
                  for k in range(-reach[2], reach[2] + 1)]
    # 原子 i ごとに全 j・周期像を配列で処理する。全 i,j,像 の配列は作らず、
    # 作業メモリを抑えながら Python の三重ループを避ける。
    shifts = np.asarray(shifts)
    for i in range(na):
        rij = cart[None, :, :] + shifts[:, None, :] - cart[i]
        distance = np.linalg.norm(rij, axis=2)
        valid = distance >= 1e-6
        if cutoff is not None:
            valid &= distance <= cutoff + 1e-12 * max(1.0, cutoff)
        weight = np.zeros_like(distance)
        weight[valid] = k0 / distance[valid] ** power
        if isotropic:
            blocks = weight.sum(axis=0)[:, None, None] * np.eye(3)
        else:
            unit = np.zeros_like(rij)
            np.divide(rij, distance[:, :, None], out=unit, where=valid[:, :, None])
            blocks = np.einsum("sj,sja,sjb->jab", weight, unit, unit)
        phi[3*i:3*i+3, :] -= blocks.transpose(1, 0, 2).reshape(3, 3 * na)
        phi[3*i:3*i+3, 3*i:3*i+3] += blocks.sum(axis=0)
    return 0.5 * (phi + phi.T)


def invariant_tiebreaker(sym_mats, seed=12345):
    """対称操作 D(g) のすべてと可換な、一般の対称行列 H を作る。

    H = (1/h) Σ_g D(g)ᵀ A D(g) (A は乱数の対称行列)。既約表現の縮退は保つ。
    H の固有空間は対称操作で
    閉じるので、網目モデルの固有値が縮退して向きが決まらないときに微小量を
    加えると、一般の場合に各モードを既約表現の 1 組 (d 本) ずつに分けられる。"""
    n = sym_mats[0].shape[0]
    a = np.random.default_rng(seed).standard_normal((n, n))
    a = a + a.T
    h = sum(D.T @ a @ D for D in sym_mats) / len(sym_mats)
    return 0.5 * (h + h.T)


def refine_in_subspace(basis, cart, masses, lattice=None, cutoff=None, mass_weighted=True,
                       return_eigenvalues=False, sym_mats=None, hessians=None, tiebreaker=None):
    """対称性適応部分空間の正規直交基底 basis (k, 3N) の中で、網目モデルの
    固有ベクトルを求めて返す (k, 3N)。低振動数順。

    対称種ごとの全コピーを含む部分空間に網目 Φ を射影して対角化し、
    モデルが予測する固有ベクトルを得る (任意の 1 コピーが Φ 不変とは限らない)。
    任意性が残るのは部分空間 *内* の向きだけで、部分空間そのものは厳密。

    mass_weighted=True : basis は質量重み座標 (分子)。標準固有値問題。
    mass_weighted=False: basis は Cartesian 座標 (結晶)。質量行列を含む
                         一般化固有値問題 K v = ω² M v を解く。
    return_eigenvalues=True なら (固有値 ω² (任意単位), ベクトル) を返す。
    sym_mats (対称操作の行列) を渡すと、網目の固有値が縮退する場合でも各モードが
    既約表現の 1 組になるよう、対称性を保つ微小項で縮退を解く (固有値は変えない)。

    網目の固有値が縮退したとき (平面分子の面外モードでは剛性がゼロ) の向きは、
    まず等方バネ網目 (幾何だけで決まり、原子の並び順によらない) の微小項で決め、
    それでも決まらなければ sym_mats の一般行列で決める。
    hessians=(中心力, 等方バネ), tiebreaker を渡すと、同じ構造の全対称種で
    計算済みの行列を再利用する (いずれも Cartesian 座標、関数内では変更しない)。"""
    mass3 = np.repeat(np.asarray(masses, float), 3)
    subreps = ([basis @ D @ basis.T for D in sym_mats]
               if sym_mats is not None else None)

    def projected(phi):
        if mass_weighted:
            phi = phi / np.sqrt(np.outer(mass3, mass3))
        k = basis @ phi @ basis.T
        if subreps is not None:
            # 数値丸めやカットオフ境界でも、モデルの群対称性を保つようにする。
            k = sum(D.T @ k @ D for D in subreps) / len(subreps)
        return 0.5 * (k + k.T)

    if hessians is None:
        hessians = (network_hessian(cart, lattice, cutoff=cutoff),
                    network_hessian(cart, lattice, cutoff=cutoff, isotropic=True))
    k, k_iso = (projected(phi) for phi in hessians)
    scale = max(np.abs(k).max(), 1e-3 * np.abs(k_iso).max(), 1e-12)
    k_solve = k + 1e-4 * scale * k_iso / max(np.abs(k_iso).max(), 1e-12)
    if sym_mats is not None:
        if tiebreaker is None:
            tiebreaker = invariant_tiebreaker(sym_mats)
        h = basis @ tiebreaker @ basis.T
        k_solve = k_solve + 1e-7 * scale * h / max(np.abs(h).max(), 1e-12)
    if mass_weighted:
        _, c = np.linalg.eigh(0.5 * (k_solve + k_solve.T))
        vecs = c.T @ basis
        lam = np.einsum("ip,ij,jp->p", c, k, c)            # 網目モデル本来の固有値
    else:
        m = basis @ (mass3[:, None] * basis.T)
        _, c = generalized_eigh(k_solve, m)
        lam = np.einsum("ip,ij,jp->p", c, k, c)
        vecs = c.T @ basis
        vecs = vecs / np.linalg.norm(vecs, axis=1, keepdims=True)
    return (lam, vecs) if return_eigenvalues else vecs


def generalized_eigh(k, m):
    """対称行列 K, 正定値行列 M の一般化固有値問題 K v = λ M v を解く。

    返り値 (λ 昇順, V)。V の列が固有ベクトル (M 直交)。"""
    k = 0.5 * (k + k.T)
    m = 0.5 * (m + m.T)
    lc = np.linalg.cholesky(m)
    li = np.linalg.inv(lc)
    lam, y = np.linalg.eigh(0.5 * ((li @ k @ li.T) + (li @ k @ li.T).T))
    return lam, li.T @ y


# ════════════════════════════════════════════════════════════════════════
# 幾何の小道具
# ════════════════════════════════════════════════════════════════════════
def rotation_to_z(axis):
    """単位ベクトル axis を +z に重ねる回転行列 R (R @ axis = ez) を返す。"""
    axis = np.asarray(axis, float) / np.linalg.norm(axis)
    z = np.array([0.0, 0.0, 1.0])
    v = np.cross(axis, z)
    s = np.linalg.norm(v)
    c = float(np.dot(axis, z))
    if s < 1e-8:                                     # すでに ±z 方向
        return np.eye(3) if c > 0 else np.diag([1.0, -1.0, -1.0])
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * ((1 - c) / s ** 2)


# ════════════════════════════════════════════════════════════════════════
# XSF (VESTA 用) 出力
# ════════════════════════════════════════════════════════════════════════
_LABEL_ASCII = str.maketrans({"⁺": "p", "⁻": "m", "Σ": "Sigma", "Π": "Pi",
                              "Δ": "Delta", "'": "p", '"': "pp"})


def safe_label(label: str) -> str:
    """既約表現ラベルをファイル名に使える ASCII に変換する (E' → Ep, Σg⁺ → Sigmagp)。"""
    return label.translate(_LABEL_ASCII)


def write_xsf(path, lattice, cart, numbers, disp, header=(), scale=1.0):
    """1 つの変位パターンを XSF で書き出す。

    lattice : 格子ベクトル (行ベクトル); 分子は大きな箱を渡す。
    disp    : (N, 3) 変位。scale 倍して矢印として書く。
    header  : 先頭に '# ' 付きで書くコメント行。
    Windows の既定文字コード (cp932) では Σ⁺ などを書けないので UTF-8 を明示する。"""
    with open(path, "w", encoding="utf-8") as f:
        for h in header:
            f.write(f"# {h}\n")
        f.write("CRYSTAL\nPRIMVEC\n")
        for v in lattice:
            f.write(f"  {v[0]:12.6f} {v[1]:12.6f} {v[2]:12.6f}\n")
        f.write(f"PRIMCOORD\n  {len(numbers)}  1\n")
        for z, r, d in zip(numbers, cart, np.asarray(disp) * scale):
            f.write(f"  {display_z(z):3d} {r[0]:12.6f} {r[1]:12.6f} {r[2]:12.6f}"
                    f"  {d[0]:10.6f} {d[1]:10.6f} {d[2]:10.6f}\n")


@lru_cache(maxsize=None)
def covalent_radius(z) -> float:
    """結合判定・表示用の原子半径 [Å] (pymatgen の計算値; 取れなければ 0.8)。"""
    z = display_z(z)
    if z <= 0:
        return 0.8
    from pymatgen.core import Element
    el = Element.from_Z(z)
    r = el.atomic_radius_calculated or el.atomic_radius
    return float(r) if r else 0.8


_VESTA_CUTOFFS = None


def _vesta_cutoffs():
    """VESTA の既定の結合探索表 {(元素1, 元素2): 最大結合長 [Å]}。

    pymatgen (MIT) が同梱している VESTA 既定値 (vesta_cutoffs.yaml) を読む。
    同じ組が両方の順で載っている場合 (H–O 1.2 と O–H 2.1 など) は短い方 (共有結合) を採る。"""
    global _VESTA_CUTOFFS
    if _VESTA_CUTOFFS is None:
        import re
        import pymatgen.core
        table = {}
        path = os.path.join(os.path.dirname(pymatgen.core.__file__), "vesta_cutoffs.yaml")
        try:
            with open(path, encoding="utf-8") as f:
                for line in f:
                    m = re.match(r"\s*\[\s*(\w+)\s*,\s*(\w+)\s*\]\s*:\s*([0-9.]+)", line)
                    if m:
                        key = tuple(sorted((m.group(1), m.group(2))))
                        table[key] = min(table.get(key, np.inf), float(m.group(3)))
        except OSError:
            pass                          # 表がなければ共有結合半径の規則だけで判定する
        _VESTA_CUTOFFS = table
    return _VESTA_CUTOFFS


def _bond_radius(z):
    """共有結合半径 [Å] (Cordero ら; 表にない元素・仮想原子は 0.8)。"""
    z = display_z(z)
    if z <= 0 or z >= PSEUDO_BASE:
        return 0.8
    from pymatgen.analysis.molecule_structure_comparator import CovalentRadius
    return float(CovalentRadius.radius.get(element_symbol(z), 0.8))


def bond_cutoff(z1, z2, molecule=False):
    """原子番号 z1, z2 の組を結合とみなす最大距離 [Å] (0 なら結合なし)。

    VESTA の既定の結合探索と同じく、元素の組ごとの表 (_vesta_cutoffs) で決める。
    表にない組は結合を引かない (例: Al–Al, Ti–Ti のような金属どうし)。ただし
      - 分子 (molecule=True) では、表にない組も共有結合半径の和 × 1.25 以内なら結合
        (H–H, P–Cl など表に載っていない分子内結合のため)
      - 表に元素そのものが載っていない (希ガス・仮想原子など) 場合も同じ規則で補う"""
    e1, e2 = (element_symbol(display_z(z)) if 0 < display_z(z) < PSEUDO_BASE else ""
              for z in (z1, z2))
    table = _vesta_cutoffs()
    key = tuple(sorted((e1, e2)))
    if key in table:
        return table[key]
    known = {e for pair in table for e in pair}
    if molecule or not table or e1 not in known or e2 not in known:
        return 1.25 * (_bond_radius(z1) + _bond_radius(z2))
    return 0.0


def bond_pairs(cart, numbers, molecule=False):
    """結合とみなす原子対 (i<j) のリスト (判定は bond_cutoff)。"""
    cart = np.asarray(cart, float)
    numbers = list(numbers)
    kinds = sorted(set(numbers))
    cut = {(a, b): bond_cutoff(a, b, molecule) for a in kinds for b in kinds}
    limit = np.array([[cut[(a, b)] for b in numbers] for a in numbers])
    d = np.linalg.norm(cart[:, None] - cart[None], axis=2)
    i, j = np.where(np.triu((d < limit) & (d > 1e-6), 1))
    return list(zip(i.tolist(), j.tolist()))


def _vesta_element(z):
    """VESTA の STRUC 欄に書く元素記号 (仮想原子は代表元素、不明なら XX)。"""
    dz = display_z(z)
    return element_symbol(dz) if 0 < dz < PSEUDO_BASE else "XX"


def write_vesta(path, cell, frac, numbers, labels, disp, title="",
                arrow_length=1.0, arrow_radius=0.15, arrow_rgb=(128, 128, 128),
                molecule=False):
    """1 つの変位パターンを VESTA 形式 (.vesta) で書き出す。

    cell   : 単位胞 (行ベクトル, Å);  frac: 分率座標 (N, 3);  title: 1 行の説明
    disp   : Cartesian 変位 (N, 3)。最大の矢印が中心から ±arrow_length [Å] になるよう規格化する。
    矢印は灰色・細め (半径 arrow_radius) で、矢印の中点を原子の中心に置く (VECTT の mode=1
    「原子を貫通」: 中心から ±ベクトル の範囲に描かれる。VESTA で開いて確認した仕様)。
    VESTA のベクトル成分は「a, b, c 軸方向の単位ベクトル」に沿った Å 単位の値なので、
    Cartesian 変位から変換して書く (VESTA で実際に開いて確認した仕様)。"""
    from pymatgen.core import Lattice
    cell = np.asarray(cell, float)
    frac = np.asarray(frac, float)
    disp = normalize_max(disp) * arrow_length
    axes = cell / np.linalg.norm(cell, axis=1, keepdims=True)      # a, b, c 方向の単位ベクトル
    comps = disp @ np.linalg.inv(axes)                            # v = u·â + v·b̂ + w·ĉ
    lat = Lattice(cell)
    title = safe_label(title or "mode").encode("ascii", "replace").decode()   # VESTA は ASCII で
    L = ["#VESTA_FORMAT_VERSION 3.5.0", "", "", "CRYSTAL", "", "TITLE", title, "",
          "GROUP", "1 1 P 1", "SYMOP",
          " 0.000000  0.000000  0.000000  1  0  0   0  1  0   0  0  1   1",
          " -1.0 -1.0 -1.0  0 0 0  0 0 0  0 0 0", "CELLP",
          "  {:.6f}  {:.6f}  {:.6f}  {:.6f}  {:.6f}  {:.6f}".format(*lat.abc, *lat.angles),
          "  0.000000   0.000000   0.000000   0.000000   0.000000   0.000000", "STRUC"]
    for i, (z, lbl, f) in enumerate(zip(numbers, labels, frac), 1):
        L.append(f"  {i} {_vesta_element(z):>2s} {lbl:>10s}  1.0000  "
                 f"{f[0]:.6f}  {f[1]:.6f}  {f[2]:.6f}    1a       1")
        L.append("                            0.000000   0.000000   0.000000  0.00")
    L += ["  0 0 0 0 0 0 0", "BOUND", "       0        1         0        1         0        1",
          "  0   0   0   0  0", "SBOND"]
    # 結合: 構造中の元素の組ごとに、GUI と同じ規則 (bond_cutoff) の最大結合長を書く
    pair_max = {}
    kinds = sorted(set(int(z) for z in numbers))
    for x in kinds:
        for y in kinds:
            dmax = bond_cutoff(x, y, molecule)
            if dmax > 0:
                a, b = sorted((_vesta_element(x), _vesta_element(y)))
                pair_max[(a, b)] = max(pair_max.get((a, b), 0.0), dmax)
    for k, ((a, b), dmax) in enumerate(sorted(pair_max.items()), 1):
        L.append(f"  {k}  {a:>4s}  {b:>4s}    0.00000    {dmax:.5f}  0  1  1  0  1  "
                 "0.120  1.000 160 160 160")
    L += ["  0 0 0 0", "VECTR"]
    for i, c in enumerate(comps, 1):
        L.append(f"   {i}  {c[0]:10.5f} {c[1]:10.5f} {c[2]:10.5f} 0")
        L.append(f"    {i}   0    0    0    0")
        L.append(" 0 0 0 0 0")
    L += [" 0 0 0 0 0", "VECTT"]
    for i in range(1, len(comps) + 1):
        L.append(f"   {i}  {arrow_radius:.3f} {arrow_rgb[0]:3d} {arrow_rgb[1]:3d} "
                 f"{arrow_rgb[2]:3d} 1")
    L += [" 0 0 0 0 0", "STYLE", "VECTS 1.000000", ""]
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(L))


def axis_aligned_components(disps, axes, sort_keys=None):
    """縮退モード (d 成分) の表示用の成分を、軸の方向を基準に選び直す。

    縮退した d 個の成分は、その d 次元空間内でどう取り直しても同じく正しい
    (対称性が決めるのは空間だけで、成分の向きは任意)。表示を物質によらず同じ規則で
    決めるため、次の手順で選ぶ。物質ごとの設定や特別扱いはない。

      1. 軸 1, 2, 3 (結晶なら a, b, c、分子なら x, y, z) の順に、縮退空間の中で
         「全原子の変位のうち、その軸方向の成分の割合」が最大になるパターンを採る
         (d×d 行列の最大固有ベクトル)。全原子の和なので原子の並び順によらない。
         - 既に採った成分とほぼ同じ向き (30° 以内) なら、それと直交する残りの中で選び直す
         - その軸方向にまったく動かない、または最大が一意に決まらない軸は飛ばす
      2. それでも d 個そろわないとき (例: 全原子が c 方向にだけ動く E モード) は、
         原子を座標順 (sort_keys の辞書順) に並べ、
         「その原子がその軸方向に動く」パターンの射影で残りを決める。
         並べ方は座標で決まるので、入力ファイルの原子の順番にはよらない。
      3. 各成分の符号は「軸方向の変位の総和が正」(総和 0 なら座標順で最初の
         動く原子の最大成分が正) にそろえる。

    手順 1 で決まった成分 k は「軸 k 方向に最もよく沿った動き」で、黒鉛の E モード
    や立方晶の T モードのように軸に沿って動ける場合は軸と完全に平行になる。
    成分どうしは直交するとは限らない (六方晶の a, b は 120°)。

    disps: d 個の (N, 3) 変位;  axes: 3 本の軸の方向 (行ベクトル)
    sort_keys: (N, 3) 原子を並べる座標 (結晶は分率座標、分子は Cartesian; 省略時は入力順)
    返り値: d 個の (N, 3) 変位 (最大変位 1 に規格化)。空間は元と同じ。
    d = 1 (縮退なし) でも符号だけ手順 3 でそろえる (出力ファイルが毎回同じになるように)。"""
    d = len(disps)
    if d == 0:
        return []
    n = len(disps[0])
    span = np.array([np.asarray(u, float).reshape(-1) for u in disps])
    q, _ = np.linalg.qr(span.T)
    q = q[:, :d]                                        # 縮退空間の正規直交基底 (3N, d)
    axes = np.asarray(axes, float)
    axes = axes / np.linalg.norm(axes, axis=1, keepdims=True)
    if sort_keys is None:
        order = list(range(n))
    else:                                               # 座標の辞書順 (入力の並び順によらない)
        keys = np.round(np.asarray(sort_keys, float), 3) + 0.0     # -0.0 → 0.0
        order = sorted(range(n), key=lambda i: tuple(keys[i]))

    chosen = []                                         # 係数ベクトル (d,) のリスト
    chosen_axis = []                                    # 各成分を決めた軸の番号

    def complement():
        """縮退空間のうち、採った成分と直交する部分の基底 (d, k)。"""
        if not chosen:
            return np.eye(d)
        _, s, vt = np.linalg.svd(np.array(chosen))
        return vt[len(chosen):].T

    def dependent(c):
        if not chosen:
            return False
        b = np.linalg.qr(np.array(chosen).T)[0]
        resid = c - b @ (b.T @ c)
        return np.linalg.norm(resid) < 0.5 * np.linalg.norm(c)     # 30° 以内

    # 手順 1: 軸方向の割合が最大のパターン
    for k, ax in enumerate(axes):
        if len(chosen) == d:
            break
        along = np.einsum("nid,i->nd", q.reshape(n, 3, d), ax)    # 各原子の軸方向成分
        m_full = along.T @ along                                   # d×d (割合の行列)
        for basis in (np.eye(d), None):
            if basis is None:                                      # 既採用と重なったら補空間で
                basis = complement()
            w, v = np.linalg.eigh(basis.T @ m_full @ basis)
            w, v = w[::-1], v[:, ::-1]
            unique = len(w) == 1 or w[0] - w[1] > 1e-4
            c = basis @ v[:, 0]
            if w[0] < 1e-4 or not unique:
                break                                              # この軸では決まらない
            if not dependent(c):
                chosen.append(c / np.linalg.norm(c))
                chosen_axis.append(k)
                break

    # 手順 2: 座標順に並べた原子の、軸方向の単独変位を射影
    for i in order:
        for k, ax in enumerate(axes):
            if len(chosen) == d:
                break
            ref = np.zeros((n, 3))
            ref[i] = ax
            c = q.T @ ref.reshape(-1)
            if np.linalg.norm(c) < 1e-3:
                continue
            basis = complement()
            c_perp = basis @ (basis.T @ c)
            if np.linalg.norm(c_perp) > 1e-3:
                chosen.append(c / np.linalg.norm(c) if not dependent(c)
                              else c_perp / np.linalg.norm(c_perp))
                chosen_axis.append(k)
    if len(chosen) < d:                                            # (念のため) 元のまま
        return [normalize_max(u) for u in disps]

    # 手順 3: 符号をそろえる (成分を決めた軸 k の方向を基準に)
    result = []
    for c, k in zip(chosen, chosen_axis):
        u = (q @ c).reshape(n, 3)
        along = u @ axes[k]                                        # 各原子の軸 k 方向の変位
        big = np.abs(along).max()
        s = along.sum()
        if abs(s) < 1e-6 * big * n:                                # 正味ゼロ → 座標順で最初の原子
            s = next(along[i] for i in order if abs(along[i]) > 0.1 * big)
        result.append(normalize_max(u if s >= 0 else -u))
    return result


def normalize_max(disp):
    """変位 (N, 3) を最大原子変位が 1 になるよう規格化する (表示用)。"""
    disp = np.asarray(disp, float)
    mx = np.max(np.linalg.norm(disp, axis=1)) if len(disp) else 0.0
    return disp / mx if mx > 1e-10 else disp


def mode_file_stem(index: int, mode: Mode, component: int) -> str:
    """XSF ファイル名 (拡張子なし): mode_03_E2g_a_Raman。"""
    suffix = f"_{chr(97 + component)}" if len(mode.displacements) > 1 else ""
    return f"mode_{index:02d}_{safe_label(mode.label)}{suffix}_{mode.activity}"


def xsf_header(index: int, mode: Mode, component: int, exact_note: str) -> list[str]:
    """XSF 先頭のコメント行 (分子・結晶共通)。"""
    lines = [f"Mode {index}: {mode.label}  Activity: {mode.activity}"]
    if mode.character:
        lines.append(f"Displacement: {mode.character_label}")
    if mode.refined:
        lines.append(f"Direction: NETWORK-MODEL GUESS (m>1; 1/d^2 springs) - {exact_note}")
    else:
        lines.append("Direction: symmetry-determined (exact)")
    if len(mode.displacements) > 1:
        lines.append(f"Degenerate component {component + 1}/{len(mode.displacements)}")
    if mode.tensor:
        lines.append(f"Raman tensor: {', '.join(mode.tensor)}")
    if mode.contributions:
        lines.append("Sites: " + ", ".join(mode.contributions))
    lines.append("VESTA: Properties -> Vectors -> Show vectors")
    return lines


def ensure_outdir(outdir):
    os.makedirs(outdir, exist_ok=True)
    return outdir
