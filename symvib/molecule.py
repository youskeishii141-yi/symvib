"""孤立分子の振動モード解析 (DFT 不要)。

CIF/xyz を入れると分子点群を判定し、内部振動モードを既約表現 (Mulliken 記号)
で分類して、各モードの Raman/IR 活性と「どの原子がどの方向に動くか」を返す。

────────────────────────────────────────────────────────────────────────
処理の流れ (analyze_molecule)
────────────────────────────────────────────────────────────────────────
 1. 重心を原点に移し、許容差 tol 以内で分子を自分自身に写す直交変換を
    すべて探す (find_point_group_operations)。
 2. その操作群で平均して構造を厳密に対称化し、pymatgen で点群記号を得る。
    操作の数が指標表の位数と一致することを確認する (不一致なら中断)。
 3. 点群の種類に応じて既約分解する:
    - 直線分子 (C∞v, D∞h): 回転角 φ の連続族を数値積分した真の連続群射影
    - 複素既約表現を持つ群 (Cn, Cnh, Sn, T, Th): 複素指標による射影
    - 正二十面体群 (I, Ih): 回転角だけでクラスが決まるので角度で照合
    - その他: 各操作を幾何学的に共役類へ割り当てて射影
      (割り当てが一意でない立方群などは、各モードの指標照合に自動切替)
 4. 同じ既約表現が複数回 (m>1) 現れる場合は 1/d² バネ網目モデルで
    部分空間内の向きを推測する (m=1 は対称性だけで厳密に決まる)。

────────────────────────────────────────────────────────────────────────
座標軸の取り方 (C2v, D2, D2h)
────────────────────────────────────────────────────────────────────────
B1/B2 (C2v) や B1/B2/B3 (D2, D2h) のラベルは軸の取り方に依存する。本モジュールは
Mulliken の推奨に従い、平面分子では分子面の法線を x にとる。C2v では z を C2 軸、
D2/D2h では z を最も多くの原子を通る C2 軸とする。採用した軸はレポートに明記する。
"""
from __future__ import annotations

import math
import os
from collections import Counter

import numpy as np

from .common import (AnalysisError, AnalysisResult, Mode, ACTIVITIES, ACTIVITY_SYMBOL,
                     axis_aligned_components,
                     atomic_mass, element_symbol, ensure_outdir, mode_file_stem,
                     invariant_tiebreaker, network_hessian, normalize_max,
                     refine_in_subspace, rotation_to_z, site_element,
                     write_vesta, write_xsf, xsf_header)
from .i18n import Msg, pad, t
from .molecule_tables import CHAR_TABLES, COMPLEX_TABLES, ICOSAHEDRAL_CLASS_KEYS

# 原子の対応づけ (対称操作で写った先の原子を探す) に使う許容差 [Å]。
# 構造は事前に厳密に対称化するので小さくてよい。
MATCH_TOL = 1e-2

# 表示ラベル: 指標表の ASCII 表記 → 文献表記
_PRIME_DISPLAY = (("pp", "''"), ("p", "'"))
_LINEAR_DISPLAY = {"Sigma_g+": "Σg⁺", "Sigma_u+": "Σu⁺", "Pi_g": "Πg", "Pi_u": "Πu",
                   "Sigma+": "Σ⁺", "Pi": "Π"}
# 直線分子の活性と Raman テンソル (標準)
_LINEAR_ACTIVITY = {
    "Sigma_g+": ("Raman", ["xx+yy", "zz"]), "Sigma_u+": ("IR", []),
    "Pi_g": ("Raman", ["xz", "yz"]), "Pi_u": ("IR", []),
    "Sigma+": ("Raman+IR", ["xx+yy", "zz"]), "Pi": ("Raman+IR", ["xz", "yz"]),
}


def display_label(name: str) -> str:
    """指標表の既約表現名を表示用に変換する (A1p → A1', E2pp → E2'')。"""
    for ascii_, prime in _PRIME_DISPLAY:
        if name.endswith(ascii_):
            return name[: -len(ascii_)] + prime
    return name


# ════════════════════════════════════════════════════════════════════════
# 1. 構造の読み込み
# ════════════════════════════════════════════════════════════════════════
def load_molecule(path, notes=None):
    """xyz / CIF / その他 pymatgen が読める形式から (cart[Å], Z, symbols) を返す。

    CIF の場合は単位胞内の原子をそのまま 1 分子として扱う
    (分子が胞の境界で切れていない CIF を入力すること)。
    notes (リスト) を渡すと、読み込み時の注記 (部分占有をまとめた等) を書き足す。"""
    from pymatgen.core import Molecule, Structure
    ext = os.path.splitext(path)[1].lower()
    if ext == ".xyz":
        return read_xyz(path)
    if ext == ".cif":
        from .cif_reader import read_cif
        c = read_cif(path)
        if notes is not None:
            notes.extend(c.notes)
        return c.frac @ c.lattice, list(c.numbers), list(c.names)
    else:
        try:
            obj = Molecule.from_file(path)
        except Exception as e_mol:
            try:
                obj = Structure.from_file(path)
            except Exception as e_st:
                raise AnalysisError(Msg("err.read_struct", path=path, e1=e_mol, e2=e_st)) from e_st
    elems = [site_element(site) for site in obj]
    return np.array(obj.cart_coords, float), [z for z, _ in elems], [s for _, s in elems]


def read_xyz(path):
    """xyz ファイル (1 行目: 原子数, 2 行目: コメント, 以降: 元素 x y z [Å]) を読む。

    複数構造が連なった xyz では最初の構造だけを使う。元素は記号 ("C", "C1") でも
    原子番号でもよい。元素を特定できない記号は質量 1 amu の仮想原子として扱う。"""
    from pymatgen.core import Element
    from .cif_reader import _element_from
    from .common import register_pseudo_species
    with open(path, encoding="utf-8-sig", errors="replace") as f:   # BOM 付きも読む
        lines = f.read().splitlines()
    try:
        n = int(lines[0].split()[0])
        rows = [ln.split() for ln in lines[2:2 + n]]
        if len(rows) < n:
            raise ValueError(t("err.xyz_rows", n=n, m=len(rows)))
        numbers, symbols, cart = [], [], []
        for row in rows:
            tok = row[0]
            if tok.isdigit():
                el = Element.from_Z(int(tok))
                numbers.append(el.Z)
                symbols.append(el.symbol)
            else:
                sym, mass = _element_from(tok)
                if sym is None:                             # 未知の記号 → 仮想原子
                    numbers.append(register_pseudo_species(tok, 1.0, 0, key=("xyz", tok)))
                    symbols.append(tok)
                elif sym == "H" and mass > 1.5:            # 重水素 D・三重水素 T
                    numbers.append(register_pseudo_species(tok[:1].upper(), mass, 1,
                                                           key=("iso", mass)))
                    symbols.append(tok[:1].upper())
                else:
                    numbers.append(Element(sym).Z)
                    symbols.append(sym)
            if len(row) < 4:
                raise ValueError(t("err.xyz_rows", n=n, m=len(cart)))
            cart.append([float(x) for x in row[1:4]])
    except (IndexError, ValueError, KeyError) as e:
        raise AnalysisError(Msg("err.xyz_format", path=path, err=e)) from e
    return np.array(cart, float), numbers, symbols


def water_demo():
    """デモ用 H2O (C2v)。"""
    return (np.array([[0, 0, 0.117], [0, 0.757, -0.467], [0, -0.757, -0.467]]),
            [8, 1, 1], ["O", "H", "H"])


# ════════════════════════════════════════════════════════════════════════
# 2. 対称操作の探索と構造の対称化
# ════════════════════════════════════════════════════════════════════════
def find_point_group_operations(cart, numbers, tol):
    """重心原点の座標 cart について、分子を自分自身に写す直交変換をすべて求める。

    非共線な 2 原子 a, b を基準に選び、同種・同距離の原子対 (a', b') への写像
    (a, b, a×b) → (a', b', ±a'×b') を候補として作り、全原子が許容差 tol 以内で
    同種原子に重なるものを採用する。
    返り値: [(R (3×3), perm)]。perm[i] は原子 i の写り先。

    pymatgen の PointGroupAnalyzer.get_symmetry_operations() は Sn 群や Th 群で
    操作を取りこぼす (例: Th で 24 個中 6 個) ため、操作は自前で求める。"""
    r = np.asarray(cart, float)
    Z = np.asarray(numbers)
    norms = np.linalg.norm(r, axis=1)
    off = [i for i in range(len(r)) if norms[i] > tol]
    if not off:                                       # 単原子
        return [(np.eye(3), np.arange(len(r)))]

    def n_equiv(i):
        return int(np.sum((Z == Z[i]) & (np.abs(norms - norms[i]) < tol)))

    def sin_angle(i, j):
        return np.linalg.norm(np.cross(r[i], r[j])) / (norms[i] * norms[j])

    # 基準の 2 原子は外側の原子から選ぶ (原点に近い原子は向きの誤差が大きく、歪んだ構造では
    # 2 原子から作る回転が外側の原子で tol を超えてしまう)。候補の少ない原子を優先する
    outer = [i for i in off if norms[i] >= max(3 * tol, 0.4 * norms.max())] or off
    a = min(outer, key=lambda i: (n_equiv(i), -norms[i]))
    cands_b = [j for j in outer if sin_angle(a, j) > 0.3] or [j for j in off if sin_angle(a, j) > 0.1]
    if not cands_b:
        raise AnalysisError(Msg("err.linear_internal"))
    b = min(cands_b, key=lambda j: (n_equiv(j), -sin_angle(a, j) * norms[j]))

    same_a = [i for i in off if Z[i] == Z[a] and abs(norms[i] - norms[a]) < tol]
    same_b = [j for j in off if Z[j] == Z[b] and abs(norms[j] - norms[b]) < tol]
    d_ab = np.linalg.norm(r[a] - r[b])
    f_inv = np.linalg.inv(np.column_stack([r[a], r[b], np.cross(r[a], r[b])]))
    same_z = Z[:, None] == Z[None, :]

    found = {}
    for i in same_a:
        for j in same_b:
            if i == j or abs(np.linalg.norm(r[i] - r[j]) - d_ab) > 2 * tol:
                continue
            for s in (1.0, -1.0):                     # 真回転 / 回映・鏡映
                R = _nearest_orthogonal(
                    np.column_stack([r[i], r[j], s * np.cross(r[i], r[j])]) @ f_inv, s)
                perm = _match_atoms(r, r @ R.T, same_z, 3 * tol)
                if perm is None:
                    continue
                # 2 原子だけで作った R は誤差を拡大しうるので、全原子で最小二乗 (Kabsch) し直す
                R = _nearest_orthogonal(r[perm].T @ r, s)
                perm = _match_atoms(r, r @ R.T, same_z, tol)
                if perm is not None:
                    # 原子の置換だけでは操作は決まらない (平面分子では分子面の鏡映が
                    # 恒等操作と同じ置換になる)。置換 + det の符号なら一意。
                    found.setdefault((tuple(perm.tolist()), s), (R, perm))
    return list(found.values())


def _nearest_orthogonal(M, det_sign):
    """M に最も近い、行列式の符号が det_sign の直交行列 (特異値分解による)。"""
    u, _, vt = np.linalg.svd(M)
    if np.linalg.det(u @ vt) * det_sign < 0:
        u[:, -1] *= -1
    return u @ vt


def _match_atoms(ref, moved, same_z, tol):
    """moved[i] と tol 以内で重なる同種原子 perm[i] を返す (全単射でなければ None)。"""
    dist = np.linalg.norm(moved[:, None, :] - ref[None, :, :], axis=2)
    dist[~same_z] = np.inf
    perm = np.argmin(dist, axis=1)
    if np.all(dist[np.arange(len(ref)), perm] < tol) and len(set(perm.tolist())) == len(ref):
        return perm
    return None


def symmetrize_by_group(cart, ops):
    """操作群で平均して厳密に対称な座標を作る: r_j = (1/h) Σ_g g·r_{g⁻¹(j)}。"""
    acc = np.zeros_like(cart)
    for R, perm in ops:
        acc[perm] += cart @ R.T
    return acc / len(ops)


def _symmetrize_linear(cart, numbers, tol):
    """直線分子を分子軸上に射影し、反転中心があれば反転について平均する。

    返り値: (対称化した座標, 反転中心の有無)。"""
    axis = _inertia_axes(cart, np.ones(len(cart)))[:, 0]   # 慣性最小軸 = 分子軸
    t = cart @ axis
    sym = np.outer(t, axis)
    Z = np.asarray(numbers)
    partner = []
    for i in range(len(t)):
        cand = [j for j in range(len(t)) if Z[j] == Z[i] and abs(t[j] + t[i]) < tol]
        partner.append(cand[0] if cand else None)
    has_inversion = all(p is not None for p in partner)
    if has_inversion:
        t = 0.5 * (t - t[partner])
        sym = np.outer(t, axis)
    return sym, has_inversion


def _inertia_axes(cart, masses):
    """慣性主軸 (列ベクトル, 慣性モーメント昇順)。"""
    m = np.asarray(masses, float)
    r = cart - (m[:, None] * cart).sum(0) / m.sum()
    inertia = sum(m[i] * (r[i] @ r[i] * np.eye(3) - np.outer(r[i], r[i]))
                  for i in range(len(r)))
    _, vecs = np.linalg.eigh(inertia)
    return vecs


def is_linear(cart, masses=None, tol=0.1):
    """直線分子か: すべての原子が、重心を通る分子軸から tol [Å] 以内にあるか。

    慣性モーメントの比ではなく距離で判定する (長い分子の末端の H のように、
    軸から 1 Å も離れた原子を「直線」とみなして軸上へ動かしてしまわないように)。"""
    cart = np.asarray(cart, float)
    na = len(cart)
    if na <= 2:
        return True
    m = np.ones(na) if masses is None else np.asarray(masses, float)
    r = cart - (m[:, None] * cart).sum(0) / m.sum()
    axis = _inertia_axes(cart, m)[:, 0]                 # 慣性最小軸 = 分子軸の候補
    dist = np.linalg.norm(r - np.outer(r @ axis, axis), axis=1)
    return bool(np.max(dist) < tol)


def _group_order(pg):
    table = CHAR_TABLES.get(pg) or COMPLEX_TABLES.get(pg)
    return int(round(sum(table["g"]))) if table else None


def _point_group_from_ops(ops):
    """操作の種類 (種別, 回転角) の個数から、収録している点群の記号を求める。

    pymatgen の判定が見つけた操作と合わないとき (細長い分子・ほぼ直線の分子など) の
    代わりに使う。一致する点群が 1 つに決まらなければ None。"""
    def fp_of(keys_counts):
        c = Counter()
        for (kind, ang), n in keys_counts:
            c[(kind, int(round(abs(ang))))] += int(round(n))
        return c

    found = fp_of((op_angle_key(R), 1) for R, _ in ops)
    hits = []
    for pg, T in CHAR_TABLES.items():
        if fp_of((_table_class_key(c), g) for c, g in zip(T["classes"], T["g"])) == found:
            hits.append(pg)
    for pg, T in COMPLEX_TABLES.items():
        if fp_of(zip(T["fine"], T["g"])) == found:
            hits.append(pg)
    return hits[0] if len(hits) == 1 else None


def mechanical_matrices(ops, na):
    """各操作の力学表現 D(R) (3N×3N): 原子 i の変位を R で回して原子 perm[i] へ。"""
    mats = []
    for R, perm in ops:
        D = np.zeros((3 * na, 3 * na))
        for i, k in enumerate(perm):
            D[3*k:3*k+3, 3*i:3*i+3] = R
        mats.append(D)
    return mats


# ════════════════════════════════════════════════════════════════════════
# 3. 操作の分類 (回転角と det から)
# ════════════════════════════════════════════════════════════════════════
def classify_op(R, tol=1e-2):
    """対称操作 R を (種別, 次数) に粗く分類: ('E',1), ('C',n), ('i',1), ('S',n), ('sigma',1)。

    次数は round(360/回転角) なので C5² (144°) と C2 (180°) の区別はできない。"""
    det = np.linalg.det(R)
    tr = np.trace(R)
    if det > 0:
        ang = math.degrees(math.acos(np.clip((tr - 1) / 2, -1, 1)))
        return ("E", 1) if ang < tol else ("C", round(360 / ang))
    ang = math.degrees(math.acos(np.clip((tr + 1) / 2, -1, 1)))
    if ang < tol:
        return "sigma", 1
    if abs(ang - 180) < tol:
        return "i", 1
    return "S", round(360 / ang)


def op_angle_key(R, tol=2.0):
    """操作 R を (種別, 回転角[度]) に分類する。C5² (144°) と C2 (180°) も区別できる。"""
    det = np.linalg.det(R)
    tr = np.trace(R)
    if det > 0:
        ang = math.degrees(math.acos(np.clip((tr - 1) / 2, -1, 1)))
        return ("E", 0.0) if ang < tol else ("C", round(ang))
    ang = math.degrees(math.acos(np.clip((tr + 1) / 2, -1, 1)))
    if ang < tol:
        return "sigma", 0.0
    if abs(ang - 180) < tol:
        return "i", 180.0
    return "S", round(ang)


def op_signed_key(R, paxis, tol=2.0):
    """操作 R を (種別, 主軸まわりの符号付き回転角[度]) に分類する。

    複素既約表現を持つ群では C3 と C3² のように回転の向きを区別する必要がある。
    回映は S(θ) = σh·C(θ) の θ で表す (σh: θ=0 → 'sigma', i: θ=180 → 'i')。
    −S(θ) = C(θ+180°) なので、真回転部分 M = −R の回転角 α から θ = α − 180° と求める。"""
    det = np.linalg.det(R)
    M = R if det > 0 else -R                     # 真回転部分
    w, v = np.linalg.eig(M)
    k = int(np.argmin(np.abs(w - 1)))            # 固有値 1 の固有ベクトル = 回転軸
    ax = np.real(v[:, k])
    ax = ax / (np.linalg.norm(ax) + 1e-30)
    if np.dot(ax, paxis) < 0:
        ax = -ax
    alpha = math.degrees(math.acos(np.clip((np.trace(M) - 1) / 2, -1, 1)))
    if alpha > tol:                              # 回転の向き (主軸に対する符号)
        perp = np.array([1., 0., 0.]) if abs(ax[0]) < 0.9 else np.array([0., 1., 0.])
        perp = perp - np.dot(perp, ax) * ax
        perp /= np.linalg.norm(perp) + 1e-30
        if np.dot(np.cross(perp, M @ perp), ax) < 0:
            alpha = -alpha
    if det > 0:
        return ("E", 0.0) if abs(alpha) < tol else ("C", float(round(alpha)))
    theta = (alpha - 180.0 + 180.0) % 360.0 - 180.0   # (−180, 180] に正規化
    if abs(theta) < tol:
        return "sigma", 0.0
    if abs(abs(theta) - 180.0) < tol:
        return "i", 180.0
    return "S", float(round(theta))


def _rotation_axis(R):
    """回転 (または回映) の軸方向。鏡映では面法線になる。"""
    M = R if np.linalg.det(R) > 0 else -R
    w, v = np.linalg.eig(M)
    for i in range(3):
        if abs(w[i] - 1) < 1e-3:
            ax = np.real(v[:, i])
            return ax / (np.linalg.norm(ax) + 1e-30)
    return np.array([0., 0., 1.])


def _mirror_normal(R):
    """鏡映の面法線 (固有値 −1 の固有ベクトル)。"""
    w, v = np.linalg.eig(R)
    for i in range(3):
        if abs(w[i] + 1) < 1e-3:
            ax = np.real(v[:, i])
            return ax / (np.linalg.norm(ax) + 1e-30)
    return np.array([1., 0., 0.])


def _principal_axis(ops):
    """最高次の真回転軸 (主軸) の次数と方向。"""
    maxn, paxis = 1, np.array([0., 0., 1.])
    for R, _ in ops:
        t, n = classify_op(R)
        if t == "C" and n > maxn:
            maxn, paxis = n, _rotation_axis(R)
    return maxn, paxis


def _count_on_line(direction, cart, tol=0.15):
    """原点を通る direction 方向の直線上にある原子の数 (距離 tol Å 以内)。"""
    n = direction / np.linalg.norm(direction)
    perp = cart - np.outer(cart @ n, n)
    return int(np.sum(np.linalg.norm(perp, axis=1) < tol))


def _count_in_plane(normal, cart, tol=0.15):
    """原点を通り法線 normal の面上にある原子の数 (距離 tol Å 以内)。"""
    n = normal / np.linalg.norm(normal)
    return int(np.sum(np.abs(cart @ n) < tol))


def _table_class_key(name):
    """指標表のクラス名を op_angle_key と同じ (種別, 符号なし回転角[度]) にする。

    'C5_2' = C5² (144°), 'S10_3' = S10³ (σh·C(108°)), 'C2p'/'C2x' = 2 回軸、
    's…' = 鏡映。回映 S(θ) の角は op_angle_key と同じく min(θ, 360°−θ)。"""
    import re
    if name == "E":
        return "E", 0.0
    if name == "i":
        return "i", 180.0
    if name.startswith("s"):
        return "sigma", 0.0
    m = re.match(r"([CS])(\d+)(?:_(\d+))?", name)
    theta = (360.0 * int(m.group(3) or 1) / int(m.group(2))) % 360.0
    return m.group(1), min(theta, 360.0 - theta)


def _group_structure(ops, tol=1e-3):
    """操作の積表 mult[a][b] = index(R_a R_b)、逆元 inv[a]、共役類 (番号のリスト) を返す。"""
    rs = [np.asarray(R, float) for R, _ in ops]
    flat = np.array([R.ravel() for R in rs])

    def index(M):
        d = np.abs(flat - M.ravel()).max(axis=1)
        k = int(np.argmin(d))
        if d[k] > tol:
            raise AnalysisError(Msg("err.not_group"))
        return k

    h = len(rs)
    mult = np.array([[index(a @ b) for b in rs] for a in rs])
    inv = np.array([index(R.T) for R in rs])
    seen, classes = set(), []
    for a in range(h):
        if a in seen:
            continue
        cl = sorted({int(mult[mult[g][a]][inv[g]]) for g in range(h)})
        seen.update(cl)
        classes.append(cl)
    return mult, inv, classes


def _valid_class_function(table, col, mult, inv):
    """列の割り当て col (操作 → 表の列) で、各既約表現の射影が正則表現上で
    本当に射影演算子 (P² = P, tr P = d²) になっているかを確かめる。

    クラスの割り当てや指標表の誤りはここで必ず検出される (直交性の検査では
    見逃す誤り、例えば同じ大きさのクラスの列の入れ替えも検出できる)。"""
    h = len(col)
    g_of = mult[:, inv]                                  # g_of[y, x] = y·x⁻¹
    for row in table:
        d = row[0]
        chi = np.array([row[c] for c in col], float)
        P = (d / h) * chi[g_of]
        if np.abs(P @ P - P).max() > 1e-6 or abs(np.trace(P) - d * d) > 1e-6:
            return False
    return True


# C2'/C2'' と σv/σd は大きさ・回転角が同じで、どちらと呼ぶかは慣例 (Mulliken):
# C2' と σv は原子を多く通る方。どちらでも群論的には正しい分解になる。
_PREFER_THROUGH_ATOMS = {"C2p": 1, "C2pp": -1, "sv": 1, "sd": -1}


def _assign_classes(pg, ops, cart):
    """各操作を指標表の列 (クラス) に割り当てる。返り値: 列番号のリスト (操作の順)。

    共役類を群の積表から求め、(種別, 回転角, 類の大きさ) で表のクラスと対応させる。
    同じ特徴のクラスが複数あれば (C2'/C2'', σv/σd)、指標表と矛盾しない割り当てのうち
    原子を通る軸・面を C2', σv とするものを選ぶ。どの割り当ても矛盾するなら例外。"""
    import itertools
    T = CHAR_TABLES[pg]
    table = np.array(T["table"], float)
    mult, inv, classes = _group_structure(ops)
    rs = [R for R, _ in ops]

    groups = {}                                          # 特徴 → ([表の列], [共役類])
    for ci, name in enumerate(T["classes"]):
        kind, ang = _table_class_key(name)
        groups.setdefault((kind, round(ang, 1), int(round(T["g"][ci]))), ([], []))[0].append(ci)
    for cl in classes:
        kind, ang = op_angle_key(rs[cl[0]])
        key = next((k for k in groups if k[0] == kind and k[2] == len(cl)
                    and abs(k[1] - ang) < 1.5), None)
        if key is None:
            raise AnalysisError(Msg("err.class_assign", pg=pg))
        groups[key][1].append(cl)
    if any(len(cols) != len(cls) for cols, cls in groups.values()):
        raise AnalysisError(Msg("err.class_assign", pg=pg))

    def through_atoms(cl):
        R = rs[cl[0]]
        if np.linalg.det(R) > 0:
            return _count_on_line(_rotation_axis(R), cart)
        return _count_in_plane(_mirror_normal(R), cart)

    keys = list(groups)
    has_c2p = any(T["classes"][c] in ("C2p", "C2pp") for k in keys for c in groups[k][0]
                  if len(groups[k][0]) > 1)
    best, best_score = None, None
    for perms in itertools.product(*(itertools.permutations(groups[k][1]) for k in keys)):
        col = [None] * len(ops)
        score = 0
        for k, perm in zip(keys, perms):
            for ci, cl in zip(groups[k][0], perm):
                for o in cl:
                    col[o] = ci
                name = T["classes"][ci]
                if len(groups[k][0]) > 1 and name in _PREFER_THROUGH_ATOMS:
                    if has_c2p and name not in ("C2p", "C2pp"):
                        continue                         # σ の向きは C2' との整合で決まる
                    score += _PREFER_THROUGH_ATOMS[name] * through_atoms(cl)
        if (best_score is None or score > best_score) and \
                _valid_class_function(table, col, mult, inv):
            best, best_score = col, score
    if best is None:
        raise AnalysisError(Msg("err.class_assign", pg=pg))
    return best


def _abelian_frame(pg, ops, cart, masses):
    """C2v / D2 / D2h の分子固定軸 (x, y, z) を決める。

    Mulliken の推奨: 平面分子では分子面の法線を x にとる。
    C2v は z = C2 軸。D2/D2h は z = 最も多くの原子を通る C2 軸。
    返り値: (frame 行列 (行が x, y, z 軸), 説明文)。"""
    c2_axes = []
    for R, _ in ops:
        if classify_op(R) == ("C", 2):
            ax = _rotation_axis(R)
            if not any(abs(abs(ax @ a) - 1) < 1e-3 for a in c2_axes):
                c2_axes.append(ax)
    # 分子面の法線 (平面分子のみ)
    centered = cart - cart.mean(0)
    _, s, vt = np.linalg.svd(centered)
    planar_normal = vt[2] if s[2] < 1e-3 * max(s[0], 1e-12) else None

    def near(v, cands):
        return next((c for c in cands if abs(abs(v @ c) - 1) < 1e-3), None)

    if pg == "C2v":
        z = c2_axes[0]
        mirror_normals = [_mirror_normal(R) for R, _ in ops if classify_op(R) == ("sigma", 1)]
        x = near(planar_normal, mirror_normals) if planar_normal is not None else None
        if x is None:
            # 非平面: 原子を多く含む鏡面を yz 面とする (その法線が x)
            # 同数なら面上の原子の質量の和が大きい方 (向きによらない決め方)
            x = max(mirror_normals, key=lambda nrm: (
                int(np.sum(np.abs(cart @ nrm) < 0.15)),
                round(float(np.sum(masses[np.abs(cart @ nrm) < 0.15])), 3)))
        note = Msg("axis.c2v_planar" if planar_normal is not None else "axis.c2v_other")
    else:                                              # D2, D2h: 互いに直交する 3 本の C2 軸
        x = near(planar_normal, c2_axes) if planar_normal is not None else None
        rest = [a for a in c2_axes if x is None or abs(a @ x) < 0.5]
        def on_axis(a):                    # 軸上の原子の数、同数なら質量の和 (向きによらない)
            near = np.linalg.norm(cart - np.outer(cart @ a, a), axis=1) < 0.15
            return int(near.sum()), round(float(np.sum(masses[near])), 3)

        z = max(rest, key=on_axis)
        if x is None:
            x = max([a for a in rest if abs(a @ z) < 0.5], key=on_axis)
        note = Msg("axis.d2_planar" if planar_normal is not None else "axis.d2_other")
    y = np.cross(z, x)
    return np.array([x, y, z]), note


def _abelian_class(R, frame, pg):
    """C2v / D2 / D2h の操作を、分子固定軸に対する向きでクラス名へ割り当てる。"""
    t, _ = classify_op(R)
    if t in ("E", "i"):
        return t
    names = "xyz"
    if t == "C":
        k = int(np.argmax(np.abs(frame @ _rotation_axis(R))))
        return "C2" if pg == "C2v" else f"C2{names[k]}"
    k = int(np.argmax(np.abs(frame @ _mirror_normal(R))))     # 法線の向き
    plane = "".join(c for i, c in enumerate(names) if i != k)  # 法線 x → yz 面
    return f"sv_{plane}" if pg == "C2v" else f"s{plane}"


# ════════════════════════════════════════════════════════════════════════
# 4. 振動部分空間と射影
# ════════════════════════════════════════════════════════════════════════
def translation_rotation_basis(cart, masses):
    """質量重み座標での並進・回転の正規直交基底 (T (3,3N), R (2or3,3N))。

    回転は重心まわりの δr_i = e×(r_i − R_com)。直線分子では軸まわりの回転が
    消えるので、ランク落ちから自動的に 2 本になる。"""
    na = len(cart)
    m = np.asarray(masses, float)
    sm = np.sqrt(m)
    r = cart - (m[:, None] * cart).sum(0) / m.sum()
    trans = np.zeros((3, 3 * na))
    for ax in range(3):
        trans[ax, ax::3] = sm
    trans /= np.linalg.norm(trans, axis=1, keepdims=True)
    rots = []
    for ax in range(3):
        e = np.zeros(3)
        e[ax] = 1.0
        v = (np.cross(e, r) * sm[:, None]).reshape(-1)
        if np.linalg.norm(v) > 1e-8:
            rots.append(v / np.linalg.norm(v))
    if not rots:
        return trans, np.zeros((0, 3 * na))
    # ピボットなし QR の先頭 rank 列は、ランク落ちした入力の像とは限らない。
    # 例: 分子軸が xy 面内だと x/y 回転が従属し、独立な z 回転が最後に来る。
    # SVD で回転空間そのものを取り出す (分子の向きに依存しない)。
    _, singular, vt = np.linalg.svd(np.array(rots), full_matrices=False)
    rank = int(np.sum(singular > 1e-8 * singular[0]))
    return trans, vt[:rank]


def vibrational_space(cart, masses):
    """質量重み座標で並進・回転を除いた振動部分空間の正規直交基底 (n_vib, 3N)。"""
    n3 = 3 * len(cart)
    trans, rots = translation_rotation_basis(cart, masses)
    tr = np.vstack([trans, rots])
    q, _ = np.linalg.qr(tr.T)
    p_vib = np.eye(n3) - q[:, :tr.shape[0]] @ q[:, :tr.shape[0]].T
    w, v = np.linalg.eigh(0.5 * (p_vib + p_vib.T))
    return v[:, w > 0.5].T


def _image_basis(P, p_vib):
    """射影演算子 P の振動空間内での像の正規直交基底 (rank, 3N)。

    正しい射影なら固有値は 0 か 1 だけ。中間の値はクラスの割り当てや指標の誤りを
    意味するので、黙って進まずに例外にする。"""
    pv = p_vib @ P @ p_vib
    w, v = np.linalg.eigh(0.5 * (pv + pv.T))
    if np.any((w > 0.02) & (w < 0.98)) or np.any(w > 1.02):
        raise AnalysisError(Msg("err.projection"))
    return v[:, w > 0.5].T


def _irrep_activity(table, irrep):
    """指標表の linear / quadratic から Raman / IR 活性を判定する。"""
    ir = irrep in table.get("linear", {})
    ra = irrep in table.get("quadratic", {})
    return "Raman+IR" if (ra and ir) else "Raman" if ra else "IR" if ir else "Silent"


# ════════════════════════════════════════════════════════════════════════
# 5. 既約分解 (点群の種類ごと)
#    いずれも [(既約表現名, d, m, 部分空間基底 (m·d, 3N; 質量重み))] を返す。
# ════════════════════════════════════════════════════════════════════════
def _decompose_linear(cart, numbers, masses, has_inversion):
    """直線分子 (C∞v, D∞h) を真の連続群射影 P = (l/h) ∫ χ*(R) D(R) で分解する。

    回転角 φ の各族 (C∞(φ), σv(φ), [S∞(φ), C2'(φ)]) を数値積分する。直線分子の
    力学表現に現れるのは Σ と Π だけ (3 次元変位で Δ 以上は出ない)。
    返り値: (blocks, 分子軸を z にする回転 r_axis, その座標系の座標 cz,
             有限部分群の力学表現 D(R) のリスト)。基底は cz 座標系のもの。"""
    na = len(cart)
    n3 = 3 * na
    axis = _inertia_axes(cart, masses)[:, 0]
    r_axis = rotation_to_z(axis)
    cz = cart @ r_axis.T
    Z = list(numbers)

    def dmat(R):
        D = np.zeros((n3, n3))
        img = cz @ R.T
        for i in range(na):
            for k in range(na):
                if Z[k] == Z[i] and np.linalg.norm(img[i] - cz[k]) < MATCH_TOL:
                    D[3*k:3*k+3, 3*i:3*i+3] = R
                    break
        return D

    def rotz(phi):
        c, s = math.cos(phi), math.sin(phi)
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])

    def reflv(phi):
        n = np.array([math.cos(phi), math.sin(phi), 0.0])
        return np.eye(3) - 2 * np.outer(n, n)

    n_int = 360
    sigma_h = np.diag([1., 1., -1.])
    phis = np.linspace(0, 2 * math.pi, n_int, endpoint=False)
    fam_c = [dmat(rotz(p)) for p in phis]
    fam_sv = [dmat(reflv(p)) for p in phis]
    if has_inversion:
        fam_s = [dmat(sigma_h @ rotz(p)) for p in phis]
        fam_c2 = [dmat(2 * np.outer([math.cos(p), math.sin(p), 0.], [math.cos(p), math.sin(p), 0.])
                       - np.eye(3)) for p in phis]

    def projector(chi_c, chi_sv, chi_s, chi_c2, dim):
        P = np.zeros((n3, n3))
        for k, p in enumerate(phis):
            P += chi_c(p) * fam_c[k] + chi_sv(p) * fam_sv[k]
            if has_inversion:
                P += chi_s(p) * fam_s[k] + chi_c2(p) * fam_c2[k]
        n_fam = 4 if has_inversion else 2
        return (dim / n_fam) * P / n_int

    one = lambda p: 1.0            # noqa: E731
    zero = lambda p: 0.0           # noqa: E731
    two_cos = lambda p: 2 * math.cos(p)   # noqa: E731
    if has_inversion:              # D∞h: Σg⁺, Σu⁺, Πg, Πu
        specs = [("Sigma_g+", 1, one, one, one, one),
                 ("Sigma_u+", 1, one, one, lambda p: -1.0, lambda p: -1.0),
                 ("Pi_g", 2, two_cos, zero, lambda p: -2 * math.cos(p), zero),
                 ("Pi_u", 2, two_cos, zero, two_cos, zero)]
    else:                          # C∞v: Σ⁺, Π
        specs = [("Sigma+", 1, one, one, zero, zero),
                 ("Pi", 2, two_cos, zero, zero, zero)]

    vib = vibrational_space(cz, masses)
    p_vib = vib.T @ vib
    sm = np.sqrt(np.asarray(masses, float))
    blocks = []
    for name, d, c_c, c_sv, c_s, c_c2 in specs:
        P = projector(c_c, c_sv, c_s, c_c2, d)
        w, v = np.linalg.eigh(0.5 * (P + P.T))
        sym = v[:, w > 0.5].T                           # この対称種の全自由度 (Cartesian)
        if sym.shape[0] == 0:
            continue
        sym_mw = (sym.reshape(-1, na, 3) * sm[None, :, None]).reshape(sym.shape[0], n3)
        q, _ = np.linalg.qr(sym_mw.T)
        sym_mw = q[:, :np.linalg.matrix_rank(sym_mw, tol=1e-6)].T
        _, s, vt = np.linalg.svd(sym_mw @ p_vib)        # 振動空間との共通部分
        n_vib = int(np.sum(s > 0.5))
        if n_vib:
            blocks.append((name, d, n_vib // d, vt[:n_vib]))
    # 縮退の解消用に有限部分群 (C12v / D12h; 30° 刻み) の操作も返す
    sub = range(0, n_int, 30)
    sym_mats = [fam_c[k] for k in sub] + [fam_sv[k] for k in sub]
    if has_inversion:
        sym_mats += [fam_s[k] for k in sub] + [fam_c2[k] for k in sub]
    return blocks, r_axis, cz, sym_mats


def _decompose_complex(pg, ops, dmats, p_vib):
    """複素既約表現を持つ群 (Cn, Cnh, Sn, T, Th) の分解。

    Abelian 群は複素指標を符号付き回転角で照合し、E (複素ペア) は片成分の指標で
    重複度を求めて物理 2 次元として扱う。T/Th (real_pseudo) の E は擬 2 次元
    (ノルム 2h) で、統合擬クラスなので符号なし角度で照合する。"""
    CT = COMPLEX_TABLES[pg]
    real_pseudo = CT.get("real_pseudo", False)
    pseudo2 = CT.get("pseudo2", set())
    _, paxis = _principal_axis(ops)
    h = len(ops)
    col = []
    for R, _ in ops:
        key = op_signed_key(R, paxis)
        ci = None
        for idx, (kind, ang) in enumerate(CT["fine"]):
            if kind != key[0]:
                continue
            if real_pseudo:
                hit = abs(abs(ang) - abs(key[1])) < 3.0
            else:
                hit = abs(((ang - key[1] + 180) % 360) - 180) < 3.0
            if hit:
                ci = idx
                break
        col.append(ci)
    chi_vib = np.array([np.trace(p_vib @ D) for D in dmats])
    blocks = []
    for name, chars in CT["ctab"].items():
        d = CT["dims"][name]
        chi = np.array([chars[c] if c is not None else 0 for c in col], dtype=complex)
        norm_h = 2 * h if (real_pseudo and name in pseudo2) else h
        m = int(round((np.sum(np.conj(chi) * chi_vib) / norm_h).real))
        if m <= 0:
            continue
        if real_pseudo:            # 実指標がすでに共役ペアの和
            weight, scale = np.real(chi), d / norm_h
        elif d == 2:               # 共役ペアの射影の和: (1/h) Σ (χ_a + χ_b) D = (1/h) Σ 2Re χ_a D
            weight, scale = 2 * np.real(chi), 1.0 / h
        else:
            weight, scale = np.real(chi), 1.0 / h
        P = scale * sum(w * D for w, D in zip(weight, dmats))
        blocks.append((name, d, m, _image_basis(P, p_vib)))
    return blocks, CT


def _decompose_icosahedral(pg, ops, dmats, p_vib):
    """正二十面体群 (I, Ih) の分解。主軸が一意でないが、回転角だけでクラスが
    決まるので、角度キーで指標を照合して射影する。"""
    T = CHAR_TABLES[pg]
    table = np.array(T["table"], float)
    keys = [op_angle_key(R) for R, _ in ops]
    h = len(ops)
    chi_vib = np.array([np.trace(p_vib @ D) for D in dmats])
    blocks = []
    for irow, name in enumerate(T["irreps"]):
        by_key = {ICOSAHEDRAL_CLASS_KEYS[c]: table[irow][ci] for ci, c in enumerate(T["classes"])}
        chi = np.array([by_key.get(k, 0.0) for k in keys])
        m = int(round(np.sum(chi * chi_vib) / h))
        if m <= 0:
            continue
        d = int(round(table[irow][0]))
        P = (d / h) * sum(c * D for c, D in zip(chi, dmats))
        blocks.append((name, d, m, _image_basis(P, p_vib)))
    return blocks


def _decompose_geometric(pg, cart, ops, dmats, p_vib, frame=None):
    """各操作を共役類 (指標表の列) に割り当てて射影する (実指標の群)。

    C2v / D2 / D2h は分子固定軸 frame に対する向きでクラス名を決め (B1/B2/B3 の慣例)、
    それ以外は群の積表から求めた共役類を指標表に対応させる (_assign_classes)。
    どちらの場合も割り当てが指標表と矛盾しないことを確かめてから射影する。"""
    T = CHAR_TABLES[pg]
    classes = T["classes"]
    table = np.array(T["table"], float)
    h = len(ops)
    if frame is not None:                                  # C2v / D2 / D2h
        names = [_abelian_class(R, frame, pg) for R, _ in ops]
        if any(n not in classes for n in names):
            raise AnalysisError(Msg("err.class_assign", pg=pg))
        col = [classes.index(n) for n in names]
        mult, inv, _ = _group_structure(ops)
        if not _valid_class_function(table, col, mult, inv):
            raise AnalysisError(Msg("err.class_assign", pg=pg))
    else:
        col = _assign_classes(pg, ops, cart)
    n_vib = int(round(np.trace(p_vib)))
    blocks, dim_sum = [], 0
    for irow, name in enumerate(T["irreps"]):
        d = int(round(table[irow][0]))
        chi = [table[irow][c] for c in col]
        P = (d / h) * sum(c * D for c, D in zip(chi, dmats))
        basis = _image_basis(P, p_vib)
        if basis.shape[0] == 0:
            continue
        if basis.shape[0] % d:
            raise AnalysisError(Msg("err.projection"))
        dim_sum += basis.shape[0]
        blocks.append((name, d, basis.shape[0] // d, basis))
    if dim_sum != n_vib:
        raise AnalysisError(Msg("err.projection"))
    return blocks


# ════════════════════════════════════════════════════════════════════════
# 6. メイン解析
# ════════════════════════════════════════════════════════════════════════
def analyze_molecule(cart, numbers, symbols=None, tol=0.1, notes=None):
    """分子の振動モードを解析して AnalysisResult を返す (ファイル出力はしない)。

    cart    : Cartesian 座標 [Å] (N, 3)
    numbers : 原子番号
    tol     : 対称性判定の許容差 [Å]。この範囲のずれは対称化して吸収する。
    notes   : 結果の注記の先頭に加える文 (読み込み時の注記など)"""
    from pymatgen.core import Molecule
    from pymatgen.symmetry.analyzer import PointGroupAnalyzer

    cart = np.array(cart, float)
    numbers = [int(z) for z in numbers]
    symbols = list(symbols) if symbols else [element_symbol(z) for z in numbers]
    na = len(numbers)
    if na < 2:
        raise AnalysisError(Msg("err.min_atoms"))
    if not (math.isfinite(tol) and tol > 0):
        raise AnalysisError(Msg("err.tol"))
    masses = [atomic_mass(z) for z in numbers]
    m_arr = np.asarray(masses)
    # 重心を原点に (対称操作は原点を通る前提)
    cart = cart - (m_arr[:, None] * cart).sum(0) / m_arr.sum()
    notes = list(notes or [])

    # ── 対称化と点群判定 ──
    linear = is_linear(cart, masses, tol)
    if linear:
        sym_cart, has_inv = _symmetrize_linear(cart, numbers, tol)
        pg = "D*h" if has_inv else "C*v"
        ops = None
    else:
        # 許容差 tol で操作を探して平均 → 対称化した構造で探し直す、を操作数が
        # 変わらなくなるまで繰り返し、最後に厳密な許容差で操作を確定する
        sym_cart, n_prev = cart, -1
        for _ in range(5):
            ops = find_point_group_operations(sym_cart, numbers, tol)
            if len(ops) == n_prev:
                break
            n_prev = len(ops)
            sym_cart = symmetrize_by_group(cart, ops)
            sym_cart -= (m_arr[:, None] * sym_cart).sum(0) / m_arr.sum()
        # 歪んだ構造から求めた操作は厳密な群ではないので、「厳密な許容差で操作を取り直し →
        # その操作で平均」を収束するまで繰り返し、構造を数値精度で厳密に対称にする
        for _ in range(20):
            ops = find_point_group_operations(sym_cart, numbers, MATCH_TOL)
            refined = symmetrize_by_group(sym_cart, ops)
            refined -= (m_arr[:, None] * refined).sum(0) / m_arr.sum()
            converged = np.abs(refined - sym_cart).max() < 1e-12
            sym_cart = refined
            if converged:
                break
        ops = find_point_group_operations(sym_cart, numbers, MATCH_TOL)
        # 点群記号は pymatgen に任せる (許容差をこちらとそろえる)
        # (仮想原子は pymatgen が扱えないので、種類の区別だけ保った実在元素に置き換える)
        kinds = sorted(set(numbers))
        stand_in = [kinds.index(z) + 1 if any(k >= 119 for k in kinds) else z for z in numbers]
        # 対称化した構造なので小さい許容差でも判定できる。見つけた操作の数 (群の位数) と
        # 合う記号が得られるまで許容差を変えて試す (細長い分子を C∞v と誤判定するなどを防ぐ)
        pg = None
        for tol_pg in (0.01, 0.05, max(tol, 1e-3), 1e-3):
            cand = PointGroupAnalyzer(Molecule(stand_in, sym_cart), tolerance=tol_pg).sch_symbol
            if pg is None:
                pg = cand
            if _group_order(cand) == len(ops):
                pg = cand
                break
        else:                                  # どの許容差でも合わなければ操作から決める
            pg = _point_group_from_ops(ops) or pg
    shift = float(np.max(np.linalg.norm(sym_cart - cart, axis=1)))
    if shift > 1e-4:
        notes.append(Msg("mol.symmetrized", pg=pg, shift=shift))
    cart = sym_cart

    n_vib_expected = 3 * na - (5 if linear else 6)
    result = AnalysisResult(kind="molecule", numbers=numbers, symbols=symbols, cart=cart,
                            modes=[], point_group=pg, is_linear=linear, notes=notes)
    if n_vib_expected <= 0:
        return result

    if linear:
        blocks, r_axis, cz, sym_mats = _decompose_linear(cart, numbers, masses, pg == "D*h")
        result.modes = _blocks_to_modes(
            [(n, d, m, b, m > 1) for n, d, m, b in blocks], cz, masses,
            label_of=lambda n: _LINEAR_DISPLAY.get(n, n),
            activity_of=lambda n: _LINEAR_ACTIVITY[n], notes=notes, sym_mats=sym_mats,
            back_rotation=r_axis)
        result.notes.append(Msg("mol.linear", g="D∞h" if pg == "D*h" else "C∞v"))
        return _finish(result)

    order = _group_order(pg)
    if order is None:
        result.notes.append(Msg(
            "mol.no_table", pg=pg, n=n_vib_expected,
            avail=", ".join(sorted(set(CHAR_TABLES) | set(COMPLEX_TABLES))) + ", C∞v, D∞h"))
        return result
    if len(ops) != order:
        raise AnalysisError(Msg("err.order_mismatch", pg=pg, order=order, n=len(ops)))

    vib = vibrational_space(cart, masses)
    p_vib = vib.T @ vib
    dmats = mechanical_matrices(ops, na)

    if pg in COMPLEX_TABLES:
        blocks, CT = _decompose_complex(pg, ops, dmats, p_vib)
        table = CT
    else:
        table = CHAR_TABLES[pg]
        frame = None
        if pg in ("C2v", "D2", "D2h"):
            frame, frame_note = _abelian_frame(pg, ops, cart, np.asarray(masses, float))
            notes.append(Msg("mol.axes", frame=frame_note))
        if pg in ("I", "Ih"):
            blocks = _decompose_icosahedral(pg, ops, dmats, p_vib)
        else:
            blocks = _decompose_geometric(pg, cart, ops, dmats, p_vib, frame)
    result.modes = _blocks_to_modes(
        [(n, d, m, b, m > 1) for n, d, m, b in blocks], cart, masses,
        label_of=display_label,
        activity_of=lambda n: (_irrep_activity(table, n), table.get("quadratic", {}).get(n, [])),
        notes=notes, sym_mats=dmats)
    if any(m.tensor for m in result.modes) and pg not in ("C1", "Ci"):
        # テンソル成分は指標表の標準の向き (z = 主軸) のもので、入力座標の向きとは無関係
        notes.append(Msg("mol.tensor_frame"))
    return _finish(result)


def _blocks_to_modes(blocks, cart, masses, label_of, activity_of, notes, sym_mats=None,
                     back_rotation=None, refine=True):
    """[(名前, d, m, 基底 (m·d, 3N; 質量重み), refined)] を Mode のリストへ。

    m>1 のブロックは網目モデルで部分空間内の向きを推測してから d 本ずつに分ける。
    向きが原理的に決まらないブロックがあれば notes に書き足す。"""
    na = len(cart)
    sm = np.sqrt(np.asarray(masses, float))
    cart_out = np.asarray(cart, float) @ back_rotation if back_rotation is not None else cart
    modes, arbitrary = [], []
    hessians, tiebreaker = None, None
    for name, d, m, basis, refined in blocks:
        if refine and m > 1 and basis.shape[0] > 0:
            if hessians is None:
                hessians = (network_hessian(cart), network_hessian(cart, isotropic=True))
                if sym_mats is not None:
                    tiebreaker = invariant_tiebreaker(sym_mats)
            lam, basis = refine_in_subspace(basis, cart, masses, return_eigenvalues=True,
                                            sym_mats=sym_mats, hessians=hessians,
                                            tiebreaker=tiebreaker)
            # 網目の剛性がゼロの部分空間 (平面分子の面外モードなど) では向きは任意
            if np.max(np.abs(lam)) < 1e-6:
                arbitrary.append(label_of(name))
        activity, tensor = activity_of(name)
        for c in range(m):
            blk = basis[c * d:(c + 1) * d]
            if blk.shape[0] < d:
                continue
            disp = [row.reshape(na, 3) / sm[:, None] for row in blk]   # 質量重みを外す
            if back_rotation is not None:
                disp = [u @ back_rotation for u in disp]               # 元の座標系へ
            # 縮退成分は x, y, z 軸を基準に選び直す (縮退なしでも符号をそろえる)
            disp = axis_aligned_components(disp, np.eye(3), sort_keys=cart_out)
            modes.append(Mode(label=label_of(name), degeneracy=d, activity=activity,
                              tensor=list(tensor), displacements=disp, refined=refined))
    if arbitrary:
        notes.append(Msg("mol.arbitrary", labels=", ".join(arbitrary)))
    return modes


def _finish(result):
    result.decomposition = _decomposition_string(result.modes)
    return result


def _decomposition_string(modes):
    counts = Counter(m.label for m in modes)
    seen = []
    for m in modes:
        if m.label not in seen:
            seen.append(m.label)
    return "Γvib = " + " + ".join(f"{counts[s]}{s}" if counts[s] > 1 else s for s in seen)


# ════════════════════════════════════════════════════════════════════════
# 7. 出力 (テキストレポート / XSF)
# ════════════════════════════════════════════════════════════════════════
def format_report(result: AnalysisResult) -> str:
    """CLI に表示するテキストレポートを現在の言語で作る。"""
    na = len(result.numbers)
    lin = result.is_linear
    lines = ["=" * 72, t("rep.mol.title"), "=" * 72, "", t("rep.mol.molecule"),
             t("rep.mol.atoms", n=na, pg=result.point_group,
               lin=t("rep.mol.linear_tag") if lin else ""),
             t("rep.mol.dof", n3=3 * na, r=2 if lin else 3, v=3 * na - (5 if lin else 6))]
    for note in result.notes:
        lines.append(t("rep.note", v=note))
    modes = result.modes
    if not modes:
        return "\n".join(lines)
    lines += ["", t("rep.mol.internal", dof=sum(m.degeneracy for m in modes)),
              f"  {result.decomposition}"]
    for a in ACTIVITIES:
        sub = [m for m in modes if m.activity == a]
        if sub:
            lines.append(t("rep.activity_line", sym=ACTIVITY_SYMBOL[a], act=a, n=len(sub),
                           dof=sum(m.degeneracy for m in sub)))
    widths = (6, 10, 6, 10, 17)
    cols = t("rep.mol.cols").split("|")
    lines += ["", "━" * 72,
              " " + "".join(pad(c, w) for c, w in zip(cols, widths)) + cols[-1], "━" * 72]
    for i, m in enumerate(modes, 1):
        cells = (i, m.label, m.degeneracy, m.activity, m.direction_label)
        lines.append(" " + "".join(pad(c, w) for c, w in zip(cells, widths))
                     + (", ".join(m.tensor) if m.tensor else "—"))
    lines += ["━" * 72] + t("rep.mol.legend").splitlines()
    return "\n".join(lines)


def molecule_box(cart, margin=15.0):
    """孤立分子を XSF で表示するための大きな箱 (格子, 箱内座標)。"""
    span = cart.max(0) - cart.min(0)
    lattice = np.diag(np.maximum(span + margin, margin))
    shifted = cart + (np.diag(lattice) - span) / 2 - cart.min(0)
    return lattice, shifted


def write_mode_files(result: AnalysisResult, outdir, scale=1.0, formats=("xsf", "vesta")):
    """各モード (縮退成分ごと) を XSF と VESTA 形式で書き出し、書いたパスのリストを返す。

    XSF は最大変位 scale Å の矢印 (見た目は VESTA 任せ)。VESTA 形式は灰色の細い矢印
    (中点が原子の中心、最大で中心から ±0.53 Å) と結合の設定まで書き込む。
    formats で書き出す形式を選べる (("xsf",), ("vesta",), 両方)。分子は大きな箱 (周期像が見えない) に入れる。"""
    ensure_outdir(outdir)
    lattice, cart = molecule_box(result.cart)
    frac = cart @ np.linalg.inv(lattice)
    labels = [f"{s}{i}" for i, s in enumerate(result.symbols, 1)]
    paths = []
    for i, m in enumerate(result.modes, 1):
        for c, disp in enumerate(m.displacements):
            stem = os.path.join(outdir, mode_file_stem(i, m, c))
            header = xsf_header(i, m, c, "true direction requires quantum chemistry")
            if "xsf" in formats:
                write_xsf(stem + ".xsf", lattice, cart, result.numbers, normalize_max(disp),
                          header, scale=scale)
                paths.append(stem + ".xsf")
            if "vesta" in formats:
                write_vesta(stem + ".vesta", lattice, frac, result.numbers, labels, disp,
                            title=" / ".join(header[:3]), arrow_length=0.53, molecule=True)
                paths.append(stem + ".vesta")
    return paths
