"""結晶の因子群解析 (Γ 点, k=0)。

CIF を入れると空間群を判定し、Γ 点の格子振動を既約表現に分解する。
Mulliken ラベル、Raman/IR 活性、偏光ラマン選択則、各モードの変位パターンを返す。

────────────────────────────────────────────────────────────────────────
処理の流れ (analyze_crystal)
────────────────────────────────────────────────────────────────────────
 1. spgrep で Γ 点の既約表現を得る (spgrep は入力セルの操作で表現を作る)。
 2. 中心格子 (I, F, C, R) は primitive cell に変換し、既約表現を primitive の
    点群操作に縮約する (純粋並進ぶんの重複操作を除く)。
 3. 複素共役ペアの 1 次元既約表現は、実空間では 1 つの 2 次元モードとして
    現れるのでまとめる (ラベル・活性の判定用)。
 4. 力学表現を各既約表現に射影して対称性適応変位の部分空間を得る。そのうち
    剛体並進の成分を音響モード、残りを光学モードとし、光学コピーが複数 (m>1)
    あれば網目モデルで部分空間内の向きを推測する。

依存: numpy, spglib, spgrep, pymatgen
"""
from __future__ import annotations

import os
from collections import Counter

import numpy as np
import spglib
from spgrep import get_spacegroup_irreps

try:                                    # spgrep >= 0.4 の新しい場所
    from spgrep.rep.representation import get_character
except ImportError:                     # 古い spgrep
    from spgrep.representation import get_character

from .common import (AnalysisError, AnalysisResult, Mode, ACTIVITIES, ACTIVITY_SYMBOL,
                     atomic_mass, element_symbol, ensure_outdir, mode_file_stem,
                     axis_aligned_components, invariant_tiebreaker, network_hessian, normalize_max,
                     refine_in_subspace, site_element,
                     write_vesta, write_xsf, xsf_header)
from .i18n import Msg, pad, t

# 網目モデルのカットオフ [Å] (±1 セルの像まで数える)
NETWORK_CUTOFF = 4.0


# ════════════════════════════════════════════════════════════════════════
# 1. 結晶構造の入力
# ════════════════════════════════════════════════════════════════════════
def load_crystal(path, notes=None):
    """CIF などから (格子 (行ベクトル), 分率座標, 原子番号, 元素名) を返す。

    CIF は自前の読み込み (cif_reader) で、部分占有・ラベル・複数ブロックなどに対応する。
    notes (リスト) を渡すと、読み込み時の注記 (部分占有をまとめた等) を書き足す。"""
    if os.path.splitext(path)[1].lower() == ".cif":
        from .cif_reader import read_cif
        c = read_cif(path)
        if notes is not None:
            notes.extend(c.notes)
        return c.lattice, c.frac, c.numbers, c.names
    from pymatgen.core import Structure
    try:
        s = Structure.from_file(path)
    except Exception as e:
        raise AnalysisError(Msg("err.read_crystal", path=path, err=e)) from e
    elems = [site_element(site) for site in s]
    return s.lattice.matrix, s.frac_coords, [z for z, _ in elems], [e for _, e in elems]


def graphite_2h():
    """デモ用の黒鉛 2H (P6_3/mmc)。"""
    a, c = 2.461, 6.708
    lat = np.array([[a, 0, 0], [-a / 2, a * np.sqrt(3) / 2, 0], [0, 0, c]])
    pos = np.array([[0, 0, .25], [0, 0, .75], [1/3, 2/3, .25], [2/3, 1/3, .75]])
    return lat, pos, [6, 6, 6, 6], ["C", "C", "C", "C"]


# ════════════════════════════════════════════════════════════════════════
# 2. 対称操作の分類
# ════════════════════════════════════════════════════════════════════════
def build_R_carts(lattice, rots):
    """分率座標の回転行列を Cartesian の回転行列に変換する。"""
    li = np.linalg.inv(lattice)
    return np.array([lattice.T @ r @ li.T for r in rots])


_ORDER_BY_ANGLE = {0: 1, np.pi: 2, 2 * np.pi / 3: 3, np.pi / 2: 4, np.pi / 3: 6}


def _classify(R, tol=0.1):
    """Cartesian 操作 R を (次数, 軸, 真回転か) に分類する。

    次数は真回転部分 (回映なら −R) の回転角から (結晶では 1,2,3,4,6 のみ)。
    鏡映は「−R が 2 回回転」なので (2, 法線, False)、反転は (1, -, False)。"""
    proper = np.linalg.det(R) > 0
    rp = R if proper else -R
    theta = np.arccos(np.clip((np.trace(rp).real - 1) / 2, -1, 1))
    order = next((o for ang, o in _ORDER_BY_ANGLE.items() if abs(theta - ang) < tol), 1)
    if order == 1:
        return order, np.array([0., 0., 1.]), proper
    ev, evc = np.linalg.eig(rp)
    axis = evc[:, np.argmin(np.abs(ev - 1))].real
    n = np.linalg.norm(axis)
    if n > 1e-10:
        axis /= n
    for k in (2, 1, 0):                  # 向きをそろえる (z, y, x の順で正に)
        if abs(axis[k]) > 0.05:
            if axis[k] < 0:
                axis = -axis
            break
    return order, axis, proper


def _par(v1, v2, tol=0.1):
    n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
    if n1 < 1e-10 or n2 < 1e-10:
        return False
    return abs(abs(np.dot(v1 / n1, v2 / n2)) - 1) < tol


def _perp(v1, v2, tol=0.1):
    n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
    if n1 < 1e-10 or n2 < 1e-10:
        return False
    return abs(np.dot(v1 / n1, v2 / n2)) < tol


def _is_cubic(ops):
    """3 回軸が 4 本以上あれば立方晶系の点群。"""
    c3 = []
    for order, ax, proper in ops:
        if proper and order == 3 and not any(_par(ax, e) for e in c3):
            c3.append(ax)
    return len(c3) >= 4


# ════════════════════════════════════════════════════════════════════════
# 3. 既約表現の前処理 (複素共役ペアの統合, primitive への縮約)
# ════════════════════════════════════════════════════════════════════════
def _is_real_character(chi, tol=1e-6):
    return all(abs(c.imag) < tol for c in chi)


def _find_conjugate(chars, i, used, tol=0.01):
    """chars[i] と複素共役な指標を持つ未使用の既約表現の番号 (なければ None)。"""
    return next((j for j in range(len(chars)) if j != i and not used[j]
                 and np.all(np.abs(chars[i] - np.conj(chars[j])) < tol)), None)


def merge_complex_pairs(irreps):
    """1 次元の複素共役ペアを 2 次元の実表現にまとめる (ラベル・活性判定用)。"""
    chars = [get_character(irr) for irr in irreps]
    used = [False] * len(irreps)
    result = []
    for i, irr in enumerate(irreps):
        if used[i]:
            continue
        used[i] = True
        if _is_real_character(chars[i]) or irr.shape[1] > 1:
            result.append(irr)
            continue
        j = _find_conjugate(chars, i, used)
        if j is None:
            result.append(irr)
            continue
        used[j] = True
        d = irr[:, 0, 0]                                     # e^{iθ} → 2×2 回転行列
        result.append(np.array([[[x.real, -x.imag], [x.imag, x.real]] for x in d],
                               dtype=complex))
    return result


def reduce_irrep_to_pg(irr, rots_conv, conv_lat, rots_prim, prim_lat):
    """spgrep の既約表現 (入力セルの操作) を primitive の点群操作に縮約する。

    中心格子の conventional cell には純粋並進ぶんの重複操作があるので、
    Cartesian の回転行列で照合して点群操作ぶんだけ取り出す。"""
    rc_conv = build_R_carts(conv_lat, rots_conv)
    rc_prim = build_R_carts(prim_lat, rots_prim)
    reduced = []
    for rp in rc_prim:
        match = next((j for j, rc in enumerate(rc_conv) if np.allclose(rp, rc, atol=1e-2)), None)
        if match is None:
            raise AnalysisError(Msg("err.reduce"))
        reduced.append(irr[match])
    return np.array(reduced)


# ════════════════════════════════════════════════════════════════════════
# 4. Mulliken ラベリング (全 32 点群対応)
# ════════════════════════════════════════════════════════════════════════
_SUFFIXES = ("''", "'", "g", "u")


def _split_suffix(label):
    """'B1g' → ('B1', 'g'), "E''" → ('E', "''")。"""
    for s in _SUFFIXES:
        if label.endswith(s):
            return label[: -len(s)], s
    return label, ""


def mulliken_labels(irreps, R_carts):
    """各既約表現に Mulliken 記号 (A1g, E', T2u, ...) を付ける。

    主軸 (最高次の回転軸) まわりの指標で A/B、主軸に垂直な C2' または σv で
    添字 1/2、反転で g/u、σh で '/'' を決め、最後に区別の要らない添字を外す。"""
    ops = [_classify(R) for R in R_carts]
    cubic = _is_cubic(ops)
    max_n = max(o for o, _, proper in ops if proper)
    max_sn = max((o for o, _, proper in ops if not proper and o > 2), default=0)
    # 主軸 = 最高次の真回転軸。ただし回映軸 (−4, −6) の次数がそれより高ければ回映軸を主軸にとる
    # (D2d や S4 では C2 軸が複数あり、操作の並び順によっては C2' を主軸と取り違えるため)
    if max_sn > max_n:
        pa = next(ax for o, ax, proper in ops if not proper and o == max_sn)
    else:
        pa = next((ax for o, ax, proper in ops if proper and o == max_n and max_n > 1),
                  np.array([0., 0., 1.]))

    def first(pred):
        return next((i for i, op in enumerate(ops) if pred(*op)), None)

    i_cn = first(lambda o, ax, p: p and o == max_n and max_n > 1 and _par(ax, pa))
    i_sn = first(lambda o, ax, p: not p and o == max_sn and max_sn > 2 and _par(ax, pa))
    i_c2p = [i for i, (o, ax, p) in enumerate(ops) if p and o == 2 and _perp(ax, pa)]
    i_sv = [i for i, (o, ax, p) in enumerate(ops) if not p and o == 2 and _perp(ax, pa)]
    i_sh = first(lambda o, ax, p: not p and o == 2 and _par(ax, pa))
    i_inv = first(lambda o, ax, p: not p and o == 1)
    has_inv = i_inv is not None

    # Cs: 回転軸がなく鏡映だけ → その鏡映を σh として扱う
    if i_sh is None and not has_inv and max_n <= 1:
        i_sh = first(lambda o, ax, p: not p and o == 2)
        if i_sh is not None:
            i_sv = []
    has_sh = i_sh is not None

    # 添字 1/2 を決める基準の操作 (Bilbao の慣例; 結晶の標準の向き x ∥ a, z ∥ c で):
    #   4 回軸の群: C2' ∥ a (x)、σ の法線 ∥ a      6 回軸の群: C2' ⊥ a (y)、σ の法線 ∥ a
    #   mm2: σ の法線 ∥ b (y)                       その他 (3 回軸など): 同値なので最初のもの
    def along(indices, k):
        return next((i for i in indices if abs(abs(ops[i][1][k]) - 1) < 0.05), None)

    ref_c2p = (along(i_c2p, 0) if max_n == 4 else along(i_c2p, 1) if max_n == 6 else None)
    if ref_c2p is None and i_c2p:
        ref_c2p = i_c2p[0]
    ref_sv = (along(i_sv, 1) if max_n == 2 else along(i_sv, 0) if max_n in (4, 6) else None)
    if ref_sv is None and i_sv:
        ref_sv = i_sv[0]

    # D2 / D2h (互いに直交する 3 本の C2; S4 を持つ D2d は除く)
    is_d2 = max_n == 2 and len(i_c2p) > 0 and not cubic and max_sn <= 2
    c2_list = []
    if is_d2:
        uniq = []
        for i, (o, ax, p) in enumerate(ops):
            if p and o == 2 and not any(_par(ax, ops[j][1]) for j, _ in uniq):
                uniq.append((i, ax))
        uniq.sort(key=lambda t: -int(np.argmax(np.abs(t[1]))))   # 標準の向きで z, y, x 軸の順
        c2_list = [i for i, _ in uniq]

    # 立方晶: 添字 1/2 を決める操作 (C4 軸以外の C2、なければ σd)
    face_c2 = []
    if cubic:
        c4_axes = [ax for o, ax, p in ops if p and o == 4]
        if c4_axes:
            face_c2 = [i for i, (o, ax, p) in enumerate(ops)
                       if p and o == 2 and not any(_par(ax, c) for c in c4_axes)]
        else:
            face_c2 = [i for i, (o, ax, p) in enumerate(ops) if not p and o == 2]

    raw = []
    for irrep in irreps:
        ch = get_character(irrep).real
        d = irrep.shape[1]
        if d == 1:
            if cubic:
                letter, sub = "A", (("1" if ch[face_c2[0]] > 0 else "2") if face_c2 else "")
            elif is_d2:
                signs = [ch[c2_list[k]] > 0 for k in range(min(3, len(c2_list)))]
                if all(signs):
                    letter, sub = "A", ""
                elif len(signs) >= 3:
                    letter, sub = "B", ("1" if signs[0] else "2" if signs[1] else "3")
                else:
                    letter, sub = ("A" if signs[0] else "B"), ""
            else:
                cn = ch[i_cn] if i_cn is not None else 1.0
                letter = "A" if (max_n <= 1 or abs(cn - 1) < 0.5) else "B"
                if max_n <= 2 and i_sn is not None and abs(ch[i_sn]) > 0.5:
                    letter = "A" if ch[i_sn] > 0 else "B"      # S4 など回映が主軸
                sub = ""
                if ref_c2p is not None:
                    sub = "1" if ch[ref_c2p] > 0 else "2"
                elif ref_sv is not None:
                    sub = "1" if ch[ref_sv] > 0 else "2"
            label = letter + sub
        elif d == 2:
            sub = ""
            if not cubic and max_n >= 4 and i_cn is not None:
                # χ(Cn) = 2cos(2πk/n) となる k が添字 (E1, E2, ...)
                k_best = min(range(1, max_n),
                             key=lambda k: abs(ch[i_cn] - 2 * np.cos(2 * np.pi * k / max_n)))
                sub = str(k_best)
            label = "E" + sub
        elif d == 3:
            if cubic and face_c2:
                sub = "1" if ch[face_c2[0]] < 0 else "2"
            elif i_c2p:
                sub = "1" if ch[i_c2p[0]] < 0 else "2"
            else:
                sub = ""
            label = "T" + sub
        else:
            label = f"H{d}"
        if has_inv:
            label += "g" if ch[i_inv] > 0 else "u"
        elif has_sh:
            label += "'" if ch[i_sh] > 0 else "''"
        raw.append(label)

    # 後処理: 同じ文字 (と g/u/'/'') の既約表現が 1 つしかなければ添字の数字を外す
    count = Counter()
    for label in raw:
        base, suffix = _split_suffix(label)
        count[(base[0], suffix)] += 1
    final = []
    for label in raw:
        base, suffix = _split_suffix(label)
        if base[0] in "ABET" and count[(base[0], suffix)] <= 1 and len(base) > 1:
            label = base[0] + suffix
        final.append(label)
    return final


# ════════════════════════════════════════════════════════════════════════
# 5. 活性とテンソル成分 (指標公式)
# ════════════════════════════════════════════════════════════════════════
def _mult(chi, target):
    """規格化重複度 Σχ*·target / Σ|χ|² (複素統合した可約表現でも正しい)。"""
    norm = np.sum(np.abs(chi) ** 2).real
    return (np.sum(np.conj(chi) * target) / norm).real


def activity(irrep, R_carts):
    """IR 活性 = ベクトル (x,y,z) の表現を含む / Raman 活性 = 対称 2 階テンソルを含む。"""
    chi = get_character(irrep)
    chi_vec = np.array([np.trace(R) for R in R_carts])
    chi_sym2 = (chi_vec ** 2 + np.array([np.trace(R @ R) for R in R_carts])) / 2
    ir = abs(_mult(chi, chi_vec)) > 0.4
    ra = abs(_mult(chi, chi_sym2)) > 0.4
    return "Raman+IR" if (ra and ir) else "Raman" if ra else "IR" if ir else "Silent"


TENSOR_COMPONENTS = ["xx", "yy", "zz", "xy", "xz", "yz"]
_TENSOR_PAIRS = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]


def _sym2_basis():
    """対称 2 階テンソル空間の正規直交基底 (xx, yy, zz, (xy+yx)/√2, ...)。"""
    basis = []
    for i, j in _TENSOR_PAIRS:
        T = np.zeros((3, 3))
        if i == j:
            T[i, i] = 1.0
        else:
            T[i, j] = T[j, i] = 1 / np.sqrt(2)
        basis.append(T)
    return basis


_SYM2 = _sym2_basis()


def sym2_rep(R):
    """対称 2 階テンソル (6 次元) 上での直交操作 R の表現行列 T(R): T_ab = ⟨e_a, R e_b Rᵀ⟩。"""
    return np.array([[np.sum(ea * (R @ eb @ R.T)) for eb in _SYM2] for ea in _SYM2])


def tensor_components(irreps, labels, R_carts):
    """各既約表現 (ラベル) に属する Raman テンソル成分 {label: ['xx', ...]}。

    射影演算子 P = (1/Σ|χ|²) Σ_g χ*(g) T(g) の対角要素が 0 でない成分 (閾値 1e-6)。"""
    T = [sym2_rep(R) for R in R_carts]
    res = {}
    for irr, lbl in zip(irreps, labels):
        chi = get_character(irr)
        P = sum(np.conj(c) * t for c, t in zip(chi, T)) / np.sum(np.abs(chi) ** 2)
        comps = [TENSOR_COMPONENTS[a] for a in range(6) if abs(P[a, a]) > 1e-6]
        if comps:
            res[lbl] = comps
    return res


def raman_tensor_forms(rows, rep, R_std):
    """モードの縮退成分 rows (d, 3N) に対する Raman テンソルの形 (文字列) のリスト。

    成分 k と同じ変換をするテンソルの空間を射影演算子
        P_k = (1/h) Σ_g ⟨u_k|D(g)|u_k⟩ T(g)      (T(g): 対称テンソル空間上の表現, 標準の向き)
    の像として求める。像の基底を既約行階段形にし、各基底テンソルに文字 a, b, c, … を当てる。
    成分どうしは直交するとは限らない (六方晶の a, b 軸) ので、成分 1 を基準に順に直交化した
    組 (x 型・y 型の partner) について求める。縮退した空間内の向きは任意なので、物理的に意味が
    あるのはこの組が張るテンソルの空間であり、文献 (Loudon 1964) の形と同じ組になる。"""
    q, _ = np.linalg.qr(np.asarray(rows, float).T)          # 成分 1 の向きを保って直交化
    T = [sym2_rep(R) for R in R_std]
    h = len(rep)
    forms = []
    for k in range(q.shape[1]):
        u = q[:, k]
        P = sum((u @ D @ u) * t for D, t in zip(rep, T)) / h
        w, v = np.linalg.eigh(0.5 * (P + P.T))
        forms.append(format_tensor(v[:, w > 1e-6].T))
    return forms


def _rref(B, tol=1e-9):
    """既約行階段形 (ピボットの成分を係数 1 にそろえる)。"""
    B = np.array(B, float)
    r = 0
    for c in range(B.shape[1]):
        if r >= B.shape[0]:
            break
        k = r + int(np.argmax(np.abs(B[r:, c])))
        if abs(B[k, c]) < tol:
            continue
        B[[r, k]] = B[[k, r]]
        B[r] /= B[r, c]
        for i in range(B.shape[0]):
            if i != r:
                B[i] -= B[i, c] * B[r]
        r += 1
    return B[:r]


def format_tensor(basis6):
    """正規直交テンソル基底 _SYM2 の座標群 → Cartesian テンソルの文字列表現。"""
    from fractions import Fraction
    if len(basis6) == 0:
        return ""
    # _SYM2 のせん断基底は (xy+yx)/√2。実際の行列要素に直してから
    # 行階段形にしないと、対角・非対角成分の係数比が √2 ずれる。
    B = _rref(np.asarray(basis6, float) / np.array([1, 1, 1, np.sqrt(2), np.sqrt(2), np.sqrt(2)]))
    letters = "abcdefgh"
    entry = {}
    for ci, comp in enumerate(TENSOR_COMPONENTS):
        terms = []
        for k in range(len(B)):
            x = B[k, ci]
            if abs(x) < 1e-9:
                continue
            f = Fraction(x).limit_denominator(12)
            coef = (f"{x:.3g}" if abs(float(f) - x) > 1e-6
                    else "" if f == 1 else "-" if f == -1 else str(f))
            terms.append(f"{coef}{letters[k]}")
        entry[comp] = "+".join(terms).replace("+-", "-") if terms else "·"
    e = entry
    return (f"[[{e['xx']},{e['xy']},{e['xz']}],[{e['xy']},{e['yy']},{e['yz']}],"
            f"[{e['xz']},{e['yz']},{e['zz']}]]")


# ════════════════════════════════════════════════════════════════════════
# 6. 力学表現と射影
# ════════════════════════════════════════════════════════════════════════
def _atom_permutation(lattice, positions, numbers, rotation, translation, tol):
    """周期境界と原子種を含めて操作の全単射を照合する (perm[i] = 写り先)。"""
    positions = np.asarray(positions, float)
    numbers = np.asarray(numbers)
    image = positions @ np.asarray(rotation).T + translation
    delta = image[:, None, :] - positions[None, :, :]
    delta -= np.round(delta)
    match = ((numbers[:, None] == numbers[None, :])
             & (np.linalg.norm(delta @ lattice, axis=2) < tol))
    if not np.all(match.any(axis=1)):
        raise AnalysisError(Msg("err.image_atom"))
    perm = np.argmax(match, axis=1)
    if len(np.unique(perm)) != len(positions):
        raise AnalysisError(Msg("err.image_atom"))
    return perm


def mech_rep(lattice, positions, numbers, rots, trans, symprec=1e-3):
    """力学表現 D(g) (操作数, 3N, 3N)。原子 i の変位を回転して写り先の原子 j へ。

    写り先の照合は Cartesian 距離で行い、許容差は対称性判定 (symprec) にそろえる。"""
    na, ns = len(numbers), len(rots)
    rcs = build_R_carts(lattice, rots)
    tol = max(2 * symprec, 1e-4)
    rep = np.zeros((ns, 3 * na, 3 * na), dtype=complex)
    for s in range(ns):
        perm = _atom_permutation(lattice, positions, numbers, rots[s], trans[s], tol)
        for i, j in enumerate(perm):
            rep[s, 3*j:3*j+3, 3*i:3*i+3] = rcs[s]
    return rep


def _match_irrep_info(irr, irr_conj, irreps_merged, labels, R_carts, tensors):
    """未統合の既約表現 (と複素ペア相手) に対応する統合済み表現の (ラベル, 活性, テンソル)。"""
    chi = get_character(irr)
    if irr_conj is not None:
        chi = chi + get_character(irr_conj)
    for k, merged in enumerate(irreps_merged):
        if np.max(np.abs(get_character(merged) - chi)) < 0.05:
            lbl = labels[k]
            return lbl, activity(merged, R_carts), tensors.get(lbl, [])
    return "?", (activity(irr, R_carts) if irr_conj is None else "?"), []


def isotypic_basis(rep, chi, dim):
    """既約表現 (指標 chi, 次元 dim) に属する変位の部分空間の実正規直交基底 (k, 3N)。

    射影演算子 P = (dim/h) Σ_g χ*(g) D(g) の像を求める。複素共役ペアは chi に
    ペアの指標の和 (実数) を渡せば、実空間での 2 次元モードの部分空間が得られる。
    結晶振動は実数なので、複素既約表現のモードは共役ペアで縮退し、実空間では
    1 つの 2 次元表現として現れる。"""
    h = len(rep)
    P = (dim / h) * sum(np.conj(chi[g]) * rep[g] for g in range(h))
    P = np.real(0.5 * (P + P.conj().T))
    w, v = np.linalg.eigh(P)
    return v[:, w > 0.5].T


def split_into_copies(basis, d, tiebreak):
    """既約表現の部分空間 (m·d 次元) を、対称操作で閉じた d 次元のコピー m 個に分ける。

    tiebreak は全対称操作と可換な一般の対称行列 (common.invariant_tiebreaker)。
    その固有空間は d 次元ずつ対称操作で閉じるので、固有値順に d 本ずつ束ねる。
    同じ既約表現が m>1 回現れるとき、各コピーの向きは対称性だけでは決まらない
    (任意性は部分空間の *中* だけで、部分空間そのものは厳密)。"""
    if basis.shape[0] <= d:
        return [basis]
    _, c = np.linalg.eigh(basis @ tiebreak @ basis.T)
    vecs = c.T @ basis
    return [vecs[k:k + d] for k in range(0, len(vecs), d)]


# ════════════════════════════════════════════════════════════════════════
# 7. 対称性適応変位基底 (因子群解析と QE 凍結フォノンで共通)
# ════════════════════════════════════════════════════════════════════════
def idealize_cell(lattice, positions, numbers, dataset, symprec=1e-3):
    """格子と座標を、spglib が見つけた対称操作で厳密に不変になるように理想化する。

    実在の CIF は格子定数や座標の丸めのため、対称操作で写した原子が 1e-4 Å 程度ずれる。
    そのまま使うと Cartesian の回転行列が直交行列からずれ、網目モデルの力定数が対称操作と
    可換でなくなって、m>1 の推測方向が既約表現のコピーの間で混ざる。
      格子: 計量 G = A Aᵀ を群で平均した G' = ⟨Wᵀ G W⟩ に、向きをほぼ保つ対称な変形
            A' = A M, M = (A⁻¹ G' A⁻ᵀ)^{1/2} で合わせる (A' A'ᵀ = G')
      座標: 各原子を、対称操作で写り合う原子の像 W⁻¹(x_j − t) の平均に置き換える
    返り値 (A', 座標', 原子の最大移動 [Å], 格子の最大ひずみ)。"""
    A = np.asarray(lattice, float)
    pos = np.asarray(positions, float)
    rots, trans = np.asarray(dataset.rotations), np.asarray(dataset.translations)
    G = A @ A.T
    G_sym = sum(W.T @ G @ W for W in rots) / len(rots)
    w, v = np.linalg.eigh(np.linalg.inv(A) @ G_sym @ np.linalg.inv(A).T)
    M = (v * np.sqrt(np.clip(w, 1e-12, None))) @ v.T        # 対称な平方根
    A_new = A @ M
    strain = float(np.abs(M - np.eye(3)).max())
    na = len(pos)
    tol = max(2 * symprec, 1e-4)
    acc = np.zeros_like(pos)
    for W, t in zip(rots, trans):
        W_inv = np.linalg.inv(W)
        perm = _atom_permutation(A, pos, numbers, W, t, tol)
        back = (pos[perm] - t) @ W_inv.T                  # 像を元に戻した位置
        acc += pos + ((back - pos + 0.5) % 1.0 - 0.5)
    pos_new = acc / len(rots)
    shift = float(np.max(np.linalg.norm((pos_new - pos) @ A, axis=1))) if na else 0.0
    return A_new, pos_new, shift, strain


def symmetry_adapted_blocks(lattice, positions, numbers, symprec=1e-3):
    """Γ 点の各既約表現について、対称性適応変位基底を作る。

    返り値 dict:
        lattice, positions, numbers : primitive cell (中心格子なら変換後)
        conv_natoms                 : 入力セルの原子数
        rots, trans, R_carts        : 対称操作
        rep                         : 力学表現 D(g) (実行列のリスト)
        dataset                     : spglib の対称性データ
        tiebreak                    : 全操作と可換な一般の対称行列 (コピーの分割用)
        blocks : list of dict(label, d, activity, tensor, chi_eff,
                               basis=(m·d, 3N) その既約表現の部分空間,
                               copies=[(d, 3N) 実配列 × m])"""
    lattice = np.array(lattice, float)
    positions = np.array(positions, float)
    numbers = [int(z) for z in numbers]
    conv_na = len(numbers)
    if not (np.isfinite(symprec) and symprec > 0):
        raise AnalysisError(Msg("err.symprec"))
    try:
        in_ds = spglib.get_symmetry_dataset((lattice, positions, numbers), symprec=symprec)
        if in_ds is None:
            raise AnalysisError(Msg("err.no_sg"))
        # 格子と座標を対称操作で厳密に不変になるよう理想化し、以降はその構造で計算する
        lattice, positions, shift, strain = idealize_cell(lattice, positions, numbers, in_ds, symprec)
        in_ds = spglib.get_symmetry_dataset((lattice, positions, numbers), symprec=symprec)
        # spgrep は conventional cell (標準の中心化) を前提にするので、入力が primitive cell
        # や菱面体軸・超格子でも、向きを変えずに (no_idealize) conventional cell に直して渡す
        conv = spglib.standardize_cell((lattice, positions, numbers), to_primitive=False,
                                       no_idealize=True, symprec=symprec)
        # Γ 点フォノンは primitive cell で定義される (中心格子の余分な並進を除く)
        prim = spglib.standardize_cell((lattice, positions, numbers), to_primitive=True,
                                       no_idealize=True, symprec=symprec)
        if conv is None or prim is None:
            raise AnalysisError(Msg("err.no_sg"))
        input_cell = lattice                      # 理想化した入力セル (表示・成分選択用)
        conv_lat = np.array(conv[0], float)
        irreps_conv, rots_conv, _, _ = get_spacegroup_irreps(
            conv_lat, np.array(conv[1], float), [int(z) for z in conv[2]], [0, 0, 0],
            symprec=symprec)
        lattice = np.array(prim[0], float)
        positions = np.array(prim[1], float)
        numbers = [int(z) for z in prim[2]]
        cell = (lattice, positions, numbers)
        dataset = spglib.get_symmetry_dataset(cell, symprec=symprec)
        sym = spglib.get_symmetry(cell, symprec=symprec)
    except AnalysisError:
        raise
    except Exception as e:                       # spglib / spgrep の内部エラー
        raise AnalysisError(Msg("err.no_sg_detail", err=e)) from e
    if dataset is None or sym is None or in_ds is None:
        raise AnalysisError(Msg("err.no_sg"))
    rots, trans = sym["rotations"], sym["translations"]
    R_carts = build_R_carts(lattice, rots)
    # ラベル・活性・テンソル成分は空間群の標準の向き (Bilbao と同じ x ∥ a, z ∥ c など) で
    # 決める。入力の向き (菱面体軸や非標準の設定) によらず同じ結果になる。
    std_rot = np.array(in_ds.std_rotation_matrix, float)
    R_std = np.array([std_rot @ R @ std_rot.T for R in R_carts])

    irreps_red = [reduce_irrep_to_pg(irr, rots_conv, conv_lat, rots, lattice)
                  for irr in irreps_conv]
    irreps = merge_complex_pairs(irreps_red)       # ラベル・活性用 (複素ペアを統合)
    labels = mulliken_labels(irreps, R_std)
    tensors = tensor_components(irreps, labels, R_std)
    rep = mech_rep(lattice, positions, numbers, rots, trans, symprec).real
    tiebreak = invariant_tiebreaker(list(rep))

    # 射影は未統合の既約表現の指標で行い、複素共役ペアは指標の和でまとめる
    chars = [get_character(irr) for irr in irreps_red]
    used = [False] * len(irreps_red)
    blocks = []
    for i, irr in enumerate(irreps_red):
        if used[i]:
            continue
        used[i] = True
        conj = None
        chi_eff = chars[i]
        if irr.shape[1] == 1 and not _is_real_character(chars[i]):
            j = _find_conjugate(chars, i, used)
            if j is not None:
                used[j] = True
                conj = irreps_red[j]
                chi_eff = chars[i] + chars[j]           # ペアの合成指標 (実)
        lbl, act, tcomp = _match_irrep_info(irr, conj, irreps, labels, R_std, tensors)
        d = 2 if (conj is not None or irr.shape[1] == 2) else irr.shape[1]
        basis = isotypic_basis(rep, chi_eff, irr.shape[1])
        if basis.shape[0]:
            blocks.append(dict(label=lbl, d=d, activity=act, tensor=list(tcomp),
                               chi_eff=chi_eff, basis=basis,
                               copies=split_into_copies(basis, d, tiebreak)))
    # 入力セルが標準の設定 (軸の取り方) か。違えばラベル・テンソル成分の軸を注記する
    std_setting = bool(np.allclose(in_ds.transformation_matrix, np.eye(3), atol=1e-6))
    return dict(lattice=lattice, positions=positions, numbers=numbers,
                conv_natoms=conv_na, rots=rots, trans=trans, R_carts=R_carts,
                R_std=R_std, std_rot=std_rot, std_setting=std_setting,
                dataset=dataset, blocks=blocks, rep=list(rep), tiebreak=tiebreak,
                input_cell=input_cell, idealized=(shift, strain))


# ════════════════════════════════════════════════════════════════════════
# 8. メイン解析
# ════════════════════════════════════════════════════════════════════════
def _translation_basis(na):
    """剛体並進 (音響) の正規直交基底 (3, 3N)。"""
    t = np.zeros((3, 3 * na))
    for ax in range(3):
        t[ax, ax::3] = 1.0 / np.sqrt(na)
    return t


def split_acoustic_optical(iso, d, n_ac, tb, mass3):
    """既約表現の部分空間 iso (m·d, 3N) を音響 (n_ac·d, 3N) と光学 (残り, 3N) に分ける。

    音響 = 部分空間のうち剛体並進 tb の成分 (この表現に含まれるベクトル表現の数 n_ac 組)。
    光学 = 残りのうち、質量で重みを付けて並進と直交する部分 (Σ M_i u_i = 0: 重心が動かない)。
    Γ 点の真の光学モードはこの条件を満たす。どちらも対称操作で閉じている。"""
    if not n_ac:
        return iso[:0], iso
    u, _, _ = np.linalg.svd(iso @ tb.T)                # 並進との重なりが大きい順
    acoustic = u[:, :n_ac * d].T @ iso
    c = (iso * mass3) @ acoustic.T                     # 質量計量での並進との内積 (k, n_ac·d)
    _, sv, vt = np.linalg.svd(c.T)
    rank = int(np.sum(sv > 1e-8 * max(sv.max(), 1e-12)))
    q, _ = np.linalg.qr((vt[rank:] @ iso).T)           # 表示・推測用に正規直交化 (空間は同じ)
    return acoustic, q.T


def analyze_crystal(lattice, positions, numbers, symbols=None, symprec=1e-3, notes=None):
    """結晶の Γ 点振動モードを解析して AnalysisResult を返す (ファイル出力はしない)。

    lattice   : 格子ベクトル (行ベクトル, Å)
    positions : 分率座標
    symprec   : 空間群判定の許容差 (spglib / spgrep 共通)
    notes     : 結果の注記の先頭に加える文 (読み込み時の注記など)"""
    sa = symmetry_adapted_blocks(lattice, positions, numbers, symprec)
    input_cell = sa["input_cell"]                  # 入力セル (対称性に合わせて理想化済み)
    lattice, positions, numbers = sa["lattice"], sa["positions"], sa["numbers"]
    ds = sa["dataset"]
    na = len(numbers)
    symbols = [element_symbol(z) for z in numbers]
    cart = positions @ lattice
    frac_in = (cart @ np.linalg.inv(input_cell)) % 1.0     # 入力セルの分率座標 (成分の選択用)
    frac_in = np.where(frac_in > 1 - 5e-4, 0.0, frac_in)
    masses = [atomic_mass(z) for z in numbers]
    wyckoffs = list(ds.wyckoffs)
    cubic = _is_cubic([_classify(R) for R in sa["R_carts"]])
    tb = _translation_basis(na)
    chi_vec = np.array([np.trace(R) for R in sa["R_carts"]])

    # サイト (対称操作で移り合う原子の組 = 軌道) の名前: 元素(Wyckoff)。
    # 同じ元素の軌道が複数あれば O1(e), O2(e) のように番号を付ける
    orbit_of = list(ds.equivalent_atoms)
    orbits = list(dict.fromkeys(orbit_of))                  # 出現順
    per_species = Counter(symbols[o] for o in orbits)
    numbering, site_name = Counter(), {}
    for o in orbits:
        s = symbols[o]
        numbering[s] += 1
        num = str(numbering[s]) if per_species[s] > 1 else ""
        site_name[o] = f"{s}{num}({wyckoffs[o]})"

    notes = list(notes or [])
    shift, strain = sa["idealized"]
    if shift > 1e-5 or strain > 1e-6:
        notes.append(Msg("cry.idealized", shift=shift, strain=strain))
    if sa["conv_natoms"] > na:
        notes.append(Msg("cry.centered", n0=sa["conv_natoms"], n=na))
    has_inv = any(not p and o == 1 for o, _, p in (_classify(R) for R in sa["R_carts"]))
    if has_inv:
        notes.append(Msg("cry.inversion"))
    if not sa["std_setting"]:
        notes.append(Msg("cry.std_axes", sg=ds.international))
    mass3 = np.repeat(np.asarray(masses, float), 3)

    modes, decomposition = [], []
    hessians = None
    for blk in sa["blocks"]:
        iso, d = blk["basis"], blk["d"]               # この既約表現の部分空間 (m·d, 3N)
        m = iso.shape[0] // d
        decomposition.append(f"{m}{blk['label']}" if m > 1 else blk["label"])

        n_ac = int(round(abs(_mult(blk["chi_eff"], chi_vec))))   # 含まれるベクトル表現の数
        acoustic_rows, optical_rows = split_acoustic_optical(iso, d, n_ac, tb, mass3)

        # 光学コピーが複数あれば、その部分空間内の向きを網目モデルで推測
        m_opt = optical_rows.shape[0] // d
        refined = m_opt > 1
        if refined:
            if hessians is None:
                hessians = (network_hessian(cart, lattice, cutoff=NETWORK_CUTOFF),
                            network_hessian(cart, lattice, cutoff=NETWORK_CUTOFF, isotropic=True))
            optical_rows = refine_in_subspace(optical_rows, cart, masses, lattice=lattice,
                                              cutoff=NETWORK_CUTOFF, mass_weighted=False,
                                              sym_mats=sa["rep"], hessians=hessians,
                                              tiebreaker=sa["tiebreak"])
        copies = ([(acoustic_rows[k:k + d], True) for k in range(0, len(acoustic_rows), d)]
                  + [(optical_rows[k:k + d], False) for k in range(0, len(optical_rows), d)])

        # 関与するサイト: この既約表現の部分空間で動きうる原子の軌道 (対称性だけで決まる厳密な情報)
        moving = [np.sum(iso[:, 3*a:3*a+3] ** 2) > 1e-8 for a in range(na)]
        sites = [site_name[o] for o in orbits
                 if any(moving[a] for a in range(na) if orbit_of[a] == o)]
        for basis, is_ac in copies:
            # 縮退成分は結晶軸 a, b, c を基準に表示用に選び直す (全物質で同じ規則)
            disps = axis_aligned_components([row.reshape(na, 3) for row in basis], input_cell,
                                            sort_keys=frac_in)
            raman = (raman_tensor_forms([u.ravel() for u in disps], sa["rep"], sa["R_std"])
                     if "Raman" in blk["activity"] and not is_ac else [])
            modes.append(Mode(
                label=blk["label"], degeneracy=d, activity=blk["activity"],
                tensor=blk["tensor"], displacements=disps,
                refined=refined and not is_ac, acoustic=is_ac,
                character=_displacement_character(
                    np.array([(row.reshape(na, 3) @ sa["std_rot"].T).ravel() for row in basis]),
                    cubic),
                contributions=sites, raman_tensors=raman))

    sg = f"{ds.international} (No. {ds.number})"
    result = AnalysisResult(kind="crystal", numbers=numbers, symbols=symbols, cart=cart,
                            modes=modes, point_group=ds.pointgroup, space_group=sg,
                            lattice=lattice, cell=input_cell,
                            decomposition="Γ = " + " + ".join(decomposition), notes=notes)
    result.meta = dict(wyckoffs=wyckoffs, n_ops=len(sa["rots"]), cubic=cubic)
    return result


def _displacement_character(basis, cubic):
    """変位が z (面外) / xy (面内) / 混合 のどれに偏っているかのキー ("z", "xy", "mixed")。
    立方晶では区別しない ("iso")。表示は Mode.character_label で現在の言語に。"""
    if cubic:
        return "iso"
    z = np.sum(basis[:, 2::3] ** 2)
    xy = np.sum(basis[:, 0::3] ** 2) + np.sum(basis[:, 1::3] ** 2)
    total = z + xy
    if total > 0 and z / total > 0.8:
        return "z"
    if total > 0 and xy / total > 0.8:
        return "xy"
    return "mixed"


# ════════════════════════════════════════════════════════════════════════
# 9. 出力 (テキストレポート / XSF)
# ════════════════════════════════════════════════════════════════════════
PORTO = [("Z(XX)Z̄", "xx"), ("Z(YY)Z̄", "yy"), ("Z(XY)Z̄", "xy"),
         ("X(YY)X̄", "yy"), ("X(YZ)X̄", "yz"), ("X(ZZ)X̄", "zz"),
         ("Y(XZ)Ȳ", "xz"), ("Y(XX)Ȳ", "xx"), ("Y(ZZ)Ȳ", "zz")]


def format_report(result: AnalysisResult, wyckoffs=None, n_ops=None, cubic=False) -> str:
    """CLI に表示するテキストレポートを現在の言語で作る。"""
    na = len(result.numbers)
    L = ["=" * 72, t("rep.cry.title"), "=" * 72, "",
         t("rep.cry.structure"), t("rep.cry.sg", v=result.space_group),
         t("rep.cry.pg", v=result.point_group),
         t("rep.cry.natoms", n=na) + (t("rep.cry.nops", n=n_ops) if n_ops else "")]
    if wyckoffs is not None:
        seen = []
        for s, w in zip(result.symbols, wyckoffs):
            if (s, w) not in seen:
                seen.append((s, w))
                n = sum(1 for s2, w2 in zip(result.symbols, wyckoffs) if (s2, w2) == (s, w))
                L.append(f"    {s} ×{n}  (Wyckoff {w})")
    for note in result.notes:
        L.append(t("rep.note", v=note))

    modes = result.modes
    opt = result.optical_modes
    aco = [m for m in modes if m.acoustic]
    L += ["", t("rep.cry.decomp", n=3 * na), "", f"  {result.decomposition}", "",
          t("rep.cry.acoustic", n=len(aco), dof=sum(m.degeneracy for m in aco)),
          t("rep.cry.optical", n=len(opt), dof=sum(m.degeneracy for m in opt))]
    for a in ACTIVITIES:
        sub = [m for m in opt if m.activity == a]
        if sub:
            L.append(t("rep.activity_line", sym=ACTIVITY_SYMBOL[a], act=a, n=len(sub),
                       dof=sum(m.degeneracy for m in sub)))
    widths = (10, 3, 10, 18, 5, 17)
    cols = t("rep.cry.cols").split("|")
    L += ["", "━" * 76,
          " " + "".join(pad(c, w) for c, w in zip(cols, widths)) + cols[-1], "━" * 76]
    for m in modes:
        atoms = ", ".join(m.contributions)
        cells = (m.label, m.degeneracy, m.activity, m.character_label,
                 "♪" if m.acoustic else "", m.direction_label)
        L.append(" " + "".join(pad(c, w) for c, w in zip(cells, widths)) + atoms)
    L += ["━" * 76] + t("rep.cry.legend").splitlines() + [""]

    raman = [m for m in opt if "Raman" in m.activity]
    if raman:
        L += ["━" * 76, t("rep.cry.raman_title"), "━" * 76, t("rep.cry.raman_note")]
        seen = set()
        for m in raman:
            if m.label in seen:
                continue
            seen.add(m.label)
            L += ["", t("rep.cry.raman_mode", label=m.label, d=m.degeneracy),
                  t("rep.cry.tensor_comps", v=", ".join(m.tensor))]
            if len(m.raman_tensors) == 1:
                L.append(t("rep.cry.tensor_form", v=m.raman_tensors[0]))
            else:
                for k, form in enumerate(m.raman_tensors):
                    L.append(t("rep.cry.tensor_form_k", k=k + 1, v=form))
            L += ["", t("rep.cry.porto")]
            for name, comp in PORTO:
                L.append(f"      {'✓' if comp in m.tensor else '✗'}  {name:12}  ← {comp}")
        L += ["", "━" * 76]
    return "\n".join(L)


def cell_atoms(result: AnalysisResult, boundary=True, tol=1e-4):
    """入力ファイルの単位胞 (result.cell) の中にある原子をすべて列挙する。

    解析は primitive cell で行うが、表示と VESTA 出力は入力の単位胞で行う。
    Γ 点の振動は全ての格子並進で同位相なので、像の原子は元の原子と同じ変位をする。
    boundary=True なら単位胞の境界 (分率座標 0 と 1) 上の像も含める。
    返り値: (Cartesian 座標 (M, 3), 分率座標 (M, 3), 元の原子の番号 (M,))"""
    cell = result.cell if result.cell is not None else result.lattice
    prim = result.lattice
    inv_cell = np.linalg.inv(cell)
    # 単位胞の 8 隅を primitive の分率座標で表し、必要な並進の範囲を決める
    corners = np.array([[i, j, k] for i in (0, 1) for j in (0, 1) for k in (0, 1)]) @ cell
    reach = np.abs(corners @ np.linalg.inv(prim)).max(axis=0).astype(int) + 2
    shifts = np.array([[i, j, k] for i in range(-reach[0], reach[0] + 1)
                       for j in range(-reach[1], reach[1] + 1)
                       for k in range(-reach[2], reach[2] + 1)]) @ prim
    upper = 1 + tol if boundary else 1 - tol
    carts, fracs, parents, seen = [], [], [], set()
    for a, r in enumerate(result.cart):
        cand = r + shifts
        f = cand @ inv_cell
        keep = np.all((f > -tol) & (f < upper), axis=1)
        for rc, fc in zip(cand[keep], f[keep]):
            key = tuple(np.round(fc + 1e-9, 3))
            if key in seen:
                continue
            seen.add(key)
            carts.append(rc)
            fracs.append(np.where(np.abs(fc) < tol, 0.0, fc))
            parents.append(a)
    return np.array(carts), np.array(fracs), np.array(parents, int)


def write_mode_files(result: AnalysisResult, outdir, scale=0.2, formats=("xsf", "vesta")):
    """光学モード (縮退成分ごと) を XSF と VESTA 形式で書き出し、書いたパスのリストを返す。

    XSF   : primitive cell。矢印の見た目は VESTA 任せ (長さは最大 scale Å に縮める)。
    VESTA : 入力の単位胞。灰色の細い矢印 (中点が原子の中心、最大で中心から ±0.67 Å) と
            結合の設定 (VESTA 既定の結合探索表) まで書き込む。
    formats で書き出す形式を選べる (("xsf",), ("vesta",), 両方)。
    ファイル番号は解析結果の全モード (音響を含む) での通し番号。"""
    ensure_outdir(outdir)
    paths = []
    _, frac_cell, parent = cell_atoms(result, boundary=False)
    cell = result.cell if result.cell is not None else result.lattice
    labels = [result.symbols[p] for p in parent]
    numbers_cell = [result.numbers[p] for p in parent]
    for i, m in enumerate(result.modes, 1):
        if m.acoustic:
            continue
        for c, disp in enumerate(m.displacements):
            stem = os.path.join(outdir, mode_file_stem(i, m, c))
            header = xsf_header(i, m, c, "true eigenvector requires DFT (Stage 2)")
            if "xsf" in formats:
                write_xsf(stem + ".xsf", result.lattice, result.cart, result.numbers,
                          normalize_max(disp), header, scale=scale)
                paths.append(stem + ".xsf")
            if "vesta" in formats:
                write_vesta(stem + ".vesta", cell, frac_cell, numbers_cell, labels,
                            np.asarray(disp)[parent], title=" / ".join(header[:3]),
                            arrow_length=0.67)
                paths.append(stem + ".vesta")
    return paths
