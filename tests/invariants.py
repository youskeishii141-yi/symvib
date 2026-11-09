"""解析結果を、コードの内部経路に頼らずに検査する「不変量」と外部基準との照合。

結晶 (analyze_crystal の結果) について
  - 各モードの変位が対称操作で閉じた既約表現 (|χ|²/h = 1; 複素ペア統合は 2) の部分空間
  - 同じラベルは同じ指標、違うラベルは直交する指標
  - 自由度の合計 3N、音響 3、光学モードは Σ M u = 0
  - IR / Raman 活性とラマンテンソル成分を、ベクトル表現・対称 2 階テンソル表現への射影から
    独立に求めて一致
  - Bilbao Crystallographic Server の点群指標表 (tests/data/bilbao/point_tables.json; 配布物には
    含めず tests/fetch_bilbao_tables.py で取得。未取得なら照合を省略) と
    Mulliken ラベル・活性・テンソル成分が一致 (操作は標準設定の軸方向で Bilbao のクラスに対応づける)
  - 相関法 (各操作で不動な原子の数 × tr R) による既約分解の重複度が一致
分子 (analyze_molecule の結果) についても同種の検査を行う。

対称操作は spglib (結晶) / find_point_group_operations (分子) で取り直し、力学表現も
ここで作り直す。解析コードの指標表・既約表現・ラベル付けの関数は一切使わない。
"""
from __future__ import annotations

import json
import os
import re
from collections import Counter

import numpy as np
import spglib

from symvib.common import atomic_mass

_BILBAO = None
COMPONENTS = ["xx", "yy", "zz", "xy", "xz", "yz"]
_PAIRS = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]


BILBAO_PATH = os.path.join(os.path.dirname(__file__), "data", "bilbao", "point_tables.json")


def bilbao_available():
    """Bilbao の指標表を取得済みか (配布物には含めない。tests/fetch_bilbao_tables.py で取得)。"""
    return os.path.exists(BILBAO_PATH)


def bilbao_tables():
    """Bilbao の 32 点群の指標表 (tests/fetch_bilbao_tables.py で取得したもの)。未取得なら None。"""
    global _BILBAO
    if _BILBAO is None and bilbao_available():
        with open(BILBAO_PATH, encoding="utf-8") as f:
            _BILBAO = json.load(f)["tables"]
    return _BILBAO


# ────────────────────────────────────────────────────────────────────────
# 幾何・表現の道具
# ────────────────────────────────────────────────────────────────────────
def cart_rotations(lattice, rots):
    """分率座標の回転行列 → Cartesian (lattice は行ベクトル)。"""
    li = np.linalg.inv(lattice)
    return np.array([lattice.T @ r @ li.T for r in rots])


def classify(R, tol=1e-3):
    """Cartesian 直交行列 → (真回転か, 次数 n, 軸の単位ベクトル (符号任意))。

    反転 (1, 偽)、鏡映は (2, 偽) で軸 = 法線、回反 -n は (n, 偽) (−R が n 回回転)。"""
    proper = np.linalg.det(R) > 0
    rp = R if proper else -R
    cos = np.clip((np.trace(rp) - 1) / 2, -1, 1)
    ang = np.arccos(cos)
    n = 1 if ang < 1e-6 else int(round(2 * np.pi / ang))
    if n == 1:
        return proper, 1, None
    w, v = np.linalg.eig(rp)
    axis = np.real(v[:, np.argmin(np.abs(w - 1))])
    return proper, n, axis / np.linalg.norm(axis)


def sym2_matrix(R):
    """対称 2 階テンソル (xx, yy, zz, xy, xz, yz の 6 次元) 上の R の表現行列。"""
    # 正規直交基底 (非対角成分は (e_ij+e_ji)/√2) のテンソルを R で変換して成分を読む
    basis = []
    for (i, j) in _PAIRS:
        T = np.zeros((3, 3))
        if i == j:
            T[i, i] = 1.0
        else:
            T[i, j] = T[j, i] = 1 / np.sqrt(2)
        basis.append(T)
    return np.array([[np.sum(basis[a] * (R @ basis[b] @ R.T)) for b in range(6)]
                     for a in range(6)])


def projected_components(chi, R_list, norm=None):
    """指標 chi の既約表現に属する対称テンソル成分 (xx..yz) の集合と、ベクトル表現の重複度。

    P = (d/h) Σ χ*(g) T(g) の対角要素が 0 でない成分を採る (閾値 1e-6)。"""
    h = len(R_list)
    chi = np.asarray(chi, complex)
    n2 = norm if norm is not None else float(np.sum(np.abs(chi) ** 2).real)
    T = sum(np.conj(c) * sym2_matrix(R) for c, R in zip(chi, R_list)) / n2
    V = sum(np.conj(c) * R for c, R in zip(chi, R_list)) / n2
    comps = {COMPONENTS[a] for a in range(6) if abs(T[a, a]) > 1e-6}
    return comps, float(np.trace(V).real), float(np.trace(T).real)


def mech_rep_from_ops(lattice, frac, numbers, rots, trans, tol=1e-3):
    """力学表現 D(g) と、各操作で不動な原子の数 (相関法用) を返す。"""
    na = len(numbers)
    R = cart_rotations(lattice, rots)
    D, fixed = [], []
    for s in range(len(rots)):
        M = np.zeros((3 * na, 3 * na))
        nfix = 0
        for i in range(na):
            img = rots[s] @ frac[i] + trans[s]
            for j in range(na):
                d = img - frac[j]
                d -= np.round(d)
                if numbers[i] == numbers[j] and np.linalg.norm(d @ lattice) < tol:
                    M[3 * j:3 * j + 3, 3 * i:3 * i + 3] = R[s]
                    nfix += (i == j)
                    break
            else:
                raise ValueError("image atom not found")
        D.append(M)
        fixed.append(nfix)
    return np.array(D), R, np.array(fixed)


# ────────────────────────────────────────────────────────────────────────
# Bilbao の表との対応づけ
# ────────────────────────────────────────────────────────────────────────
def _parse_class(symbol):
    """Bilbao のクラス記号 → (真回転か, 次数, 方向 [uvw] or None)。'4+' の ± は無視。"""
    improper = symbol.startswith("-") or symbol.startswith("m")
    s = symbol[1:] if symbol.startswith("-") else symbol
    if s.endswith(("+", "-")):                      # 回転の向き (4+, -3- など) は無視
        s = s[:-1]
    if s.startswith("m"):
        order, rest = 2, s[1:]
    else:
        m = re.match(r"(\d)(.*)", s)
        order, rest = int(m.group(1)), m.group(2)
    direction = None
    if rest in ("x", "y", "z"):
        direction = {"x": (1, 0, 0), "y": (0, 1, 0), "z": (0, 0, 1)}[rest]
    elif rest:
        direction = tuple(int(x) for x in re.findall(r"-?\d", rest))
    return (not improper), order, direction


def _label_norm(label):
    """Bilbao の 1E / 2E を外し、プライムの位置をそろえる。"""
    return re.sub(r"^[12](?=E)", "", label)


def merged_bilbao_rows(table):
    """複素共役ペア (1E, 2E) を 1 つの実表現にまとめた行 [(label, chars(実), functions)]。"""
    rows = {}
    for ir in table["irreps"]:
        lbl = _label_norm(ir["label"])
        ch = np.array([complex(a, b) for a, b in ir["chars"]])
        if lbl in rows:
            rows[lbl][0] += ch
            rows[lbl][1] = rows[lbl][1] or ir["functions"]
        else:
            rows[lbl] = [ch, ir["functions"]]
    return [(lbl, np.real_if_close(ch, tol=1e-6), fn) for lbl, (ch, fn) in rows.items()]


# Bilbao の点群表は主軸を z にとる。単斜晶の標準設定 (unique axis b) では 2 回軸が y なので、
# 表の (x, y, z) を標準設定の (x, z, y) に読み替える
_MONOCLINIC_FRAME = {"x": "x", "y": "z", "z": "y"}


def functions_to_expected(fn, frame=None):
    """Bilbao の functions 欄 → (IR 活性か, Raman テンソル成分の集合)。frame で軸名を読み替える。"""
    frame = frame or {}
    ir = bool(re.search(r"(?<![a-zA-Z0-9_(])(x|y|z)(?![a-zA-Z0-9^])", fn))
    text = fn.replace("^", "")
    comps = set()
    for tok in re.findall(r"[xyz]2|[xyz]{2}", text):
        letters = [tok[0], tok[0]] if tok.endswith("2") else [tok[0], tok[1]]
        letters = [frame.get(c, c) for c in letters]
        comps.add("".join(sorted(letters, key="xyz".index)))
    return ir, comps


def assign_bilbao_classes(R_std, std_lattice, table):
    """標準設定の Cartesian 操作 R_std を Bilbao の表のクラス番号に対応づける。

    返り値 (クラス番号のリスト, 注記)。方向で区別できないときは注記に残す。"""
    classes = [_parse_class(c["symbol"]) for c in table["classes"]]
    A = np.asarray(std_lattice, float)
    R_std = np.asarray(R_std)
    h = len(R_std)

    def index(M):
        return int(np.argmin([np.abs(M - X).max() for X in R_std]))

    # 共役類 (数値的に)
    seen, conj = set(), []
    for a in range(h):
        if a in seen:
            continue
        cl = sorted({index(R_std[g] @ R_std[a] @ R_std[g].T) for g in range(h)})
        seen.update(cl)
        conj.append(cl)

    def direction_score(k, members):
        """Bilbao のクラス k の方向と、共役類 members の軸との一致 (最大の |cos|)。"""
        d = classes[k][2]
        if d is None:
            return 0.5
        v = np.asarray(d, float) @ A
        v /= np.linalg.norm(v)
        best = 0.0
        for a in members:
            axis = classify(R_std[a])[2]
            if axis is not None:
                best = max(best, abs(float(np.dot(v, axis))))
        return best

    col, notes, used = [None] * h, [], []
    for cl in conj:
        proper, n, _ = classify(R_std[cl[0]])
        cands = [k for k, (p, o, _) in enumerate(classes) if p == proper and o == n]
        if not cands:
            raise ValueError(f"no Bilbao class for operation (proper={proper}, n={n})")
        if len(cands) > 1:
            scored = sorted(cands, key=lambda k: (direction_score(k, cl), k not in used),
                            reverse=True)
            best = scored[0]
            if direction_score(best, cl) < 0.999 and any(classes[k][2] is not None for k in cands):
                # 方向が Bilbao の表の向きと違う設定 (例: -4m2 の 2 回軸は [110])。
                # 同じ型・同じ大きさのクラスがまだ残っていればそれに割り当てる
                free = [k for k in cands if k not in used and table["classes"][k]["mult"] == len(cl)]
                best = free[0] if free else best
                notes.append(f"{'' if proper else '-'}{n} axes not along the Bilbao table "
                             f"directions (assigned to class {table['classes'][best]['symbol']})")
            cands = [best]
        used.append(cands[0])
        for a in cl:
            col[a] = cands[0]
    counts = Counter(col)
    for k, c in enumerate(table["classes"]):
        if counts.get(k, 0) != c["mult"]:
            notes.append(f"class {c['symbol']} has {counts.get(k, 0)} ops (Bilbao: {c['mult']})")
    return col, notes


def bilbao_label(chi_std, col, table):
    """標準設定での指標 chi_std (操作順) に一致する Bilbao の行 → (ラベル, functions)。"""
    chi = np.asarray(chi_std)
    for lbl, ch, fn in merged_bilbao_rows(table):
        if np.allclose(np.asarray([ch[c] for c in col]), chi, atol=1e-4):
            return lbl, fn
    return None, None


# ────────────────────────────────────────────────────────────────────────
# 結晶の検査
# ────────────────────────────────────────────────────────────────────────
def check_crystal(result, lattice_in, frac_in, numbers_in, symprec=1e-3, bilbao=True):
    """analyze_crystal の結果 result を独立に検査し、問題の文字列リストを返す (空なら合格)。

    lattice_in / frac_in / numbers_in は解析に渡した入力 (標準設定への回転を spglib で取り直す)。
    返り値 (問題のリスト, 情報 dict: bilbao_labels={code label: Bilbao label}, ...)。"""
    probs, info = [], {}
    numbers = list(result.numbers)
    na = len(numbers)
    lattice = np.asarray(result.lattice, float)
    frac = np.linalg.solve(lattice.T, np.asarray(result.cart).T).T
    masses = np.array([atomic_mass(z) for z in numbers])

    # 自由度・音響
    if sum(m.degeneracy for m in result.modes) != 3 * na:
        probs.append(f"DOF {sum(m.degeneracy for m in result.modes)} != 3N={3 * na}")
    if sum(m.degeneracy for m in result.modes if m.acoustic) != 3:
        probs.append("acoustic DOF != 3")
    # 光学モード: 重心不動
    for m in result.modes:
        if not m.acoustic:
            for u in m.displacements:
                if np.linalg.norm((masses[:, None] * u).sum(0)) > 1e-6 * masses.sum():
                    probs.append(f"{m.label}: centre of mass moves")
                    break

    # 対称操作 (spglib で取り直す) と力学表現
    sym = spglib.get_symmetry((lattice, frac, numbers), symprec=symprec)
    rots, trans = sym["rotations"], sym["translations"]
    D, R, fixed = mech_rep_from_ops(lattice, frac, numbers, rots, trans, max(2 * symprec, 1e-4))
    h = len(rots)
    ds_in = spglib.get_symmetry_dataset((np.asarray(lattice_in, float), np.asarray(frac_in, float),
                                         [int(z) for z in numbers_in]), symprec=symprec)
    std_rot = np.asarray(ds_in.std_rotation_matrix, float)
    R_std = np.array([std_rot @ r @ std_rot.T for r in R])
    chi_vec = np.array([np.trace(r) for r in R])
    chi_sym2 = np.array([np.trace(sym2_matrix(r)) for r in R])

    by_label = {}
    for m in result.modes:
        B = np.array([np.asarray(u).reshape(-1) for u in m.displacements])
        q, _ = np.linalg.qr(B.T)
        if q.shape[1] != m.degeneracy:
            probs.append(f"{m.label}: components are linearly dependent")
            continue
        # 対称操作で閉じているか
        for Dg in D:
            moved = Dg @ q
            if not np.allclose(q @ (q.T @ moved), moved, atol=1e-6):
                probs.append(f"{m.label}: not closed under symmetry")
                break
        chi = np.array([np.trace(q.T @ Dg @ q) for Dg in D])
        norm = chi @ chi / h
        if not (abs(norm - 1) < 1e-4 or abs(norm - 2) < 1e-4):
            probs.append(f"{m.label}: reducible (|chi|^2/h = {norm:.3f})")
        if m.label in by_label:
            if not np.allclose(by_label[m.label][0], chi, atol=1e-4):
                probs.append(f"{m.label}: inconsistent characters between modes")
        else:
            by_label[m.label] = (chi, norm)
        # 活性 (独立計算)
        ir = abs(chi @ chi_vec) > 1e-4
        ra = abs(chi @ chi_sym2) > 1e-4
        act = "Raman+IR" if ir and ra else "Raman" if ra else "IR" if ir else "Silent"
        if act != m.activity:
            probs.append(f"{m.label}: activity {m.activity} (independent: {act})")
        # テンソル成分 (標準設定の向きで)
        comps, _, _ = projected_components(chi, R_std, norm=norm * h)
        if comps != set(m.tensor) and ra:
            probs.append(f"{m.label}: tensor {sorted(m.tensor)} (independent: {sorted(comps)})")
        if not ra and m.tensor:
            probs.append(f"{m.label}: tensor given but not Raman active")
    labels = list(by_label)
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            if abs(by_label[labels[i]][0] @ by_label[labels[j]][0]) > 1e-4:
                probs.append(f"{labels[i]} / {labels[j]}: same irrep with different labels")

    # 相関法: 既約表現の重複度 n = (1/h) Σ χ*(g) χ_mech(g), χ_mech(g) = (不動原子数) tr R
    chi_mech = fixed * chi_vec
    counts = Counter(m.label for m in result.modes)
    for lbl, (chi, norm) in by_label.items():
        n = (chi @ chi_mech) / (norm * h)
        if abs(n - counts[lbl]) > 1e-6:
            probs.append(f"{lbl}: multiplicity {counts[lbl]} (correlation method: {n:.3f})")
    if abs(sum(by_label[l][0][0] * counts[l] for l in by_label) - 3 * na) > 1e-6:
        probs.append("sum of d·n != 3N")

    # Bilbao の指標表との照合 (表を取得していなければ省略し、その旨を info に残す)
    if bilbao and not bilbao_available():
        info["bilbao_skipped"] = True
    elif bilbao:
        pg = ds_in.pointgroup.strip()
        table = bilbao_tables().get(pg)
        if table is None:
            probs.append(f"no Bilbao table for point group {pg}")
        else:
            try:
                col, notes = assign_bilbao_classes(R_std, ds_in.std_lattice, table)
            except ValueError as e:
                probs.append(str(e))
                col, notes = None, []
            info["bilbao_notes"] = notes
            if col is not None:
                bl = {}
                for lbl, (chi, norm) in by_label.items():
                    blabel, fn = bilbao_label(chi, col, table)
                    bl[lbl] = blabel
                    if blabel is None:
                        probs.append(f"{lbl}: characters match no Bilbao irrep")
                        continue
                    if blabel != lbl:
                        probs.append(f"{lbl}: Bilbao label is {blabel}")
                    if fn is None or pg == "1":            # C1 の表には functions 欄がない
                        continue
                    frame = _MONOCLINIC_FRAME if 3 <= ds_in.number <= 15 else None
                    ir_b, comps_b = functions_to_expected(fn, frame)
                    mode = next(m for m in result.modes if m.label == lbl)
                    if ir_b != ("IR" in mode.activity):
                        probs.append(f"{lbl}: IR activity differs from Bilbao functions ({fn})")
                    if bool(comps_b) != ("Raman" in mode.activity):
                        probs.append(f"{lbl}: Raman activity differs from Bilbao functions ({fn})")
                    if comps_b and comps_b != set(mode.tensor):
                        probs.append(f"{lbl}: tensor {sorted(mode.tensor)} vs Bilbao {sorted(comps_b)}")
                info["bilbao_labels"] = bl
    return probs, info


# ────────────────────────────────────────────────────────────────────────
# 分子の検査
# ────────────────────────────────────────────────────────────────────────
def check_molecule(result):
    """analyze_molecule の結果を独立に検査し、問題の文字列リストを返す。"""
    from symvib.molecule import find_point_group_operations
    probs = []
    numbers = list(result.numbers)
    na = len(numbers)
    n_vib = 3 * na - (5 if result.is_linear else 6)
    if sum(m.degeneracy for m in result.modes) != n_vib:
        probs.append(f"DOF {sum(m.degeneracy for m in result.modes)} != {n_vib}")
    masses = np.array([atomic_mass(z) for z in numbers])
    sm = np.sqrt(masses)
    cart = np.asarray(result.cart)
    com = (masses[:, None] * cart).sum(0) / masses.sum()
    for m in result.modes:
        for u in m.displacements:
            if np.linalg.norm((masses[:, None] * u).sum(0)) > 1e-6 * masses.sum():
                probs.append(f"{m.label}: momentum not zero")
            L = np.cross(cart - com, u) * masses[:, None]
            if np.linalg.norm(L.sum(0)) > 1e-6 * masses.sum():
                probs.append(f"{m.label}: angular momentum not zero")
    if result.is_linear:
        return probs
    ops = find_point_group_operations(cart, numbers, 1e-2)
    h = len(ops)
    D = []
    for R, perm in ops:
        M = np.zeros((3 * na, 3 * na))
        for i, k in enumerate(perm):
            M[3 * k:3 * k + 3, 3 * i:3 * i + 3] = R
        D.append(M)
    Rs = [R for R, _ in ops]
    chi_vec = np.array([np.trace(R) for R in Rs])
    chi_sym2 = np.array([np.trace(sym2_matrix(R)) for R in Rs])
    by_label = {}
    for m in result.modes:
        B = np.array([(u * sm[:, None]).reshape(-1) for u in m.displacements])
        q, _ = np.linalg.qr(B.T)
        for Dg in D:
            moved = Dg @ q
            if not np.allclose(q @ (q.T @ moved), moved, atol=1e-6):
                probs.append(f"{m.label}: not closed under symmetry")
                break
        chi = np.array([np.trace(q.T @ Dg @ q) for Dg in D])
        norm = chi @ chi / h
        if not (abs(norm - 1) < 1e-4 or abs(norm - 2) < 1e-4):
            probs.append(f"{m.label}: reducible (|chi|^2/h = {norm:.3f})")
        if m.label in by_label and not np.allclose(by_label[m.label], chi, atol=1e-4):
            probs.append(f"{m.label}: inconsistent characters")
        by_label.setdefault(m.label, chi)
        ir, ra = abs(chi @ chi_vec) > 1e-4, abs(chi @ chi_sym2) > 1e-4
        act = "Raman+IR" if ir and ra else "Raman" if ra else "IR" if ir else "Silent"
        if act != m.activity:
            probs.append(f"{m.label}: activity {m.activity} (independent: {act})")
    labels = list(by_label)
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            if abs(by_label[labels[i]] @ by_label[labels[j]]) > 1e-4:
                probs.append(f"{labels[i]} / {labels[j]}: same irrep with different labels")
    return probs
