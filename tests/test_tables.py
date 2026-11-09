"""分子の指標表 (molecule_tables) を群論的に検証する。

各点群について理想的な操作の群を生成元から作り、
  - 複素指標表 (COMPLEX_TABLES) も含めて、表の各行が正則表現上で射影演算子 (P² = P, tr P = d²)
    になること (クラスの対応は解析コードと同じ fine キー / 共役類割り当てで行う)
  - linear (IR) / quadratic (Raman) の割り当てが、ベクトル表現・対称 2 階テンソル表現を
    表の指標で分解した結果と一致すること。記載された成分 (xz, yz, xx-yy, 2zz-xx-yy, …) も
    射影演算子で得られる成分と一致すること
を確かめる。実表現の表 (CHAR_TABLES) の射影検査は test_point_groups.py にある。
"""
import math
import re

import numpy as np
import pytest

from structures import _closure, _rot
from symvib import molecule as M
from symvib.molecule_tables import CHAR_TABLES, COMPLEX_TABLES

SH = np.diag([1., 1., -1.])
Z = (0, 0, 1)
PHI = (1 + 5 ** 0.5) / 2
SPECIAL = {
    "Td": [_rot((1, 1, 1), 2 * math.pi / 3), _rot(Z, math.pi / 2) @ SH],
    "Oh": [_rot((1, 1, 1), 2 * math.pi / 3), _rot(Z, math.pi / 2), -np.eye(3)],
    "O": [_rot((1, 1, 1), 2 * math.pi / 3), _rot(Z, math.pi / 2)],
    "T": [_rot((1, 1, 1), 2 * math.pi / 3), _rot(Z, math.pi)],
    "Th": [_rot((1, 1, 1), 2 * math.pi / 3), _rot(Z, math.pi), -np.eye(3)],
    "I": [_rot((0, 1, PHI), 2 * math.pi / 5), _rot(Z, math.pi)],
    "Ih": [_rot((0, 1, PHI), 2 * math.pi / 5), _rot(Z, math.pi), -np.eye(3)],
    "Cs": [SH], "Ci": [-np.eye(3)], "C1": [np.eye(3)],
}


def generators(name):
    if name in SPECIAL:
        return SPECIAL[name]
    kind, n, suf = re.fullmatch(r"([CDS])(\d+)([vhd]?)", name).groups()
    n = int(n)
    if kind == "S":
        return [_rot(Z, 2 * math.pi / n) @ SH]
    g = [_rot(Z, 2 * math.pi / n)]
    if kind == "D":
        g.append(np.diag([1., -1., -1.]))                   # C2' along x
    if suf == "v":
        g.append(np.diag([1., -1., 1.]))                    # σv: xz plane (normal y)
    if suf == "h":
        g.append(SH)
    if suf == "d":
        g.append(_rot(Z, math.pi / n) @ SH)
    return g


def _components(expr):
    """'xx-yy', '2zz-xx-yy', 'xx+yy', 'xz' → 含まれる成分の集合 {'xx', 'yy'}。"""
    return {"".join(sorted(tok, key="xyz".index)) for tok in re.findall(r"[xyz]{2}", expr)}


def _sym2(R):
    basis = []
    for i, j in [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]:
        T = np.zeros((3, 3))
        if i == j:
            T[i, i] = 1
        else:
            T[i, j] = T[j, i] = 1 / np.sqrt(2)
        basis.append(T)
    return np.array([[np.sum(a * (R @ b @ R.T)) for b in basis] for a in basis])


def _characters_on_ops(pg, ops):
    """表の各既約表現について、操作順の (複素) 指標と次元 d を返す [(名前, chi, d)]。"""
    if pg in COMPLEX_TABLES:
        CT = COMPLEX_TABLES[pg]
        _, paxis = M._principal_axis(ops)
        col = []
        for R, _ in ops:
            key = M.op_signed_key(R, paxis)
            for idx, (kind, ang) in enumerate(CT["fine"]):
                if kind == key[0] and (abs(abs(ang) - abs(key[1])) < 3 if CT.get("real_pseudo")
                                       else abs(((ang - key[1] + 180) % 360) - 180) < 3):
                    col.append(idx)
                    break
            else:
                raise AssertionError(f"{pg}: unassigned operation {key}")
        out = []
        for name, chars in CT["ctab"].items():
            chi = np.array([chars[c] for c in col], complex)
            if CT["dims"][name] == 2 and not CT.get("real_pseudo"):
                chi = 2 * chi.real                          # 共役ペアの和 (実 2 次元)
            out.append((name, chi, CT["dims"][name]))
        return out
    T = CHAR_TABLES[pg]
    # 一般位置の軌道に加え、x 軸上・xz 面上の原子を置く: C2' = x 軸、σv = xz 面 (原子を通る方を
    # C2', σv と呼ぶ Mulliken の慣例) となり、表の x²−y² (B1) などの成分の向きと一致する
    cart = np.array([R @ p for R, _ in ops
                     for p in (np.array([0.31, 0.17, 0.53]), np.array([0.6, 0, 0]),
                               np.array([0.7, 0, 0.2]))])
    if pg in ("C2v", "D2", "D2h"):
        frame = np.eye(3)                                   # 理想の向き: z = 主軸, x = C2' / σ 法線
        col = [T["classes"].index(M._abelian_class(R, frame, pg)) for R, _ in ops]
    elif pg in ("I", "Ih"):
        from symvib.molecule_tables import ICOSAHEDRAL_CLASS_KEYS
        keys = [M.op_angle_key(R) for R, _ in ops]
        inv = {v: k for k, v in ICOSAHEDRAL_CLASS_KEYS.items()}
        col = [T["classes"].index(inv[k]) for k in keys]
    else:
        col = M._assign_classes(pg, ops, cart)
    table = np.array(T["table"], float)
    return [(name, table[i][col].astype(complex), int(round(table[i][0])))
            for i, name in enumerate(T["irreps"])]


@pytest.mark.parametrize("pg", sorted(set(CHAR_TABLES) | set(COMPLEX_TABLES)))
def test_table_rows_are_projectors_and_activities_are_right(pg):
    group = _closure(generators(pg))
    ops = [(R, np.arange(1)) for R in group]
    h = len(group)
    rows = _characters_on_ops(pg, ops)
    order = M._group_order(pg)
    assert h == order, (pg, h, order)
    mult, inv, _ = M._group_structure(ops)
    g_of = mult[:, inv]
    dim_sum = 0
    for name, chi, d in rows:
        norm = np.vdot(chi, chi).real / h
        physical_pair = abs(norm - 2) < 1e-6                # 複素ペアをまとめた実 2 次元表現
        assert abs(norm - 1) < 1e-6 or physical_pair, (pg, name, norm)
        # 複素ペアは 1 次元既約表現 2 つの和: 正則表現での射影の階数は 1 + 1 = 2
        P = (d / h if not physical_pair else 1.0 / h) * chi[g_of]
        assert np.abs(P @ P - P).max() < 1e-6, (pg, name)
        expected_rank = 2 if physical_pair else d * d
        assert abs(np.trace(P).real - expected_rank) < 1e-6, (pg, name, np.trace(P))
        dim_sum += expected_rank
    assert dim_sum == h, pg                                 # Σ d² = h (全既約表現がそろっている)
    # 活性: ベクトル表現・対称 2 階テンソル表現の分解
    table = CHAR_TABLES.get(pg) or COMPLEX_TABLES[pg]
    chi_vec = np.array([np.trace(R) for R in group])
    T6 = [_sym2(R) for R in group]
    chi_sym2 = np.array([np.trace(t) for t in T6])
    linear = table.get("linear", {})
    quadratic = table.get("quadratic", {})
    for name, chi, d in rows:
        n2 = np.vdot(chi, chi).real
        n_vec = np.vdot(chi, chi_vec).real / n2
        n_q = np.vdot(chi, chi_sym2).real / n2
        assert (n_vec > 0.5) == (name in linear), (pg, name, "IR")
        assert (n_q > 0.5) == (name in quadratic), (pg, name, "Raman")
        if name in quadratic:
            P = sum(np.conj(c) * t for c, t in zip(chi, T6)) / n2
            comps = {["xx", "yy", "zz", "xy", "xz", "yz"][a] for a in range(6) if abs(P[a, a]) > 1e-6}
            listed = set().union(*(_components(e) for e in quadratic[name]))
            assert comps == listed, (pg, name, comps, listed)
        if name in linear:
            V = sum(np.conj(c) * R for c, R in zip(chi, group)) / n2
            comps = {"xyz"[a] for a in range(3) if abs(V[a, a]) > 1e-6}
            assert comps == set(linear[name]), (pg, name, comps, linear[name])
