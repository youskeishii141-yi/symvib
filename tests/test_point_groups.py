"""点群の分解を、コードのクラス割り当てに頼らない独立な検査で確かめる。

各モードのブロックについて
  - 指標のノルム |χ|²/h が 1 (既約) であること
  - 同じラベルのブロックは同じ指標、違うラベルは直交する指標であること
  - IR / Raman 活性が、ベクトル表現・対称 2 階テンソル表現との重なりと一致すること
を、見つけた対称操作と変位から直接計算して確かめる。大きめの分子・任意の向き・
原子の順番の入れ替えでも結果が正しいことを確認する。
"""
import math

import numpy as np
import pytest

from structures import _closure, _mol, _orbit, _ring, _rot
from symvib.common import atomic_mass
from symvib import molecule as M

SH = np.diag([1., 1., -1.])


def _td():
    return _closure([_rot((1, 1, 1), 2 * math.pi / 3), _rot((0, 0, 1), math.pi / 2) @ SH])


def _oh():
    return _closure([_rot((1, 1, 1), 2 * math.pi / 3), _rot((0, 0, 1), math.pi / 2), -np.eye(3)])


def _molecules():
    td, oh = _td(), _oh()
    a, b = 0.89, 1.78
    tw = math.radians(15)
    d2d = _closure([_rot((0, 0, 1), math.pi / 2) @ SH, np.diag([1., -1., 1.])])
    return {
        # 立方群 (以前は網目モデルのフォールバックで誤った分解になっていた)
        "adamantane": ("Td", _mol(_orbit("C", (a, a, a), td) + _orbit("C", (b, 0, 0), td)
                                  + _orbit("H", (a + .63, a + .63, a + .63), td)
                                  + _orbit("H", (b + .63, .77, -.77), td)),
                       "5A1 + A2 + 6E + 7T1 + 11T2"),
        "Ni(CO)4": ("Td", _mol([("Ni", (0, 0, 0))] + _orbit("C", (1.05, 1.05, 1.05), td)
                               + _orbit("O", (1.71, 1.71, 1.71), td)),
                    "2A1 + 2E + T1 + 4T2"),
        "Mo(CO)6": ("Oh", _mol([("Mo", (0, 0, 0))] + _orbit("C", (2.06, 0, 0), oh)
                               + _orbit("O", (3.21, 0, 0), oh)),
                    "2A1g + 2Eg + T1g + 2T2g + 4T1u + 2T2u"),
        "cubane": ("Oh", _mol(_orbit("C", (.78, .78, .78), oh)
                              + _orbit("H", (1.41, 1.41, 1.41), oh)), None),
        # 主軸が 5, 7, 8, 12 回の群・D5 (以前はクラス名の対応や指標表の誤りで失敗)
        "ferrocene D5": ("D5", _mol([("Fe", (0, 0, 0))]
                                    + _ring("C", 5, 1.22, 1.66, tw / 2)
                                    + _ring("C", 5, 1.22, -1.66, -tw / 2)
                                    + _ring("H", 5, 2.30, 1.66, tw / 2)
                                    + _ring("H", 5, 2.30, -1.66, -tw / 2)), None),
        "tropylium": ("D7h", _mol(_ring("C", 7, 1.61) + _ring("H", 7, 2.70)),
                      "2A1' + A2' + 3E1' + 4E2' + 4E3' + A2'' + E1'' + 2E2'' + 2E3''"),
        "C8H8 2-": ("D8h", _mol(_ring("C", 8, 1.84) + _ring("H", 8, 2.93)), None),
        "Cr(C6H6)2 D6d": ("D6d", _mol([("Cr", (0, 0, 0))] + _ring("C", 6, 1.42, 1.61)
                                      + _ring("C", 6, 1.42, -1.61, math.pi / 6)
                                      + _ring("H", 6, 2.50, 1.61)
                                      + _ring("H", 6, 2.50, -1.61, math.pi / 6)), None),
        "IF7": ("D5h", _mol([("I", (0, 0, 0)), ("F", (0, 0, 1.79)), ("F", (0, 0, -1.79))]
                            + _ring("F", 5, 1.86)), "2A1' + 3E1' + 2E2' + 2A2'' + E1'' + E2''"),
        "spiropentane": ("D2d", _mol([("C", (0, 0, 0))] + _orbit("C", (0.755, 0, 1.30), d2d)
                                     + _orbit("H", (1.30, 0.92, 1.58), d2d)),
                         "5A1 + 3A2 + 4B1 + 5B2 + 8E"),
        # Cs: 鏡面が z に垂直でない向きでも正しく
        "formic acid": ("Cs", _mol([("C", (0, 0.42, 0)), ("O", (1.09, 0.92, 0)),
                                    ("O", (-1.15, 1.13, 0)), ("H", (-0.02, -0.68, 0)),
                                    ("H", (-1.88, 0.49, 0))]), "7A' + 2A''"),
        # ほぼ直線だが直線ではない分子 (以前は C∞v と誤判定し、原子を軸上へ動かしていた)
        "HNCS": ("Cs", _mol([("H", (-0.73, 0.68, 0)), ("N", (0, 0, 0)), ("C", (1.21, 0, 0)),
                             ("S", (2.77, 0.03, 0))]), "5A' + A''"),
        "CH3(CC)3H": ("C3v", _mol([("C", (0, 0, 0))]
                                  + [("C", (0, 0, z)) for k in range(3)
                                     for z in (1.46 + 2.58 * k, 2.67 + 2.58 * k)]
                                  + [("H", (0, 0, 8.86))] + _ring("H", 3, 1.03, -0.36)),
                      "9A1 + 9E"),
    }


def _validate(r, numbers):
    """独立な検査。問題のリストを返す (空なら正しい)。"""
    probs = []
    na = len(numbers)
    if sum(m.degeneracy for m in r.modes) != 3 * na - (5 if r.is_linear else 6):
        probs.append("degrees of freedom")
    ops = M.find_point_group_operations(r.cart, numbers, 1e-2)
    h = len(ops)
    sm = np.sqrt([atomic_mass(z) for z in numbers])
    dmats = M.mechanical_matrices(ops, na)
    chiv = np.array([np.trace(R) for R, _ in ops])
    chiq = np.array([(np.trace(R) ** 2 + np.trace(R @ R)) / 2 for R, _ in ops])
    by_label = {}
    for m in r.modes:
        B = np.array([(u * sm[:, None]).reshape(-1) for u in m.displacements])
        q, _ = np.linalg.qr(B.T)
        chi = np.array([np.trace(q.T @ D @ q) for D in dmats])
        norm = chi @ chi / h
        if not (abs(norm - 1) < 1e-4 or abs(norm - 2) < 1e-4):
            probs.append(f"{m.label} reducible ({norm:.3f})")
        if m.label in by_label and not np.allclose(by_label[m.label], chi, atol=1e-4):
            probs.append(f"{m.label} inconsistent")
        by_label.setdefault(m.label, chi)
        ir, ra = abs(chi @ chiv) > 1e-4, abs(chi @ chiq) > 1e-4
        act = "Raman+IR" if ir and ra else "Raman" if ra else "IR" if ir else "Silent"
        if act != m.activity:
            probs.append(f"{m.label} activity {m.activity} (true {act})")
    labels = list(by_label)
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            if abs(by_label[labels[i]] @ by_label[labels[j]]) > 1e-4:
                probs.append(f"{labels[i]} / {labels[j]} are the same irrep")
    return probs


@pytest.mark.parametrize("name", list(_molecules()))
def test_decomposition_is_correct_in_any_orientation_and_order(name):
    pg, (cart, Z, sym), expected = _molecules()[name]
    rng = np.random.default_rng(sum(map(ord, name)))
    for trial in range(3):
        c, z, s = np.asarray(cart, float), list(Z), list(sym)
        if trial:                                   # 任意の回転と原子の並べ替え
            q, _ = np.linalg.qr(rng.standard_normal((3, 3)))
            perm = rng.permutation(len(z))
            c, z, s = (c @ q.T)[perm], [z[i] for i in perm], [s[i] for i in perm]
        r = M.analyze_molecule(c, z, s)
        assert r.point_group == pg
        assert _validate(r, r.numbers) == []
        if expected:
            assert r.decomposition.split("= ")[1] == expected


def test_every_character_table_is_consistent_with_its_group():
    """全指標表について、理想的な群を作ってクラスの割り当てが矛盾なくできること
    (列の入れ替えのような、直交性の検査では見つからない表の誤りも検出する)。"""
    from symvib.molecule_tables import CHAR_TABLES
    import re
    phi = (1 + 5 ** 0.5) / 2

    def gens(name):
        z = (0, 0, 1)
        special = {
            "Td": [_rot((1, 1, 1), 2 * math.pi / 3), _rot(z, math.pi / 2) @ SH],
            "Oh": [_rot((1, 1, 1), 2 * math.pi / 3), _rot(z, math.pi / 2), -np.eye(3)],
            "O": [_rot((1, 1, 1), 2 * math.pi / 3), _rot(z, math.pi / 2)],
            "I": [_rot((0, 1, phi), 2 * math.pi / 5), _rot(z, math.pi)],
            "Ih": [_rot((0, 1, phi), 2 * math.pi / 5), _rot(z, math.pi), -np.eye(3)],
            "Cs": [SH], "Ci": [-np.eye(3)], "C1": [np.eye(3)],
        }
        if name in special:
            return special[name]
        kind, n, suf = re.fullmatch(r"([CDS])(\d+)([vhd]?)", name).groups()
        n = int(n)
        cn = _rot(z, 2 * math.pi / n)
        g = [cn] if kind == "C" else [cn, np.diag([1., -1., -1.])]
        if suf == "v":
            g.append(np.diag([1., -1., 1.]))
        if suf == "h":
            g.append(SH)
        if suf == "d" and kind == "D":
            g.append(_rot(z, math.pi / n) @ SH)
        return g

    for pg, T in CHAR_TABLES.items():
        if pg in ("C2v", "D2", "D2h"):
            continue                                # 軸の慣例で割り当てる (別の検査で確認)
        group = _closure(gens(pg))
        ops = [(R, np.arange(1)) for R in group]
        assert len(ops) == int(round(sum(T["g"]))), pg
        # 一般の位置の原子 (正則表現) で割り当てが矛盾なくできる
        point = np.array([0.31, 0.17, 0.53])
        cart = np.array([R @ point for R in group])
        assert M._point_group_from_ops(ops) == pg
        M._assign_classes(pg, ops, cart)            # 矛盾があれば例外
