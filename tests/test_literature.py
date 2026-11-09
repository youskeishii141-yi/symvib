"""文献値との照合: よく研究された結晶の Γ 点既約分解・Raman / IR 活性モードの数。

構造は Wyckoff 位置の文献値 (座標はおおよそでよい; 対称性で決まる) から spglib の
空間群データベースの操作で展開する。期待値の出典は各項目のコメントに示す。
ラベルは空間群の標準設定 (Bilbao と同じ) で表すので、非標準設定 (Pbnm など) の文献値は
軸の読み替えを明記して変換している。
"""
import re
from collections import Counter

import numpy as np
import pytest
import spglib
from pymatgen.core import Lattice

from symvib.crystal import analyze_crystal


def _hall(symbol, choice=""):
    """H-M 記号 (ITA full / short) と設定記号から spglib の Hall 番号を引く。"""
    norm = lambda s: re.sub(r"[\s_]", "", s).lower()                 # noqa: E731
    for h in range(1, 531):
        t = spglib.get_spacegroup_type(h)
        if norm(symbol) in (norm(t.international_full), norm(t.international_short),
                            norm(t.international)) and (choice == "" or t.choice == choice):
            return h
    raise KeyError(symbol)


def structure(symbol, cell, atoms, choice=""):
    """(格子, 分率座標, 原子番号) を Wyckoff 位置の代表座標から作る。atoms: [(元素, (x,y,z))]。"""
    from pymatgen.core import Element
    hall = _hall(symbol, choice)
    sym = spglib.get_symmetry_from_database(hall)
    lattice = Lattice.from_parameters(*cell).matrix
    frac, numbers = [], []
    for el, p in atoms:
        orb = []
        for R, t in zip(sym["rotations"], sym["translations"]):
            f = (R @ np.asarray(p, float) + t) % 1.0
            f[f > 1 - 1e-8] = 0.0
            if not any(np.linalg.norm((((f - q) + 0.5) % 1 - 0.5) @ lattice) < 1e-3 for q in orb):
                orb.append(f)
        frac += orb
        numbers += [Element(el).Z] * len(orb)
    return lattice, np.array(frac), numbers


def _counts(result):
    opt = result.optical_modes
    return (dict(Counter(m.label for m in opt if "Raman" in m.activity)),
            dict(Counter(m.label for m in opt if "IR" in m.activity)),
            dict(Counter(m.label for m in opt if m.activity == "Silent")),
            dict(Counter(m.label for m in result.modes if m.acoustic)))


C = {
    # ── 立方晶 ──
    # ダイヤモンド: Raman T2g (F2g) 1 本 (Solin & Ramdas, PRB 1, 1687 (1970))
    "diamond": (("F d -3 m", (3.567,) * 3 + (90,) * 3, [("C", (0, 0, 0))], "1"),
                {"T2g": 1}, {}, {}),
    # NaCl: 光学は T1u (IR) のみ
    "NaCl": (("F m -3 m", (5.64,) * 3 + (90,) * 3, [("Na", (0, 0, 0)), ("Cl", (.5, .5, .5))]),
             {}, {"T1u": 1}, {}),
    # CsCl: T1u (IR)
    "CsCl": (("P m -3 m", (4.12,) * 3 + (90,) * 3, [("Cs", (0, 0, 0)), ("Cl", (.5, .5, .5))]),
             {}, {"T1u": 1}, {}),
    # 蛍石 CaF2: Raman T2g, IR T1u (Russell, Proc. Phys. Soc. 85, 194 (1965))
    "CaF2": (("F m -3 m", (5.46,) * 3 + (90,) * 3, [("Ca", (0, 0, 0)), ("F", (.25, .25, .25))]),
             {"T2g": 1}, {"T1u": 1}, {}),
    # 閃亜鉛鉱 ZnS: T2 (Raman+IR)
    "ZnS": (("F -4 3 m", (5.41,) * 3 + (90,) * 3, [("Zn", (0, 0, 0)), ("S", (.25, .25, .25))]),
            {"T2": 1}, {"T2": 1}, {}),
    # SrTiO3 (立方ペロブスカイト): 3T1u (IR) + T2u (silent) (Cowley, Phys. Rev. 134, A981 (1964))
    "SrTiO3": (("P m -3 m", (3.905,) * 3 + (90,) * 3,
                [("Sr", (0, 0, 0)), ("Ti", (.5, .5, .5)), ("O", (.5, .5, 0))]),
               {}, {"T1u": 3}, {"T2u": 1}),
    # スピネル MgAl2O4: Raman A1g + Eg + 3T2g, IR 4T1u (Fraas, Moore, Salzberg, JCP 58, 3585 (1973))
    "spinel": (("F d -3 m", (8.08,) * 3 + (90,) * 3,
                [("Mg", (.125, .125, .125)), ("Al", (.5, .5, .5)), ("O", (.2624,) * 3)], "2"),
               {"A1g": 1, "Eg": 1, "T2g": 3}, {"T1u": 4}, {"A2u": 2, "Eu": 2, "T1g": 1, "T2u": 2}),
    # YAG Y3Al5O12: Raman 3A1g + 8Eg + 14T2g, IR 17T1u (Hurrell et al., Phys. Rev. 173, 851 (1968))
    "YAG": (("I a -3 d", (12.0,) * 3 + (90,) * 3,
             [("Y", (.125, 0, .25)), ("Al", (0, 0, 0)), ("Al", (.375, 0, .25)),
              ("O", (-.0306, .0512, .15))]),
            {"A1g": 3, "Eg": 8, "T2g": 14}, {"T1u": 17},
            {"A2g": 5, "T1g": 14, "A1u": 5, "A2u": 5, "Eu": 10, "T2u": 16}),
    # 黄鉄鉱 FeS2 (Pa-3): Raman Ag + Eg + 3Tg, IR 5Tu (Vogt, Chattopadhyay, Stolz, JPCS 44, 869 (1983))
    "pyrite": (("P a -3", (5.418,) * 3 + (90,) * 3, [("Fe", (0, 0, 0)), ("S", (.385,) * 3)]),
               {"Ag": 1, "Eg": 1, "Tg": 3}, {"Tu": 5}, {"Au": 2, "Eu": 2}),
    # Cu2O (Pn-3m): Raman T2g, IR 2T1u, silent A2u + Eu + T2u (Carabatos & Prevot, pss(b) 44, 701 (1971))
    "Cu2O": (("P n -3 m", (4.27,) * 3 + (90,) * 3, [("O", (0, 0, 0)), ("Cu", (.25, .25, .25))], "1"),
             {"T2g": 1}, {"T1u": 2}, {"A2u": 1, "Eu": 1, "T2u": 1}),
    # ── 正方晶 ──
    # ルチル TiO2: Raman A1g + B1g + B2g + Eg, IR A2u + 3Eu, silent A2g + 2B1u
    # (Porto, Fleury, Damen, Phys. Rev. 154, 522 (1967))
    "rutile": (("P 42/m n m", (4.594, 4.594, 2.959, 90, 90, 90),
                [("Ti", (0, 0, 0)), ("O", (.305, .305, 0))]),
               {"A1g": 1, "B1g": 1, "B2g": 1, "Eg": 1}, {"A2u": 1, "Eu": 3}, {"A2g": 1, "B1u": 2}),
    # アナターゼ TiO2: Raman A1g + 2B1g + 3Eg, IR A2u + 2Eu, silent B2u
    # (Ohsaka, Izumi, Fujiki, J. Raman Spectrosc. 7, 321 (1978))
    "anatase": (("I 41/a m d", (3.785, 3.785, 9.514, 90, 90, 90),
                 [("Ti", (0, .75, .125)), ("O", (0, .25, .0816))], "2"),
                {"A1g": 1, "B1g": 2, "Eg": 3}, {"A2u": 1, "Eu": 2}, {"B2u": 1}),
    # ジルコン ZrSiO4: Raman 2A1g + 4B1g + B2g + 5Eg, IR 3A2u + 4Eu
    # (Dawson, Hargreave, Wilkinson, J. Phys. C 4, 240 (1971))
    "zircon": (("I 41/a m d", (6.607, 6.607, 5.982, 90, 90, 90),
                [("Zr", (0, .75, .125)), ("Si", (0, .25, .375)), ("O", (0, .066, .1951))], "2"),
               {"A1g": 2, "B1g": 4, "B2g": 1, "Eg": 5}, {"A2u": 3, "Eu": 4},
               {"A2g": 1, "A1u": 1, "B1u": 1, "B2u": 2}),
    # シーライト CaWO4 (I4_1/a): Raman 3Ag + 5Bg + 5Eg, IR 4Au + 4Eu, silent 3Bu
    # (Porto & Scott, Phys. Rev. 157, 716 (1967))
    "scheelite": (("I 41/a", (5.243, 5.243, 11.376, 90, 90, 90),
                   [("W", (0, .25, .125)), ("Ca", (0, .25, .625)), ("O", (.1504, .0085, .2111))], "2"),
                  {"Ag": 3, "Bg": 5, "Eg": 5}, {"Au": 4, "Eu": 4}, {"Bu": 3}),
    # 黄銅鉱 CuFeS2 (I-42d): Raman A1 + 3B1 + 3B2 + 6E, IR 3B2 + 6E, silent 2A2
    # (Koschel, Sorger, Baars, J. Phys. Colloq. 36, C3-177 (1975))
    "chalcopyrite": (("I -4 2 d", (5.29, 5.29, 10.42, 90, 90, 90),
                      [("Cu", (0, 0, 0)), ("Fe", (0, 0, .5)), ("S", (.257, .25, .125))]),
                     {"A1": 1, "B1": 3, "B2": 3, "E": 6}, {"B2": 3, "E": 6}, {"A2": 2}),
    # PbTiO3 (P4mm): 3A1 + B1 + 4E (A1, E: Raman+IR; B1: Raman) (Burns & Scott, PRB 7, 3088 (1973))
    "PbTiO3": (("P 4 m m", (3.90, 3.90, 4.15, 90, 90, 90),
                [("Pb", (0, 0, 0)), ("Ti", (.5, .5, .54)), ("O", (.5, .5, .11)), ("O", (.5, 0, .62))]),
               {"A1": 3, "B1": 1, "E": 4}, {"A1": 3, "E": 4}, {}),
    # ── 三方晶 ──
    # α-石英 (P3_121): Raman 4A1 + 8E, IR 4A2 + 8E (Scott & Porto, Phys. Rev. 161, 903 (1967))
    "quartz": (("P 31 2 1", (4.913, 4.913, 5.405, 90, 90, 120),
                [("Si", (.4697, 0, 1 / 3)), ("O", (.4135, .2669, .1191))]),
               {"A1": 4, "E": 8}, {"A2": 4, "E": 8}, {}),
    # コランダム Al2O3 (R-3c): Raman 2A1g + 5Eg, IR 2A2u + 4Eu (Porto & Krishnan, JCP 47, 1009 (1967))
    "corundum": (("R -3 c", (4.759, 4.759, 12.991, 90, 90, 120),
                  [("Al", (0, 0, .3522)), ("O", (.3062, 0, .25))], "H"),
                 {"A1g": 2, "Eg": 5}, {"A2u": 2, "Eu": 4}, {"A1u": 2, "A2g": 3}),
    # 方解石 CaCO3 (R-3c): Raman A1g + 4Eg, IR 3A2u + 5Eu (Porto, Giordmaine, Damen, Phys. Rev. 147, 608 (1966))
    "calcite": (("R -3 c", (4.99, 4.99, 17.06, 90, 90, 120),
                 [("Ca", (0, 0, 0)), ("C", (0, 0, .25)), ("O", (.257, 0, .25))], "H"),
                {"A1g": 1, "Eg": 4}, {"A2u": 3, "Eu": 5}, {"A1u": 2, "A2g": 3}),
    # LiNbO3 (R3c): 4A1 + 9E (Raman+IR) + 5A2 (silent) (Barker & Loudon, Phys. Rev. 158, 433 (1967))
    "LiNbO3": (("R 3 c", (5.148, 5.148, 13.863, 90, 90, 120),
                [("Li", (0, 0, .2829)), ("Nb", (0, 0, 0)), ("O", (.0492, .3446, .0647))], "H"),
               {"A1": 4, "E": 9}, {"A1": 4, "E": 9}, {"A2": 5}),
    # Bi2Se3 (R-3m): Raman 2A1g + 2Eg, IR 2A2u + 2Eu (Richter, Köhler, Becker, pss(b) 84, 619 (1977))
    "Bi2Se3": (("R -3 m", (4.14, 4.14, 28.64, 90, 90, 120),
                [("Se", (0, 0, 0)), ("Se", (0, 0, .2109)), ("Bi", (0, 0, .4008))], "H"),
               {"A1g": 2, "Eg": 2}, {"A2u": 2, "Eu": 2}, {}),
    # イルメナイト FeTiO3 (R-3): Raman 5Ag + 5Eg, IR 4Au + 4Eu (White, in Infrared & Raman Spectroscopy of Lunar and Terrestrial Minerals (1975))
    "ilmenite": (("R -3", (5.088, 5.088, 14.085, 90, 90, 120),
                  [("Fe", (0, 0, .3553)), ("Ti", (0, 0, .1464)), ("O", (.3172, .0234, .2448))], "H"),
                 {"Ag": 5, "Eg": 5}, {"Au": 4, "Eu": 4}, {}),
    # ── 六方晶 ──
    # 黒鉛 2H: Raman 2E2g, IR A2u + E1u, silent 2B2g (Nemanich & Solin, PRB 20, 392 (1979))
    "graphite": (("P 63/m m c", (2.461, 2.461, 6.708, 90, 90, 120),
                  [("C", (0, 0, .25)), ("C", (1 / 3, 2 / 3, .25))]),
                 {"E2g": 2}, {"A2u": 1, "E1u": 1}, {"B2g": 2}),
    # h-BN: Raman 2E2g, IR A2u + E1u, silent 2 本 (Geick, Perry, Rupprecht, Phys. Rev. 146, 543 (1966))。
    # 文献では silent を B1g と書くが、これは C2' / C2'' の呼び方 (どちらを添字 1 の基準にするか) の
    # 違いで、Bilbao の標準設定 (添字 1 = [120] 方向の 2 回軸に対称) では黒鉛と同じく B2g になる
    # (指標: 2_[100] → +1, 2_[120] (z=1/4) → −1)。
    "hBN": (("P 63/m m c", (2.504, 2.504, 6.661, 90, 90, 120),
             [("B", (1 / 3, 2 / 3, .25)), ("N", (1 / 3, 2 / 3, .75))]),
            {"E2g": 2}, {"A2u": 1, "E1u": 1}, {"B2g": 2}),
    # 2H-MoS2: Raman A1g + E1g + 2E2g, IR A2u + E1u, silent 2B2g + B1u + E2u
    # (Verble & Wieting, PRL 25, 362 (1970))
    "MoS2": (("P 63/m m c", (3.16, 3.16, 12.29, 90, 90, 120),
              [("Mo", (1 / 3, 2 / 3, .25)), ("S", (1 / 3, 2 / 3, .621))]),
             {"A1g": 1, "E1g": 1, "E2g": 2}, {"A2u": 1, "E1u": 1}, {"B2g": 2, "B1u": 1, "E2u": 1}),
    # ウルツ鉱 ZnO (P6_3mc): A1 + E1 (Raman+IR), 2E2 (Raman), 2B1 (silent)
    # (Damen, Porto, Tell, Phys. Rev. 142, 570 (1966))
    "ZnO": (("P 63 m c", (3.25, 3.25, 5.207, 90, 90, 120),
             [("Zn", (1 / 3, 2 / 3, 0)), ("O", (1 / 3, 2 / 3, .382))]),
            {"A1": 1, "E1": 1, "E2": 2}, {"A1": 1, "E1": 1}, {"B1": 2}),
    # ── 直方晶 ──
    # フォルステライト Mg2SiO4 (Pbnm): Raman 11Ag + 11B1g + 7B2g + 7B3g, IR 9B1u + 13B2u + 13B3u,
    # silent 10Au (Iishi, Am. Mineral. 63, 1198 (1978))。
    # 標準設定 Pnma への読み替え (Pbnm の a,b,c = Pnma の c,a,b): 添字 1→2, 2→3, 3→1
    "forsterite": (("P b n m", (4.754, 10.197, 5.980, 90, 90, 90),
                    [("Mg", (0, 0, 0)), ("Mg", (.9915, .2774, .25)), ("Si", (.4262, .0940, .25)),
                     ("O", (.7657, .0913, .25)), ("O", (.2215, .4474, .25)), ("O", (.2777, .1628, .0331))],
                    "cab"),
                   {"Ag": 11, "B1g": 7, "B2g": 11, "B3g": 7}, {"B1u": 13, "B2u": 9, "B3u": 13},
                   {"Au": 10}),
    # KTiOPO4 (Pna2_1): 47A1 + 48A2 + 47B1 + 47B2 (Kugel et al., J. Phys. C 21, 5565 (1988))
    "KTP": (("P n a 21", (12.814, 6.404, 10.616, 90, 90, 90),
             [("K", (.378, .780, .312)), ("K", (.105, .699, .066)), ("Ti", (.373, .500, 0)),
              ("Ti", (.247, .269, .251)), ("P", (.498, .337, .260)), ("P", (.181, .502, .512)),
              ("O", (.486, .487, .150)), ("O", (.510, .466, .383)), ("O", (.400, .199, .281)),
              ("O", (.594, .193, .241)), ("O", (.112, .311, .541)), ("O", (.112, .692, .487)),
              ("O", (.253, .540, .627)), ("O", (.253, .461, .399)), ("O", (.224, .041, .609)),
              ("O", (.224, .967, .404))]),
            {"A1": 47, "A2": 48, "B1": 47, "B2": 47}, {"A1": 47, "B1": 47, "B2": 47}, {}),
    # ── 単斜晶 ──
    # バデレアイト ZrO2 (P2_1/c): Raman 9Ag + 9Bg, IR 8Au + 7Bu (Feinberg & Perry, JPCS 42, 513 (1981))
    "baddeleyite": (("P 1 21/c 1", (5.15, 5.21, 5.32, 90, 99.2, 90),
                     [("Zr", (.2754, .0395, .2083)), ("O", (.0700, .3317, .3447)),
                      ("O", (.4496, .7569, .4792))]),
                    {"Ag": 9, "Bg": 9}, {"Au": 8, "Bu": 7}, {}),
    # β-Ga2O3 (C2/m): Raman 10Ag + 5Bg, IR 4Au + 8Bu (Dohy, Lucazeau, Revcolevschi, JSSC 45, 180 (1982))
    "Ga2O3": (("C 1 2/m 1", (12.23, 3.04, 5.80, 90, 103.7, 90),
               [("Ga", (.0905, 0, .7948)), ("Ga", (.3414, 0, .6857)), ("O", (.1674, 0, .1011)),
                ("O", (.4957, 0, .2553)), ("O", (.8279, 0, .4365))]),
              {"Ag": 10, "Bg": 5}, {"Au": 4, "Bu": 8}, {}),
    # 透輝石 CaMgSi2O6 (C2/c): Raman 14Ag + 16Bg, IR 13Au + 14Bu (Rutstein & White, Am. Mineral. 56, 877 (1971))
    "diopside": (("C 1 2/c 1", (9.746, 8.899, 5.251, 90, 105.63, 90),
                  [("Ca", (0, .3015, .25)), ("Mg", (0, .9082, .25)), ("Si", (.2862, .0933, .2293)),
                   ("O", (.1156, .0873, .1422)), ("O", (.3611, .2500, .3180)), ("O", (.3505, .0176, .9953))]),
                 {"Ag": 14, "Bg": 16}, {"Au": 13, "Bu": 14}, {}),
    # ── 三斜晶 ──
    # カイヤナイト Al2SiO5 (P-1, 32 原子): Raman 48Ag, IR 45Au (Mernagh & Liu, Phys. Chem. Minerals 18, 126 (1991))
    "kyanite": (("P -1", (7.126, 7.852, 5.572, 89.99, 101.11, 106.03),
                 [("Al", (.3254, .7040, .4582)), ("Al", (.2974, .6989, .9505)), ("Al", (.0998, .3862, .6403)),
                  ("Al", (.1120, .9175, .1649)), ("Si", (.2962, .0649, .7066)), ("Si", (.2910, .3317, .1892)),
                  ("O", (.1095, .1468, .1288)), ("O", (.1230, .6856, .1812)), ("O", (.2747, .4545, .9547)),
                  ("O", (.2831, .9354, .9353)), ("O", (.1081, .1520, .6669)), ("O", (.1219, .6307, .6389)),
                  ("O", (.2822, .4453, .4288)), ("O", (.2915, .9467, .4659)), ("O", (.5008, .2749, .2440)),
                  ("O", (.5015, .2312, .7553))]),
                {"Ag": 48}, {"Au": 45}, {}),
}


@pytest.mark.parametrize("name", sorted(C))
def test_literature_mode_counts(name):
    (symbol, cell, atoms, *choice), raman, ir, silent = C[name]
    lat, frac, num = structure(symbol, cell, atoms, choice[0] if choice else "")
    r = analyze_crystal(lat, frac, num)
    got_raman, got_ir, got_silent, acoustic = _counts(r)
    assert got_raman == raman, (name, r.decomposition)
    assert got_ir == ir, (name, r.decomposition)
    assert got_silent == silent, (name, r.decomposition)
    assert sum(m.degeneracy for m in r.modes if m.acoustic) == 3


def test_b2g_of_graphite_and_hbn_follow_bilbao_convention():
    """黒鉛 (2b, 2c) と h-BN (2c, 2d) の面外 silent モードはどちらも Bilbao の慣例で B2g:
    [100] 方向の 2 回軸で +1、z=1/4 を通る [120] 方向の 2 回軸で −1 (直接計算で確かめる)。"""
    import spglib
    from invariants import cart_rotations
    for name in ("graphite", "hBN"):
        (sym, cell, atoms), *_ = C[name]
        lat, frac, num = structure(sym, cell, atoms)
        r = analyze_crystal(lat, frac, num)
        silent = [m for m in r.modes if m.activity == "Silent"]
        assert {m.label for m in silent} == {"B2g"}
        ops = spglib.get_symmetry((lat, frac, num), symprec=1e-3)
        R = cart_rotations(np.asarray(lat), ops["rotations"])
        for m in silent:
            u = m.displacements[0]
            for W, t, Rc in zip(ops["rotations"], ops["translations"], R):
                if not (np.linalg.det(W) > 0 and np.trace(W) == -1):   # 面内の 2 回軸だけ
                    continue
                axis = np.linalg.eig(W)[1][:, np.argmin(np.abs(np.linalg.eigvals(W) - 1))].real
                axis = axis / axis[np.argmax(np.abs(axis))]
                moved = np.zeros_like(u)
                for i in range(len(num)):
                    img = (W @ frac[i] + t) % 1.0
                    j = next(k for k in range(len(num)) if num[k] == num[i]
                             and np.linalg.norm((((img - frac[k]) + .5) % 1 - .5) @ lat) < 1e-3)
                    moved[j] = Rc @ u[i]
                chi = float(np.sum(moved * u) / np.sum(u * u))
                if np.allclose(np.abs(axis), [1, 0, 0]) or np.allclose(np.abs(axis), [0, 1, 0]):
                    assert np.isclose(chi, +1), (name, axis)        # 2_[100], 2_[010]
                elif np.allclose(np.abs(axis), [1, 2, 0]) or np.allclose(np.abs(axis), [2, 1, 0]):
                    assert np.isclose(chi, -1), (name, axis)        # 2_[120], 2_[210]
