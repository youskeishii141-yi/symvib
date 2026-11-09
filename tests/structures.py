"""回帰テスト用の分子・結晶構造ジェネレータ。

分子は (cart[Å], Z, symbols)、結晶は pymatgen の Structure を返す。
いずれも理想幾何から作るので、対称性は厳密に成り立つ。
"""
import itertools
import math

import numpy as np

# ── 分子 ────────────────────────────────────────────────────────────────


def _mol(atoms):
    """[(symbol, (x, y, z)), ...] → (cart, Z, symbols)。"""
    from pymatgen.core import Element
    sym = [a[0] for a in atoms]
    cart = np.array([a[1] for a in atoms], float)
    Z = [Element(s).Z for s in sym]
    return cart, Z, sym


def _ring(sym, n, r, z=0.0, phase=0.0):
    return [(sym, (r * math.cos(phase + 2 * math.pi * k / n),
                   r * math.sin(phase + 2 * math.pi * k / n), z)) for k in range(n)]


def _rot(axis, ang):
    ax = np.asarray(axis, float) / np.linalg.norm(axis)
    x, y, z = ax
    c, s = math.cos(ang), math.sin(ang)
    return np.array([[c + x*x*(1-c), x*y*(1-c) - z*s, x*z*(1-c) + y*s],
                     [y*x*(1-c) + z*s, c + y*y*(1-c), y*z*(1-c) - x*s],
                     [z*x*(1-c) - y*s, z*y*(1-c) + x*s, c + z*z*(1-c)]])


def _closure(gens):
    """生成元の行列から有限群の全要素を作る。"""
    group = [np.eye(3)]
    frontier = [np.eye(3)]
    while frontier:
        new = []
        for g in frontier:
            for h in gens:
                p = h @ g
                if not any(np.allclose(p, q, atol=1e-8) for q in group):
                    group.append(p)
                    new.append(p)
        frontier = new
    return group


def _orbit(sym, point, group):
    pts = []
    for g in group:
        p = g @ np.asarray(point, float)
        if not any(np.allclose(p, q, atol=1e-6) for q in pts):
            pts.append(p)
    return [(sym, tuple(p)) for p in pts]


def molecules():
    """名前 → (cart, Z, symbols)。"""
    phi = (1 + math.sqrt(5)) / 2
    m = {}
    m["H2O"] = _mol([("O", (0, 0, 0.117)), ("H", (0, 0.757, -0.467)),
                     ("H", (0, -0.757, -0.467))])
    m["NH3"] = _mol([("N", (0, 0, 0.38))] + _ring("H", 3, 0.94))
    t = 0.63
    m["CH4"] = _mol([("C", (0, 0, 0)), ("H", (t, t, t)), ("H", (-t, -t, t)),
                     ("H", (-t, t, -t)), ("H", (t, -t, -t))])
    m["CO2"] = _mol([("C", (0, 0, 0)), ("O", (0, 0, 1.16)), ("O", (0, 0, -1.16))])
    m["HCN"] = _mol([("H", (0, 0, -1.066)), ("C", (0, 0, 0)), ("N", (0, 0, 1.156))])
    m["C2H2"] = _mol([("C", (0, 0, 0.6)), ("C", (0, 0, -0.6)),
                      ("H", (0, 0, 1.66)), ("H", (0, 0, -1.66))])
    m["C6H6"] = _mol(_ring("C", 6, 1.39) + _ring("H", 6, 2.47))
    s = 1.56
    m["SF6"] = _mol([("S", (0, 0, 0))] + [("F", tuple(v)) for v in
                    (s*np.vstack([np.eye(3), -np.eye(3)]))])
    m["BF3"] = _mol([("B", (0, 0, 0))] + _ring("F", 3, 1.31))
    m["XeF4"] = _mol([("Xe", (0, 0, 0))] + _ring("F", 4, 1.95))
    # エタン (D3d, 互い違い)
    m["C2H6"] = _mol([("C", (0, 0, 0.77)), ("C", (0, 0, -0.77))]
                     + _ring("H", 3, 1.02, 1.16) + _ring("H", 3, 1.02, -1.16, math.pi / 3))
    # アレン (D2d)
    m["C3H4"] = _mol([("C", (0, 0, 0)), ("C", (0, 0, 1.31)), ("C", (0, 0, -1.31)),
                      ("H", (0.93, 0, 1.87)), ("H", (-0.93, 0, 1.87)),
                      ("H", (0, 0.93, -1.87)), ("H", (0, -0.93, -1.87))])
    # トランス-ジクロロエチレン (C2h)
    m["C2H2Cl2"] = _mol([("C", (0.67, 0, 0)), ("C", (-0.67, 0, 0)),
                         ("Cl", (1.5, 1.45, 0)), ("Cl", (-1.5, -1.45, 0)),
                         ("H", (1.2, -0.95, 0)), ("H", (-1.2, 0.95, 0))])
    # フェロセン重なり形 (D5h)
    m["FeCp2_D5h"] = _mol([("Fe", (0, 0, 0))]
                          + _ring("C", 5, 1.22, 1.66) + _ring("C", 5, 1.22, -1.66)
                          + _ring("H", 5, 2.30, 1.66) + _ring("H", 5, 2.30, -1.66))
    # 一般位置軌道から作る複素表現群・高次群
    z = np.array([0, 0, 1.0])
    sh = np.diag([1.0, 1.0, -1.0])
    S4 = _closure([_rot(z, math.pi / 2) @ sh])
    m["S4"] = _mol([("Si", (0, 0, 0))] + _orbit("C", (1.4, 0.4, 0.5), S4)
                   + _orbit("H", (2.0, 1.2, 1.1), S4))
    C3h = _closure([_rot(z, 2 * math.pi / 3), sh])
    m["B(OH)3_C3h"] = _mol([("B", (0, 0, 0))] + _orbit("O", (1.36, 0, 0), C3h)
                           + _orbit("H", (1.75, 0.9, 0), C3h))
    S6 = _closure([_rot(z, math.pi / 3) @ sh])
    m["S6"] = _mol([("Fe", (0, 0, 0))] + _orbit("C", (1.5, 0.3, 0.7), S6)
                   + _orbit("H", (2.1, 1.0, 1.2), S6))
    S8 = _closure([_rot(z, math.pi / 4) @ sh])
    m["S8"] = _mol(_orbit("S", (1.6, 0.3, 0.45), S8) + _orbit("O", (2.2, 1.1, 1.0), S8))
    C4h = _closure([_rot(z, math.pi / 2), sh])
    m["C4h"] = _mol([("Pt", (0, 0, 0))] + _orbit("N", (2.0, 0.6, 0), C4h)
                    + _orbit("H", (2.5, 1.2, 0.8), C4h))
    C5h = _closure([_rot(z, 2 * math.pi / 5), sh])
    m["C5h"] = _mol(_orbit("C", (1.2, 0, 0), C5h) + _orbit("H", (2.0, 0.7, 0), C5h))
    C6h = _closure([_rot(z, math.pi / 3), sh])
    m["C6h"] = _mol(_orbit("C", (1.4, 0, 0), C6h) + _orbit("Cl", (2.5, 0.9, 0), C6h))
    C3 = _closure([_rot(z, 2 * math.pi / 3)])
    m["C3"] = _mol([("P", (0, 0, 0.3))] + _orbit("C", (1.5, 0.4, -0.3), C3)
                   + _orbit("H", (2.0, 1.3, 0.2), C3))
    T = _closure([_rot(z, math.pi), _rot((1, 1, 1), 2 * math.pi / 3)])
    m["T"] = _mol([("C", (0, 0, 0))] + _orbit("H", (1.1, 0.35, 0.2), T))
    Th = _closure(T + [-np.eye(3)])
    m["Th"] = _mol([("Co", (0, 0, 0))] + _orbit("N", (1.9, 0, 0), Th)
                   + _orbit("O", (2.5, 0.9, 0), Th))
    ico = []
    for a, b in itertools.product((1, -1), repeat=2):
        ico += [(0, a, b * phi), (a, b * phi, 0), (b * phi, 0, a)]
    ico = np.array(ico) / np.linalg.norm(ico[0])
    m["B12H12_Ih"] = _mol([("B", tuple(1.7 * v)) for v in ico]
                          + [("H", tuple(2.9 * v)) for v in ico])
    # エチレン (D2h; 平面)、ホルムアルデヒド (C2v; B1 と B2 が両方現れる)
    m["C2H4"] = _mol([("C", (0, 0, 0.667)), ("C", (0, 0, -0.667)),
                      ("H", (0, 0.923, 1.232)), ("H", (0, -0.923, 1.232)),
                      ("H", (0, 0.923, -1.232)), ("H", (0, -0.923, -1.232))])
    m["H2CO"] = _mol([("C", (0, 0, 0)), ("O", (0, 0, 1.21)),
                      ("H", (0, 0.94, -0.59)), ("H", (0, -0.94, -0.59))])
    # 純回転二十面体群 I (キラル)
    I = _closure([_rot((0, 1, phi), 2 * math.pi / 5), _rot((0, 0, 1), math.pi)])
    m["I"] = _mol([("C", (0, 0, 0))] + _orbit("H", (1.119, -0.264, 1.637), I))
    # C60 (切頂二十面体; 頂点は (0, ±1, ±3φ), (±1, ±(2+φ), ±2φ), (±φ, ±2, ±(2φ+1)) の巡回置換)
    pts = set()
    for base in [(0, 1, 3 * phi), (1, 2 + phi, 2 * phi), (phi, 2, 2 * phi + 1)]:
        for s in itertools.product((1, -1), repeat=3):
            v = tuple(b * si for b, si in zip(base, s))
            for k in range(3):
                pts.add((v[k], v[(k + 1) % 3], v[(k + 2) % 3]))
    pts = np.array(sorted(pts))
    m["C60"] = _mol([("C", tuple(p)) for p in pts / np.linalg.norm(pts[0]) * 3.55])
    # PCl5 (D3h 三方両錐), XeF2 (D∞h), N2O (C∞v), フェロセン互い違い形 (D5d)
    m["PCl5"] = _mol([("P", (0, 0, 0)), ("Cl", (0, 0, 2.12)), ("Cl", (0, 0, -2.12))]
                     + _ring("Cl", 3, 2.02))
    m["XeF2"] = _mol([("Xe", (0, 0, 0)), ("F", (0, 0, 2.0)), ("F", (0, 0, -2.0))])
    m["N2O"] = _mol([("N", (0, 0, 0)), ("N", (0, 0, 1.13)), ("O", (0, 0, 2.32))])
    m["FeCp2_D5d"] = _mol([("Fe", (0, 0, 0))]
                          + _ring("C", 5, 1.22, 1.66) + _ring("C", 5, 1.22, -1.66, math.pi / 5)
                          + _ring("H", 5, 2.30, 1.66) + _ring("H", 5, 2.30, -1.66, math.pi / 5))
    return m


def distorted(cart, amplitude=0.03, seed=1):
    """実験・計算構造を模した小さなランダムな歪みを加える。"""
    rng = np.random.default_rng(seed)
    return cart + rng.uniform(-amplitude, amplitude, cart.shape)


# ── 結晶 ────────────────────────────────────────────────────────────────


def crystals():
    """名前 → pymatgen Structure (conventional cell)。"""
    from pymatgen.core import Lattice, Structure
    c = {}
    c["NaCl"] = Structure.from_spacegroup(
        "Fm-3m", Lattice.cubic(5.64), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    c["diamond"] = Structure.from_spacegroup(
        "Fd-3m", Lattice.cubic(3.567), ["C"], [[0, 0, 0]])
    c["rutile"] = Structure.from_spacegroup(
        "P4_2/mnm", Lattice.tetragonal(4.594, 2.959), ["Ti", "O"],
        [[0, 0, 0], [0.305, 0.305, 0]])
    c["ZnO"] = Structure.from_spacegroup(
        "P6_3mc", Lattice.hexagonal(3.25, 5.207), ["Zn", "O"],
        [[1/3, 2/3, 0], [1/3, 2/3, 0.382]])
    c["SrTiO3"] = Structure.from_spacegroup(
        "Pm-3m", Lattice.cubic(3.905), ["Sr", "Ti", "O"],
        [[0, 0, 0], [0.5, 0.5, 0.5], [0.5, 0.5, 0]])
    c["quartz"] = Structure.from_spacegroup(
        "P3_121", Lattice.hexagonal(4.913, 5.405), ["Si", "O"],
        [[0.4697, 0, 1/3], [0.4135, 0.2669, 0.1191]])
    c["calcite"] = Structure.from_spacegroup(
        "R-3c", Lattice.hexagonal(4.99, 17.06), ["Ca", "C", "O"],
        [[0, 0, 0], [0, 0, 0.25], [0.257, 0, 0.25]])
    return c
