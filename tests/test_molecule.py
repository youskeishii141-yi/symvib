"""分子振動解析の回帰テスト (文献値との一致と群論的な整合性)。"""
from collections import Counter

import numpy as np
import pytest

from structures import distorted, molecules
from symvib.common import atomic_mass
from symvib.molecule import (analyze_molecule, find_point_group_operations,
                                         write_mode_files)

MOLECULES = molecules()

# 文献値 (Mulliken の軸の取り方; 平面分子は分子面の法線が x)
EXPECTED = {
    "H2O": ("C2v", {"A1": 2, "B2": 1}),
    "H2CO": ("C2v", {"A1": 3, "B1": 1, "B2": 2}),
    "NH3": ("C3v", {"A1": 2, "E": 2}),
    "CH4": ("Td", {"A1": 1, "E": 1, "T2": 2}),
    "CO2": ("D*h", {"Σg⁺": 1, "Σu⁺": 1, "Πu": 1}),
    "HCN": ("C*v", {"Σ⁺": 2, "Π": 1}),
    "C2H2": ("D*h", {"Σg⁺": 2, "Σu⁺": 1, "Πg": 1, "Πu": 1}),
    "C2H4": ("D2h", {"Ag": 3, "Au": 1, "B1u": 2, "B2u": 2, "B3u": 1, "B2g": 1, "B3g": 2}),
    "C6H6": ("D6h", {"A1g": 2, "A2g": 1, "B2g": 2, "E1g": 1, "E2g": 4, "A2u": 1,
                     "B1u": 2, "B2u": 2, "E1u": 3, "E2u": 2}),
    "SF6": ("Oh", {"A1g": 1, "Eg": 1, "T1u": 2, "T2g": 1, "T2u": 1}),
    "BF3": ("D3h", {"A1'": 1, "E'": 2, "A2''": 1}),
    "XeF4": ("D4h", {"A1g": 1, "B1g": 1, "B2g": 1, "A2u": 1, "B2u": 1, "Eu": 2}),
    "C2H6": ("D3d", {"A1g": 3, "Eg": 3, "A1u": 1, "A2u": 2, "Eu": 3}),
    "C3H4": ("D2d", {"A1": 3, "B1": 1, "B2": 3, "E": 4}),
    "FeCp2_D5h": ("D5h", {"A1'": 4, "A2'": 1, "E1'": 6, "E2'": 6, "A1''": 2, "A2''": 4,
                          "E1''": 5, "E2''": 6}),
    "B(OH)3_C3h": ("C3h", {"A'": 3, "E'": 4, "A''": 2, "E''": 1}),
    "B12H12_Ih": ("Ih", {"Ag": 2, "T1g": 1, "Gg": 2, "Hg": 4, "T1u": 3, "T2u": 2,
                         "Gu": 2, "Hu": 2}),
    # C60: Dresselhaus, Dresselhaus, Eklund, Science of Fullerenes (1996) Table 11.1
    "C60": ("Ih", {"Ag": 2, "T1g": 3, "T2g": 4, "Gg": 6, "Hg": 8, "Au": 1, "T1u": 4, "T2u": 5,
                   "Gu": 6, "Hu": 7}),
    # PCl5 (Herzberg II): 2A1' + 3E' + 2A2'' + E'';  XeF2 / N2O: 直線 3 原子分子
    "PCl5": ("D3h", {"A1'": 2, "E'": 3, "A2''": 2, "E''": 1}),
    "XeF2": ("D*h", {"Σg⁺": 1, "Σu⁺": 1, "Πu": 1}),
    "N2O": ("C*v", {"Σ⁺": 2, "Π": 1}),
    # フェロセン互い違い形 (D5d; Nakamoto): 重なり形 (D5h) との相関 A'→Ag, E1'→E1u, E1''→E1g など
    "FeCp2_D5d": ("D5d", {"A1g": 4, "A2g": 1, "E1g": 5, "E2g": 6, "A1u": 2, "A2u": 4,
                          "E1u": 6, "E2u": 6}),
}


def _counts(result):
    return dict(Counter(m.label for m in result.modes))


def _mass_weighted(mode, sm):
    return np.array([(u * sm[:, None]).reshape(-1) for u in mode.displacements])


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_literature_decomposition(name):
    cart, Z, sym = MOLECULES[name]
    result = analyze_molecule(cart, Z, sym)
    pg, expected = EXPECTED[name]
    assert result.point_group == pg
    assert _counts(result) == expected


@pytest.mark.parametrize("name", sorted(MOLECULES))
def test_modes_are_symmetry_invariant_vibrations(name):
    """全モード: 自由度が 3N-6 (直線 3N-5)、並進を含まず、互いに直交し、
    縮退ブロックが対称操作で閉じている (= 既約表現の部分空間になっている)。"""
    cart, Z, sym = MOLECULES[name]
    r = analyze_molecule(cart, Z, sym)
    na = len(Z)
    assert sum(m.degeneracy for m in r.modes) == 3 * na - (5 if r.is_linear else 6)
    sm = np.sqrt([atomic_mass(z) for z in Z])
    for m in r.modes:
        for u in m.displacements:                       # 運動量 Σ m u = 0
            assert np.allclose((u * (sm ** 2)[:, None]).sum(0), 0, atol=1e-8)
    # 異なるモードどうしは質量重みで直交、各モードの成分は線形独立 (縮退成分は軸にそろえる
    # ので互いに直交とは限らない)
    blocks = [np.linalg.qr(_mass_weighted(m, sm).T)[0].T for m in r.modes]
    for m, blk in zip(r.modes, blocks):
        assert np.linalg.matrix_rank(_mass_weighted(m, sm), tol=1e-8) == m.degeneracy
    rows = np.vstack(blocks)
    assert np.allclose(rows @ rows.T, np.eye(len(rows)), atol=1e-6)
    if r.is_linear:
        return
    ops = find_point_group_operations(r.cart, Z, 1e-2)
    for m in r.modes:
        block = _mass_weighted(m, sm)
        q, _ = np.linalg.qr(block.T)
        for R, perm in ops:
            moved = np.zeros_like(block)
            for i, k in enumerate(perm):
                moved[:, 3*k:3*k+3] = block[:, 3*i:3*i+3] @ R.T
            assert np.allclose(q @ (q.T @ moved.T), moved.T, atol=1e-6), m.label


@pytest.mark.parametrize("name", ["CH4", "C6H6", "SF6", "C2H4", "CO2"])
def test_noisy_geometry_is_symmetrized(name):
    """小さな歪み (実験・計算構造の誤差) を含んでも同じ結果になる。"""
    cart, Z, sym = MOLECULES[name]
    clean = analyze_molecule(cart, Z, sym)
    noisy = analyze_molecule(distorted(cart, 0.03), Z, sym)
    assert noisy.point_group == clean.point_group
    assert _counts(noisy) == _counts(clean)
    assert any("対称化" in n for n in noisy.notes)


def test_distorted_molecule_with_atom_near_centre_of_mass():
    """重心のすぐ近くに原子 (N) がある歪んだ C3 分子。以前は原点に近い原子を基準に回転を作って
    いたため、0.02 Å の歪みで対称操作を見失い C1 になっていた (一般位置の軌道から作った分子で発見)。"""
    cart = np.array([[0.0122, 0.0123, 0.0006], [-0.7672, -1.146, 0.4846], [0.6091, 0.8482, 0.9625],
                     [-0.6471, -1.2773, 1.8678], [-1.4581, -0.0971, -1.8519],
                     [-0.4098, 0.5937, -1.2455], [0.8327, 2.1476, 0.5031]])
    r = analyze_molecule(cart, [7, 6, 6, 1, 1, 6, 1], None)
    assert r.point_group == "C3" and _counts(r) == {"A": 5, "E": 5}


def test_strict_tolerance_keeps_lower_symmetry():
    cart, Z, sym = MOLECULES["CH4"]
    r = analyze_molecule(distorted(cart, 0.03), Z, sym, tol=1e-3)
    assert r.point_group == "C1"


def test_planar_out_of_plane_guess_is_flagged():
    cart, Z, sym = MOLECULES["C6H6"]
    r = analyze_molecule(cart, Z, sym)
    assert any("任意" in n for n in r.notes)


def test_xsf_output_is_utf8_and_ascii_named(tmp_path):
    cart, Z, sym = MOLECULES["CO2"]
    r = analyze_molecule(cart, Z, sym)
    paths = write_mode_files(r, tmp_path)
    assert len(paths) == 8                      # Σg⁺, Σu⁺, Πu (2 成分) × (xsf, vesta)
    for p in paths:
        assert str(p).isascii()
        with open(p, encoding="utf-8") as f:
            text = f.read()
        assert ("PRIMCOORD" if p.endswith(".xsf") else "VECTR") in text
        assert p.endswith(".xsf") or text.isascii()      # VESTA 形式は ASCII のみ


def test_read_xyz(tmp_path):
    from symvib.molecule import load_molecule
    path = tmp_path / "water.xyz"
    path.write_text("3\nwater\nO 0 0 0.117\nH1 0 0.757 -0.467\n1 0 -0.757 -0.467\n",
                    encoding="utf-8")
    cart, Z, sym = load_molecule(str(path))
    assert Z == [8, 1, 1] and sym == ["O", "H", "H"]
    assert analyze_molecule(cart, Z, sym).decomposition == "Γvib = 2A1 + B2"


def test_bad_xyz_reports_clear_error(tmp_path):
    from symvib.common import AnalysisError
    from symvib.molecule import load_molecule
    path = tmp_path / "bad.xyz"
    path.write_text("5\ncomment\nO 0 0 0\n", encoding="utf-8")
    with pytest.raises(AnalysisError, match="xyz"):
        load_molecule(str(path))


@pytest.mark.parametrize("name", ["CH4", "C6H6", "C5h", "C6h", "SF6", "FeCp2_D5h"])
def test_degenerate_components_do_not_depend_on_atom_order(name):
    """縮退成分 (と、網目で決まらない面外モードの分け方) は原子の順番によらない。"""
    cart, Z, sym = molecules()[name]
    cart, Z = np.asarray(cart, float), np.asarray(Z)
    perm = np.random.default_rng(11).permutation(len(Z))
    r1 = analyze_molecule(cart, list(Z), list(sym))
    r2 = analyze_molecule(cart[perm], list(Z[perm]), [sym[i] for i in perm])
    inv = np.argsort(perm)
    for m1, m2 in zip(r1.modes, r2.modes):
        assert m1.label == m2.label
        for u, v in zip(m1.displacements, m2.displacements):
            a, b = u.ravel(), v[inv].ravel()
            assert np.dot(a, b) / np.linalg.norm(a) / np.linalg.norm(b) > 0.999


def test_molecule_bonds_fall_back_to_covalent_radii():
    """分子では VESTA の表にない組 (H–H) も共有結合半径で結ぶ。"""
    from symvib.common import bond_pairs
    assert bond_pairs([[0, 0, 0], [0, 0, 0.74]], [1, 1], molecule=True) == [(0, 1)]
    cart, Z, _ = molecules()["C6H6"]
    pairs = bond_pairs(cart, Z, molecule=True)
    assert len(pairs) == 12                                 # C–C 6 本 + C–H 6 本
