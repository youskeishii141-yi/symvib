"""結晶因子群解析の回帰テスト。"""
from collections import Counter

import pytest

from structures import crystals
from symvib.crystal import analyze_crystal, graphite_2h, load_crystal, write_mode_files

# 音響モードを含む Γ 点の既約分解 (文献値)
EXPECTED = {
    "NaCl": ("Fm-3m", {"T1u": 2}),
    "diamond": ("Fd-3m", {"T1u": 1, "T2g": 1}),
    "rutile": ("P4_2/mnm", {"A1g": 1, "A2g": 1, "A2u": 2, "Eg": 1, "Eu": 4, "B1g": 1,
                            "B1u": 2, "B2g": 1}),
    "SrTiO3": ("Pm-3m", {"T1u": 4, "T2u": 1}),
    "quartz": ("P3_121", {"A1": 4, "A2": 5, "E": 9}),
    "calcite": ("R-3c", {"A1g": 1, "A1u": 2, "A2g": 3, "A2u": 4, "Eg": 4, "Eu": 6}),
    # 文献では 2B1 と書くことが多い (B1/B2 は σv/σd の取り方に依存)
    "ZnO": ("P6_3mc", {"A1": 2, "B1": 2, "E1": 2, "E2": 2}),
}


def _cif(tmp_path, name, structure):
    path = tmp_path / f"{name}.cif"
    structure.to(filename=str(path))
    return str(path)


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_literature_decomposition(tmp_path, name):
    r = analyze_crystal(*load_crystal(_cif(tmp_path, name, crystals()[name])))
    sg, expected = EXPECTED[name]
    assert r.space_group.startswith(sg + " ")
    assert dict(Counter(m.label for m in r.modes)) == expected
    assert sum(m.degeneracy for m in r.modes) == 3 * len(r.numbers)
    assert sum(m.degeneracy for m in r.modes if m.acoustic) == 3


@pytest.mark.parametrize("name", ["rutile", "SrTiO3", "calcite", "quartz"])
def test_mode_blocks_are_symmetry_invariant(tmp_path, name):
    """各モード (縮退ブロック) が力学表現で閉じている (= 既約表現の 1 組になっている)。"""
    import numpy as np
    from symvib.crystal import symmetry_adapted_blocks
    lat, pos, num, _ = load_crystal(_cif(tmp_path, name, crystals()[name]))
    r = analyze_crystal(lat, pos, num)
    rep = symmetry_adapted_blocks(lat, pos, num)["rep"]
    for m in r.modes:
        block = np.array([u.reshape(-1) for u in m.displacements])
        q, _ = np.linalg.qr(block.T)
        for D in rep:
            moved = D @ block.T
            assert np.allclose(q @ (q.T @ moved), moved, atol=1e-6), m.label


def test_graphite_demo(tmp_path):
    r = analyze_crystal(*graphite_2h())
    assert r.decomposition == "Γ = 2A2u + 2B2g + 2E2g + 2E1u"
    raman = [m for m in r.optical_modes if m.activity == "Raman"]
    assert [m.label for m in raman] == ["E2g", "E2g"]
    assert len(write_mode_files(r, tmp_path)) == 18          # 9 成分 × (xsf, vesta)


def test_oxidation_state_cif(tmp_path):
    """酸化数付き CIF (Ti4+, O2-) も元素記号として読める。"""
    st = crystals()["rutile"].copy()
    st.add_oxidation_state_by_element({"Ti": 4, "O": -2})
    lat, pos, num, sym = load_crystal(_cif(tmp_path, "rutile_ox", st))
    assert set(sym) == {"Ti", "O"}
    assert analyze_crystal(lat, pos, num).space_group.startswith("P4_2/mnm")


def _degenerate(result):
    return [m for m in result.modes if m.degeneracy > 1]


@pytest.mark.parametrize("name", ["quartz", "calcite", "rutile", "ZnO", "SrTiO3"])
def test_degenerate_components_do_not_depend_on_atom_order(name):
    """縮退成分の選び方は物質ごとの設定なしの共通規則で、CIF の原子の順番にもよらない。"""
    import numpy as np
    s = crystals()[name]
    lat, frac = s.lattice.matrix, s.frac_coords
    z = np.array([site.specie.Z for site in s])
    perm = np.random.default_rng(7).permutation(len(z))
    r1 = analyze_crystal(lat, frac, z)
    r2 = analyze_crystal(lat, frac[perm], z[perm])

    def keyed(r, u):                          # 座標で並べ直した変位 (順番によらない比較用)
        f = np.round(np.linalg.solve(r.lattice.T, r.cart.T).T % 1.0, 3) % 1.0
        order = sorted(range(len(f)), key=lambda i: (r.numbers[i], *f[i]))
        return u[order].ravel()

    for m1, m2 in zip(_degenerate(r1), _degenerate(r2)):
        assert m1.label == m2.label
        for u, v in zip(m1.displacements, m2.displacements):
            a, b = keyed(r1, u), keyed(r2, v)
            assert np.dot(a, b) / np.linalg.norm(a) / np.linalg.norm(b) > 0.999


@pytest.mark.parametrize("name", ["graphite", "ZnO", "NaCl", "SrTiO3"])
def test_components_follow_crystal_axes_when_possible(name):
    """軸方向に動けるモード (六方晶の E1/E1u、立方晶の T1u) では、成分 k が軸 k に平行。"""
    import numpy as np
    if name == "graphite":
        lat, frac, z, _ = graphite_2h()
    else:
        s = crystals()[name]
        lat, frac, z = s.lattice.matrix, s.frac_coords, [site.specie.Z for site in s]
    r = analyze_crystal(lat, frac, z)
    axes = np.asarray(lat, float) / np.linalg.norm(lat, axis=1, keepdims=True)
    polar = [m for m in _degenerate(r) if m.label.startswith(("E1", "T1u"))]
    assert polar
    for m in polar:
        for k, u in enumerate(m.displacements):
            moving = np.linalg.norm(u, axis=1) > 1e-6
            cos = np.abs(u[moving] @ axes[k]) / np.linalg.norm(u[moving], axis=1)
            assert np.all(cos > 0.999)


def _corundum():
    from pymatgen.core import Lattice, Structure
    s = Structure.from_spacegroup("R-3c", Lattice.hexagonal(4.759, 12.991), ["Al", "O"],
                                  [[0, 0, 0.3522], [0.3062, 0, 0.25]])
    return s.lattice.matrix, s.frac_coords, [site.specie.Z for site in s]


def test_bonds_follow_vesta_defaults_no_metal_metal():
    """結合は VESTA 既定の元素対の表で判定: コランダムは Al–O だけ (Al–Al を引かない)。"""
    import numpy as np
    from symvib.common import bond_pairs
    from symvib.crystal import cell_atoms
    r = analyze_crystal(*_corundum())
    cart, _, parent = cell_atoms(r)
    numbers = [r.numbers[p] for p in parent]
    pairs = bond_pairs(cart, numbers)
    kinds = {tuple(sorted((numbers[i], numbers[j]))) for i, j in pairs}
    assert kinds == {(8, 13)}
    assert max(np.linalg.norm(cart[i] - cart[j]) for i, j in pairs) < 2.11


def test_modes_list_involved_sites_and_formats_are_selectable(tmp_path):
    """関与サイトは元素(Wyckoff) の名前で、書き出す形式は選べる。"""
    r = analyze_crystal(*_corundum())
    assert r.decomposition == "Γ = 2A1g + 2A1u + 3A2g + 3A2u + 5Eg + 5Eu"
    assert all(m.contributions == ["Al(c)", "O(e)"] for m in r.modes if m.label == "Eg")
    vesta = write_mode_files(r, tmp_path / "v", formats=("vesta",))
    xsf = write_mode_files(r, tmp_path / "x", formats=("xsf",))
    assert vesta and all(p.endswith(".vesta") for p in vesta)
    assert xsf and all(p.endswith(".xsf") for p in xsf)
    text = open(vesta[0], encoding="utf-8").read()
    assert " Al     O " in text and " Al    Al " not in text          # SBOND は Al–O だけ


def test_primitive_and_rhombohedral_inputs_work():
    """primitive cell・菱面体軸・超格子でも同じ分解とラベル (標準の向き) になる。"""
    from pymatgen.core import Lattice, Structure
    nacl = Structure.from_spacegroup("Fm-3m", Lattice.cubic(5.64), ["Na", "Cl"],
                                     [[0, 0, 0], [.5, .5, .5]])
    for s in (nacl.get_primitive_structure(), nacl * (2, 1, 1)):
        r = analyze_crystal(s.lattice.matrix, s.frac_coords, [x.specie.Z for x in s])
        assert r.decomposition == "Γ = 2T1u"
    lat, frac, z = _corundum()
    cor = Structure(lat, z, frac).get_primitive_structure()        # 菱面体軸
    r = analyze_crystal(cor.lattice.matrix, cor.frac_coords, [x.specie.Z for x in cor])
    assert r.decomposition == "Γ = 2A1g + 2A1u + 3A2g + 3A2u + 5Eg + 5Eu"
    assert all(set(m.tensor) == {"xx", "yy", "zz"} for m in r.modes if m.label == "A1g")


def test_orthorhombic_labels_follow_bilbao():
    """D2h: B2g = xz, B3g = yz。mm2: B1 = xz (Bilbao の標準の設定)。"""
    from pymatgen.core import Lattice, Structure
    s = Structure.from_spacegroup("Pnma", Lattice.orthorhombic(5.61, 7.67, 5.36),
                                  ["Gd", "Fe", "O", "O"],
                                  [[0.06, 0.25, 0.985], [0, 0, 0.5], [0.47, 0.25, 0.1],
                                   [0.3, 0.05, 0.7]])
    r = analyze_crystal(s.lattice.matrix, s.frac_coords, [x.specie.Z for x in s])
    tensors = {m.label: m.tensor for m in r.modes}
    assert tensors["B1g"] == ["xy"] and tensors["B2g"] == ["xz"] and tensors["B3g"] == ["yz"]
    for sg in ("Pmm2", "Pna2_1", "Cmc2_1"):
        s = Structure.from_spacegroup(sg, Lattice.orthorhombic(5.1, 6.3, 7.2), ["Ti", "O"],
                                      [[0.11, 0.23, 0.31], [0.37, 0.13, 0.07]])
        r = analyze_crystal(s.lattice.matrix, s.frac_coords, [x.specie.Z for x in s])
        assert {m.label: m.tensor for m in r.modes}["B1"] == ["xz"], sg


@pytest.mark.parametrize("name", ["NaCl", "ZnO", "SrTiO3", "calcite"])
def test_optical_modes_keep_the_centre_of_mass_fixed(name):
    """Γ 点の光学モードは Σ M_i u_i = 0 (重心が動かない)。"""
    import numpy as np
    from symvib.common import atomic_mass
    s = crystals()[name]
    r = analyze_crystal(s.lattice.matrix, s.frac_coords, [x.specie.Z for x in s])
    masses = np.array([atomic_mass(z) for z in r.numbers])
    for m in r.optical_modes:
        for u in m.displacements:
            assert np.linalg.norm((masses[:, None] * u).sum(0)) < 1e-6 * masses.sum()
