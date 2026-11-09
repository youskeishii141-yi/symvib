"""コードの内部経路に頼らない独立検査 (tests/invariants.py) を、代表的な結晶・分子に対して回す。

全 230 空間群 (530 設定) の網羅は tests/sweep_space_groups.py (数分かかるので pytest から分離)。
ここでは各結晶系の代表と、Bilbao の点群表との照合が設定に依存する群 (-42m / -4m2, 321 / 312,
3m1 / 31m, -6m2 / -62m, mm2 の向き) を通常のテストとして残す。
"""
import numpy as np
import pytest
from pymatgen.core import Lattice, Structure

from invariants import bilbao_available, check_crystal, check_molecule
from structures import crystals, molecules
from sweep_space_groups import synthetic_structure
from symvib.crystal import analyze_crystal
from symvib.molecule import analyze_molecule

# 各結晶系・設定の代表 (Hall 番号): P1, P-1, P2_1/c (b), C2/c (c-unique), Pnma, Pbnm, Cmcm,
# P4/m, P-42m, P-4m2, I4_1/amd (origin 1, 2), P3_121, P3_112, P3m1, P31m, R-3c (H, R),
# P6_3/mmc, P-6m2, P-62m, Pa-3, Fd-3m (origin 1, 2), Pm-3m
REPRESENTATIVE_HALLS = [1, 2, 81, 97, 292, 294, 280, 357, 388, 392, 426, 427, 442, 443, 447,
                        448, 460, 461, 488, 481, 483, 501, 525, 526, 517]


@pytest.mark.parametrize("hall", REPRESENTATIVE_HALLS)
def test_synthetic_structure_passes_independent_checks(hall):
    lat, frac, num = synthetic_structure(hall)
    r = analyze_crystal(lat, frac, num)
    probs, info = check_crystal(r, lat, frac, num)
    # 設定の向きが Bilbao の点群表と異なる群 (-4m2, 312, 31m, -62m) では、B1/B2 や
    # 成分の割り当てが「2 回軸に対して対称なものを 1」という規則で決まる (Mulliken)。
    # その違いだけが報告されることを確かめ、それ以外は一切の不一致を認めない
    orientation_note = any("not along the Bilbao" in n for n in info.get("bilbao_notes", []))
    if orientation_note:
        probs = [p for p in probs if "Bilbao label" not in p and "vs Bilbao" not in p]
    assert probs == []


@pytest.mark.parametrize("name", sorted(crystals()))
def test_literature_crystals_pass_independent_checks(name):
    s = crystals()[name]
    lat, frac, num = s.lattice.matrix, s.frac_coords, [x.specie.Z for x in s]
    r = analyze_crystal(lat, frac, num)
    probs, info = check_crystal(r, lat, frac, num)
    assert probs == []
    if bilbao_available():
        assert all(v == k for k, v in info["bilbao_labels"].items())


@pytest.mark.parametrize("name", sorted(molecules()))
def test_molecules_pass_independent_checks(name):
    cart, Z, sym = molecules()[name]
    assert check_molecule(analyze_molecule(cart, Z, sym)) == []


def test_mulliken_labels_do_not_depend_on_operation_order():
    """D2d (-42m) など C2 軸が複数ある群で、操作の並び順を変えてもラベルが変わらない。"""
    import spglib
    from spgrep import get_spacegroup_irreps
    from symvib.crystal import (merge_complex_pairs, mulliken_labels,
                                            reduce_irrep_to_pg, symmetry_adapted_blocks)
    for sg in ("P-42m", "P-4m2", "P-4", "P4_2/mnm", "Pmmm", "P6_3/mmc", "P-62m"):
        lat = (Lattice.tetragonal(5.0, 7.0) if "4" in sg else Lattice.hexagonal(5.0, 7.0)
               if "6" in sg else Lattice.orthorhombic(5.0, 6.0, 7.0))
        s = Structure.from_spacegroup(sg, lat, ["Ti", "O"], [[0.11, 0.23, 0.31], [0.5, 0.5, 0.0]])
        cell = (s.lattice.matrix, s.frac_coords, [x.specie.Z for x in s])
        sa = symmetry_adapted_blocks(*cell)
        conv = spglib.standardize_cell(cell, to_primitive=False, no_idealize=True, symprec=1e-3)
        irr_conv, rots_conv, _, _ = get_spacegroup_irreps(np.array(conv[0]), np.array(conv[1]),
                                                          [int(z) for z in conv[2]], [0, 0, 0])
        irr = merge_complex_pairs([reduce_irrep_to_pg(i, rots_conv, np.array(conv[0]),
                                                      sa["rots"], sa["lattice"]) for i in irr_conv])
        ref = mulliken_labels(irr, sa["R_std"])
        for seed in range(4):
            perm = np.random.default_rng(seed).permutation(len(sa["R_std"]))
            got = mulliken_labels([i[perm] for i in irr], sa["R_std"][perm])
            assert got == ref, (sg, ref, got)


def test_raman_tensor_forms_match_loudon():
    """ラマンテンソルの形 (Loudon 1964 の表) — 独立パラメータが 2 つある 4/m, -3, 6/m を含む。"""
    def forms(sg, lat, species, coords):
        s = Structure.from_spacegroup(sg, lat, species, coords)
        r = analyze_crystal(s.lattice.matrix, s.frac_coords, [x.specie.Z for x in s])
        out = {}
        for m in r.optical_modes:
            if m.raman_tensors and m.label not in out:
                out[m.label] = m.raman_tensors
        return out

    f = forms("P4/m", Lattice.tetragonal(5.0, 7.0), ["Ti", "O"], [[0.11, 0.23, 0.31], [0.5, 0.5, 0]])
    assert f["Ag"] == ["[[a,·,·],[·,a,·],[·,·,b]]"]
    assert f["Bg"] == ["[[a,b,·],[b,-a,·],[·,·,·]]"]
    assert f["Eg"] == ["[[·,·,a],[·,·,b],[a,b,·]]"] * 2          # 2 つの独立パラメータ
    f = forms("P4/mmm", Lattice.tetragonal(5.0, 7.0), ["Ti", "O"], [[0.11, 0.23, 0.31], [0.5, 0.5, 0]])
    assert f["B1g"] == ["[[a,·,·],[·,-a,·],[·,·,·]]"] and f["B2g"] == ["[[·,a,·],[a,·,·],[·,·,·]]"]
    assert f["Eg"] == ["[[·,·,a],[·,·,·],[a,·,·]]", "[[·,·,·],[·,·,a],[·,a,·]]"]
    f = forms("P-3", Lattice.hexagonal(5.0, 7.0), ["Ti", "O"], [[0.11, 0.23, 0.31], [0, 0, 0]])
    assert f["Eg"] == ["[[a,b,c],[b,-a,d],[c,d,·]]"] * 2         # 4 つの独立パラメータ
    f = forms("P321", Lattice.hexagonal(5.0, 7.0), ["Si", "O"], [[0.47, 0, 0], [0.41, 0.27, 0.12]])
    assert set(f["E"]) == {"[[a,·,·],[·,-a,b],[·,b,·]]", "[[·,a,b],[a,·,·],[b,·,·]]"}
    f = forms("P6/m", Lattice.hexagonal(5.0, 7.0), ["Ti", "O"], [[0.11, 0.23, 0.31], [0, 0, 0]])
    assert f["E1g"] == ["[[·,·,a],[·,·,b],[a,b,·]]"] * 2 and f["E2g"] == ["[[a,b,·],[b,-a,·],[·,·,·]]"] * 2
    f = forms("P2_1/c", Lattice.monoclinic(5.1, 5.2, 5.3, 99.0), ["Zr", "O", "O"],
              [[0.275, 0.04, 0.208], [0.07, 0.33, 0.345], [0.45, 0.76, 0.48]])
    assert f["Ag"] == ["[[a,·,d],[·,b,·],[d,·,c]]"] and f["Bg"] == ["[[·,a,·],[a,·,b],[·,b,·]]"]
    f = forms("Pa-3", Lattice.cubic(5.4), ["Fe", "S"], [[0, 0, 0], [0.385, 0.385, 0.385]])
    assert f["Ag"] == ["[[a,·,·],[·,a,·],[·,·,a]]"]
    assert f["Eg"] == ["[[a,·,·],[·,b,·],[·,·,-a-b]]"] * 2      # 複素ペア → 対角トレースレス全体
    assert sorted(f["Tg"]) == sorted(["[[·,a,·],[a,·,·],[·,·,·]]", "[[·,·,·],[·,·,a],[·,a,·]]",
                                      "[[·,·,a],[·,·,·],[a,·,·]]"])
