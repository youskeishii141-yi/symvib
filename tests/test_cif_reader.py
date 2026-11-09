"""CIF 読み込みの回帰テスト (実在の CIF によくある「くせ」)。"""
from collections import Counter

import numpy as np
import pytest

from symvib.cif_reader import read_cif
from symvib.common import AnalysisError, atomic_mass, element_symbol, is_pseudo
from symvib.crystal import analyze_crystal, cell_atoms, load_crystal, write_mode_files

QUARTZ = """data_8103513
_symmetry_space_group_name_H-M   'P 31 2 1'
_cell_angle_alpha                90
_cell_angle_beta                 90
_cell_angle_gamma                120
_cell_length_a                   4.583(2)
_cell_length_b                   4.583(2)
_cell_length_c                   5.189(3)
loop_
_symmetry_equiv_pos_as_xyz
x,y,z
-y,x-y,z+1/3
-x+y,-x,z+2/3
y,x,-z
-x,-x+y,-z+1/3
x-y,-y,-z+2/3
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
O1 O-2 0.4006(3) 0.3052(3) 0.2423(2) 1
Si1 Si+4 0.4476(2) 0 0.3333 1
"""


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text, encoding="utf-8")
    return str(p)


def test_quartz_with_uncertainties_and_charged_types(tmp_path):
    c = read_cif(_write(tmp_path, "q.cif", QUARTZ))
    assert Counter(c.names) == {"O": 6, "Si": 3}
    r = analyze_crystal(c.lattice, c.frac, c.numbers)
    assert r.space_group.startswith("P3_121")
    assert r.decomposition == "Γ = 4A1 + 5A2 + 9E"


def test_first_of_multiple_blocks_and_labels_only(tmp_path):
    """2 つの構造を含む CIF は最初を使う。型の記載がなければラベル (Na1, CL1) から元素を推定。"""
    nacl = """data_first
_cell_length_a 5.64
_cell_length_b 5.64
_cell_length_c 5.64
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'F m -3 m'
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Na1 0 0 0
CL1 0.5 0.5 0.5
"""
    text = "data_global\n_publ_section_title 'two structures'\n\n" + nacl + \
           "\n" + QUARTZ.replace("data_8103513", "data_second")
    c = read_cif(_write(tmp_path, "multi.cif", text))
    assert c.block == "first"
    assert Counter(c.names) == {"Na": 4, "Cl": 4}            # 対称操作を空間群名から補う
    assert any("空間群" in n for n in c.notes)
    r = analyze_crystal(c.lattice, c.frac, c.numbers)
    assert r.decomposition == "Γ = 2T1u"


def test_partial_occupancy_merged_into_pseudo_atom(tmp_path):
    """同じ座標に Ti1 と V2 → 1 つの仮想原子 'Ti1V2' (質量は占有率で平均)。"""
    text = """data_rutile_like
_cell_length_a 4.594
_cell_length_b 4.594
_cell_length_c 2.959
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'P 42/m n m'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Ti1 Ti 0 0 0 0.6
V2 V 0 0 0 0.4
O1 O 0.305 0.305 0 1.0
"""
    path = _write(tmp_path, "tiv.cif", text)
    c = read_cif(path)
    pseudo = [z for z in c.numbers if is_pseudo(z)]
    assert len(pseudo) == 2 and element_symbol(pseudo[0]) == "Ti1V2"
    assert atomic_mass(pseudo[0]) == pytest.approx(0.6 * 47.867 + 0.4 * 50.9415, rel=1e-3)
    assert any("Ti1V2" in n for n in c.notes)
    notes = []
    r = analyze_crystal(*load_crystal(path, notes), notes=notes)
    assert r.space_group.startswith("P4_2/mnm")
    assert dict(Counter(m.label for m in r.modes))["Eu"] == 4   # ルチルと同じ分解
    assert any("Ti1V2" in n for n in r.notes)
    paths = write_mode_files(r, tmp_path / "out")            # 仮想原子でも書き出せる
    assert any(p.endswith(".vesta") for p in paths)


def test_unknown_element_becomes_pseudo_atom(tmp_path):
    text = QUARTZ.replace("Si1 Si+4", "Q1 Qx")
    c = read_cif(_write(tmp_path, "unk.cif", text))
    assert any("Q1" in n for n in c.notes)
    q = [z for z in c.numbers if is_pseudo(z)]
    assert len(q) == 3 and atomic_mass(q[0]) == 1.0
    assert analyze_crystal(c.lattice, c.frac, c.numbers).space_group.startswith("P3_121")


def test_missing_coordinates_is_clear_error(tmp_path):
    with pytest.raises(AnalysisError, match="分率座標"):
        read_cif(_write(tmp_path, "bad.cif", "data_x\n_cell_length_a 1\n"))


def test_cell_atoms_include_boundary_images(tmp_path):
    """表示用の原子: 境界上の像を含む。NaCl の慣用胞では 14 Na + 13 Cl = 27 個。"""
    from structures import crystals
    p = str(tmp_path / "nacl.cif")
    crystals()["NaCl"].to(filename=p)
    r = analyze_crystal(*load_crystal(p))
    cart, frac, parent = cell_atoms(r, boundary=True)
    names = Counter(r.symbols[i] for i in parent)
    assert names == {"Na": 14, "Cl": 13}
    _, frac_in, parent_in = cell_atoms(r, boundary=False)
    assert len(parent_in) == 8 and np.all((frac_in >= 0) & (frac_in < 1))


def _cif(tmp_path, name, body):
    p = tmp_path / name
    p.write_text(body, encoding="utf-8-sig")              # BOM 付きでも読めること
    return str(p)


_DIAMOND_O2 = """data_d
_cell_length_a 3.567
_cell_length_b 3.567
_cell_length_c 3.567
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'F d -3 m'
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
C1 1/8 1/8 1/8
"""


def test_space_group_symbol_without_symops(tmp_path):
    """対称操作の記載がなくても、記号から設定 (原点など) を選んで正しく展開する。"""
    from symvib.cif_reader import read_cif
    c = read_cif(_cif(tmp_path, "d.cif", _DIAMOND_O2))          # 原点 2 の座標・分数表記
    assert len(c.numbers) == 8
    p21n = _DIAMOND_O2.replace("'F d -3 m'", "'P 21/n'").replace(
        "3.567\n_cell_length_b 3.567\n_cell_length_c 3.567", "5.1\n_cell_length_b 6.2\n"
        "_cell_length_c 7.3").replace("_cell_angle_beta 90", "_cell_angle_beta 101.2").replace(
        "C1 1/8 1/8 1/8", "C1 0.11 0.21 0.31")
    c = read_cif(_cif(tmp_path, "p.cif", p21n))
    assert len(c.numbers) == 4
    pos = {tuple(np.round(f, 3)) for f in c.frac}
    assert (0.61, 0.29, 0.81) in pos                          # n 映進 (x+1/2, -y+1/2, z+1/2)


@pytest.mark.parametrize("hm, angles, expected", [
    ("P 21/c", (90, 101, 90), "P 1 2_1/c 1"),            # β ≠ 90 → unique axis b
    ("P 21/c", (90, 90, 101), "P 1 1 2_1/a"),            # γ ≠ 90 → unique axis c
    ("P 21/c", (101, 90, 90), "P 2_1/b 1 1"),            # α ≠ 90 → unique axis a
    ("C 2/c", (90, 90, 101), "A 1 1 2/a"),
    ("P 1 1 21/b", (90, 90, 101), "P 1 1 2_1/b"),        # 完全な記号はそのまま
])
def test_monoclinic_unique_axis_follows_the_cell_angles(tmp_path, hm, angles, expected):
    """短い記号 ("P 21/c") だけでは unique axis が決まらないので、90° でない角から選ぶ。"""
    text = _DIAMOND_O2.replace("'F d -3 m'", f"'{hm}'").replace(
        "_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90",
        f"_cell_angle_alpha {angles[0]}\n_cell_angle_beta {angles[1]}\n_cell_angle_gamma {angles[2]}"
    ).replace("3.567\n_cell_length_b 3.567\n_cell_length_c 3.567", "5.1\n_cell_length_b 6.2\n"
              "_cell_length_c 7.3").replace("C1 1/8 1/8 1/8", "C1 0.11 0.21 0.31")
    c = read_cif(_cif(tmp_path, "m.cif", text))
    assert len(c.numbers) == (8 if hm.startswith("C") else 4)
    assert any(expected in str(n) for n in c.notes), [str(n) for n in c.notes]


def test_element_guessing_is_conservative():
    from symvib.cif_reader import _element_from
    assert _element_from("Wat")[0] is None                    # W (タングステン) ではない
    assert _element_from("Uiso")[0] is None
    assert _element_from("NH4")[0] == "N"                     # Nh (ニホニウム) ではない
    assert _element_from("SI1")[0] == "Si"
    assert _element_from("Ti4+")[0] == "Ti"
    assert _element_from("D")[0] == "H"


_GRAPHITE_ROUNDED = """data_g
_cell_length_a 2.46
_cell_length_b 2.46
_cell_length_c 6.71
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 120
_symmetry_space_group_name_H-M 'P 63/m m c'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
C1 C 0 0 1/4
C2 C {x} {y} 0.25
"""


@pytest.mark.parametrize("x, y", [("1/3", "2/3"), ("0.333", "0.667"), ("0.3333", "0.6667"),
                                  ("-1/3", "-2/3")])
def test_fractional_coordinates_like_one_third(tmp_path, x, y):
    """1/3 などの分数や、それを丸めた小数 (0.333) でも対称性が下がらない。"""
    from symvib.crystal import analyze_crystal, load_crystal
    path = _cif(tmp_path, "g.cif", _GRAPHITE_ROUNDED.format(x=x, y=y))
    r = analyze_crystal(*load_crystal(path))
    assert r.space_group.startswith("P6_3/mmc")


def test_refined_parameters_are_not_rounded_to_fractions():
    from symvib.cif_reader import _coordinate
    assert _coordinate("0.333") == (1 / 3, True)
    assert _coordinate("0.3330")[1] is False          # 末尾の 0 まで書いた値
    assert _coordinate("0.3330(2)")[1] is False       # 不確かさ付き = 精密化された値
    assert _coordinate("0.334")[1] is False
    assert _coordinate("0.25")[1] is False            # 小数で正確
