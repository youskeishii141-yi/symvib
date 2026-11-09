"""COD の実在 CIF (小さいもの) の回帰テスト。CIF は配布物に含めず、tests/fetch_cod.py で取得する。

全 227 件は tests/check_cod.py で (数分)。ここでは単位胞の原子数が少ないものだけを通常のテストとして
回し、読み込み → 解析 → 独立検査 (invariants) が通ること、検出した空間群が COD の空間群番号と
一致する (または擬対称で高くなる既知の例である) ことを確かめる。
"""
import json
import os

import pytest

from invariants import check_crystal
from symvib.cif_reader import read_cif
from symvib.crystal import analyze_crystal

CODDIR = os.path.join(os.path.dirname(__file__), "data", "cod")
with open(os.path.join(CODDIR, "index.json"), encoding="utf-8") as _f:
    _INDEX = {int(k): v for k, v in json.load(_f)["entries"].items() if v}
SMALL = sorted(sg for sg, e in _INDEX.items() if e["atoms_in_cell"] <= 12)

# CIF は配布物に含めない (tests/fetch_cod.py で index.json の ID と同じものを取得する)


@pytest.mark.parametrize("sg", SMALL)
def test_small_cod_cif(sg):
    e = _INDEX[sg]
    path = os.path.join(CODDIR, f"{e['id']}.cif")
    if not os.path.exists(path):
        pytest.skip(f"COD {e['id']} の CIF が未取得 (python tests/fetch_cod.py で取得)")
    c = read_cif(path)
    r = analyze_crystal(c.lattice, c.frac, c.numbers)
    detected = int(r.space_group.rsplit("No. ", 1)[1].rstrip(")"))
    # 宣言より高い対称性が見つかることはある (理想化された座標の擬対称); 低くなってはいけない
    assert _order(detected) >= _order(sg), (e["id"], r.space_group, f"COD: No. {sg}")
    probs, info = check_crystal(r, c.lattice, c.frac, c.numbers)
    if any("not along the Bilbao" in n for n in info.get("bilbao_notes", [])):
        probs = [p for p in probs if "Bilbao label" not in p and "vs Bilbao" not in p]
    assert probs == [], (e["id"], probs)


_PG_ORDER = {"1": 1, "-1": 2, "2": 2, "m": 2, "2/m": 4, "222": 4, "mm2": 4, "mmm": 8, "4": 4, "-4": 4,
             "4/m": 8, "422": 8, "4mm": 8, "-42m": 8, "4/mmm": 16, "3": 3, "-3": 6, "32": 6, "3m": 6,
             "-3m": 12, "6": 6, "-6": 6, "6/m": 12, "622": 12, "6mm": 12, "-6m2": 12, "6/mmm": 24,
             "23": 12, "m-3": 24, "432": 24, "-43m": 24, "m-3m": 48}


def _order(number):
    """空間群番号 → 点群の位数 (対称性の高さの比較用)。"""
    import spglib
    hall = next(h for h in range(1, 531) if spglib.get_spacegroup_type(h).number == number)
    return _PG_ORDER[spglib.get_spacegroup_type(hall).pointgroup_international.strip()]
