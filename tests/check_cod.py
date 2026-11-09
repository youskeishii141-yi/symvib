"""COD から取得した実在の CIF (tests/data/cod/) を全部解析して検査する (数分かかる)。

    python tests/fetch_cod.py          (CIF を取得; 配布物には含めない)
    python tests/check_cod.py [--md docs/validation_cod.ja.md]

各 CIF について
  1. cif_reader で読み、analyze_crystal で解析する (読み込み時の注記も記録)
  2. 検出した空間群が COD の空間群番号と一致するか (違えば symprec の問題か CIF の問題かを記録)
  3. tests/invariants.check_crystal の独立検査 (既約性・活性・テンソル成分・相関法・Bilbao のラベル)
  4. 同じ空間群・同じ Wyckoff 位置の合成構造 (元素を入れ替え、自由パラメータをずらし、格子定数を
     変えたもの) を解析し、既約分解が一致するか
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import traceback
from collections import Counter

import numpy as np
import spglib

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from invariants import check_crystal                             # noqa: E402
from sweep_space_groups import fmt_decomp                        # noqa: E402
from symvib.cif_reader import read_cif                # noqa: E402
from symvib.crystal import analyze_crystal            # noqa: E402

CODDIR = os.path.join(os.path.dirname(__file__), "data", "cod")


def same_wyckoff_synthetic(lattice, frac, numbers, symprec=1e-3, seed=0):
    """同じ空間群・同じ Wyckoff 位置をもつ別の構造 (自由パラメータ・元素・格子定数を変える)。

    各軌道の代表原子を乱数でずらし、そのサイト対称群 (原子を不変にする操作) で平均して
    Wyckoff 位置の自由度の中に戻す。軌道ごとに別の元素を割り当てる。"""
    ds = spglib.get_symmetry_dataset((lattice, frac, numbers), symprec=symprec)
    rots, trans = ds.rotations, ds.translations
    lattice = np.asarray(lattice, float)
    rng = np.random.default_rng(seed)
    # 格子はそのまま使う (P222 や Pmmm のように記号が軸の入れ替えで変わらない群では、標準設定の
    # 軸の名前 (x, y, z) が格子定数の順序で決まるため、格子を変えると B1/B2/B3 が入れ替わる)
    new_lat = lattice
    new_frac, new_num = [], []
    for k, rep in enumerate(dict.fromkeys(ds.equivalent_atoms)):
        x0 = np.asarray(frac[rep], float)
        site = [(R, t) for R, t in zip(rots, trans)
                if np.linalg.norm((((R @ x0 + t - x0) + 0.5) % 1.0 - 0.5) @ lattice) < 2 * symprec]
        x = x0 + rng.uniform(-0.05, 0.05, 3)
        for _ in range(3):                                    # サイト対称群で平均 → 不変部分へ
            imgs = [R @ x + t for R, t in site]
            imgs = [x0 + ((i - x0 + 0.5) % 1.0 - 0.5) for i in imgs]
            x = np.mean(imgs, axis=0)
        orb = []
        for R, t in zip(rots, trans):
            f = (R @ x + t) % 1.0
            if not any(np.linalg.norm((((f - p) + 0.5) % 1.0 - 0.5) @ new_lat) < 1e-3 for p in orb):
                orb.append(f)
        new_frac += orb
        new_num += [k + 1] * len(orb)
    return new_lat, np.array(new_frac), new_num


def check_file(path, symprec=1e-3):
    out = dict(file=os.path.basename(path), problems=[], notes=[])
    try:
        c = read_cif(path)
        out["notes"] = [str(n) for n in c.notes]
        out["natoms"] = len(c.numbers)
        out["pseudo"] = any(z >= 1000 for z in c.numbers)
        r = analyze_crystal(c.lattice, c.frac, c.numbers, symprec=symprec)
        out["space_group"] = r.space_group
        out["decomposition"] = dict(Counter(m.label for m in r.modes))
        out["wyckoffs"] = sorted(Counter(r.meta.get("wyckoffs", [])).items())
        probs, info = check_crystal(r, c.lattice, c.frac, c.numbers, symprec=symprec)
        out["bilbao_notes"] = info.get("bilbao_notes", [])
        if any("not along the Bilbao" in n for n in out["bilbao_notes"]):
            probs = [p for p in probs if "Bilbao label" not in p and "vs Bilbao" not in p]
        out["problems"] += probs
        lat2, frac2, num2 = same_wyckoff_synthetic(c.lattice, c.frac, c.numbers, symprec)
        r2 = analyze_crystal(lat2, frac2, num2, symprec=symprec)
        out["synthetic_sg"] = r2.space_group
        d2 = dict(Counter(m.label for m in r2.modes))
        if d2 != out["decomposition"] or r2.space_group != r.space_group:
            out["problems"].append(f"synthetic same-Wyckoff structure differs: {r2.space_group} "
                                   f"{fmt_decomp(d2)}")
    except Exception as e:                                        # noqa: BLE001
        out["problems"].append(f"EXCEPTION {type(e).__name__}: {e}")
        out["trace"] = traceback.format_exc()
    out["ok"] = not out["problems"]
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--md", default=None)
    ap.add_argument("--json", default=None)
    ap.add_argument("--only", default=None, help="COD id (カンマ区切り)")
    a = ap.parse_args(argv)
    with open(os.path.join(CODDIR, "index.json"), encoding="utf-8") as f:
        index = json.load(f)["entries"]
    missing = [e["id"] for e in index.values()
               if e and not os.path.exists(os.path.join(CODDIR, f"{e['id']}.cif"))]
    if missing:
        print(f"{len(missing)} CIF files are missing (not distributed with the repository). "
              "Run  python tests/fetch_cod.py  first.")
        return 1
    results = []
    for sg in sorted(index, key=int):
        e = index[sg]
        if e is None:
            results.append(dict(sg_number=int(sg), file=None, ok=None, problems=["not available"]))
            print(f"sg {sg:>3}: not available in COD (ordered, ≤ 4 elements, small)")
            continue
        if a.only and e["id"] not in a.only.split(","):
            continue
        out = check_file(os.path.join(CODDIR, f"{e['id']}.cif"))
        out.update(sg_number=int(sg), cod=e)
        results.append(out)
        detected = out.get("space_group", "?")
        mismatch = "" if detected.endswith(f"(No. {sg})") else f"  !! COD says No. {sg}"
        print(f"{'ok' if out['ok'] else 'NG'} sg {sg:>3} COD {e['id']} {e['formula']:14s} "
              f"{detected:22s} {fmt_decomp(out.get('decomposition', {}))}{mismatch}")
        for p in out["problems"]:
            print("       ", p)
    if a.json:
        with open(a.json, "w", encoding="utf-8") as f:
            json.dump(results, f, ensure_ascii=False, indent=1, default=str)
    if a.md:
        write_markdown(results, a.md)
    bad = [r for r in results if r["ok"] is False]
    print(f"\n{sum(1 for r in results if r['ok'])} ok, {len(bad)} with problems, "
          f"{sum(1 for r in results if r['ok'] is None)} not available")
    return 1 if bad else 0


def write_markdown(results, path):
    L = ["# COD の実在 CIF による検査 (各空間群 1 件)", "",
         "`tests/check_cod.py` の結果。出典: Crystallography Open Database "
         "(https://www.crystallography.net/cod/, パブリックドメイン)。CIF は `tests/data/cod/` に保存。", "",
         "| No. | COD ID | 物質 | 検出した空間群 | 原子数 | 既約分解 | 読み込み時の注記 | 結果 |",
         "|---|---|---|---|---|---|---|---|"]
    for r in results:
        if r.get("file") is None:
            L.append(f"| {r['sg_number']} | — | (入手不可) | | | | | |")
            continue
        e = r["cod"]
        name = e.get("mineral") or e.get("chemname") or ""
        notes = "; ".join(r.get("notes", []))
        status = "OK" if r["ok"] else "; ".join(r["problems"][:3])
        L.append(f"| {r['sg_number']} | [{e['id']}]({e['url']}) | {e['formula']} {name} | "
                 f"{r.get('space_group', '?')} | {r.get('natoms', '')} | "
                 f"{fmt_decomp(r.get('decomposition', {}))} | {notes} | {status} |")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")


if __name__ == "__main__":
    sys.exit(main())
