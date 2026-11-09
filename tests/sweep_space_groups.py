"""全 230 空間群・530 Hall 設定の合成構造で結晶解析を網羅的に検査する (数分かかる)。

通常の pytest からは外してあるので、直接実行する:
    python tests/sweep_space_groups.py [--hall 1-530] [--md docs/validation_space_groups.ja.md]

各 Hall 設定について
  1. spglib の対称操作から、その設定で不変な格子 (計量を群で平均) と、一般位置 + 特殊位置の
     合成構造 (軌道ごとに別の元素) を作る
  2. analyze_crystal で解析し、tests/invariants.check_crystal の独立検査 (既約性・活性・
     テンソル成分・相関法の重複度・Bilbao 指標表のラベル) を通す
  3. 同じ構造を spglib の標準設定 (idealized) に直して解析し、既約分解 (ラベルの多重集合) が
     一致することを確かめる
結果は空間群ごとに集計し、--md でマークダウンの表に書き出す。
"""
from __future__ import annotations

import argparse
import os
import sys
import time
import traceback
from collections import Counter

import numpy as np
import spglib

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from invariants import check_crystal                      # noqa: E402
from symvib.crystal import analyze_crystal     # noqa: E402

# 特殊位置の候補 (この点を含む軌道を作る; 設定によっては一般位置になることもある)
SPECIAL_POINTS = [
    (0, 0, 0), (0.5, 0.5, 0.5), (0, 0, 0.5), (0.5, 0, 0), (0, 0.5, 0), (0.5, 0.5, 0),
    (0.25, 0.25, 0.25), (0, 0, 0.3137), (0.2731, 0, 0), (0.3137, 0.3137, 0),
    (0.3137, 0.3137, 0.3137), (1 / 3, 2 / 3, 0), (1 / 3, 2 / 3, 0.25), (1 / 3, 2 / 3, 0.1713),
    (0.1713, 0.25, 0.4123), (0.25, 0.1713, 0.4123), (0, 0.1713, 0.25), (0.1713, 0, 0.25),
    (0.1713, 0.4123, 0.25), (0.1713, 0.1713, 0.4123), (0.1713, 0.4123, 0.4123),
    (0.1713, 0, 0.4123), (0, 0.1713, 0.4123), (0.1713, 0.4123, 0),
]
# 一般位置は 2 軌道 (別元素) 置く: 1 軌道だけだと偶然の対称性 (P2_1 の 2 原子は P2_1/m) が出る
GENERAL_POINTS = [(0.1234, 0.2171, 0.3417), (0.4123, 0.0871, 0.2311)]


def invariant_lattice(rots, seed=0):
    """全回転 W で W^T G W = G となる計量 G の格子 (行ベクトル, Å)。三斜なら一般の格子。"""
    rng = np.random.default_rng(seed)
    g0 = np.diag([5.0, 6.1, 7.3]) ** 2 + 0.9 * rng.standard_normal((3, 3))
    g0 = 0.5 * (g0 + g0.T) + 15 * np.eye(3)
    g = sum(np.asarray(w, float).T @ g0 @ np.asarray(w, float) for w in rots) / len(rots)
    return np.linalg.cholesky(0.5 * (g + g.T))              # L L^T = G → 行ベクトルが基底


def orbit(point, rots, trans, lattice, tol=1e-3):
    pts = []
    for w, t in zip(rots, trans):
        f = (np.asarray(w) @ np.asarray(point, float) + t) % 1.0
        f[f > 1 - 1e-9] = 0.0
        if not any(np.linalg.norm((((f - p) + 0.5) % 1.0 - 0.5) @ lattice) < tol for p in pts):
            pts.append(f)
    return np.array(pts)


def synthetic_structure(hall, n_special=3, min_dist=1.2):
    """Hall 番号の設定で、一般位置 1 軌道 + 特殊位置 n_special 軌道の構造を作る。"""
    sym = spglib.get_symmetry_from_database(hall)
    rots, trans = sym["rotations"], sym["translations"]
    lattice = invariant_lattice(rots, seed=hall)
    frac, numbers = [], []
    def add(orb, z):
        nonlocal frac, numbers
        if len(frac):
            d = orb[:, None, :] - np.asarray(frac)[None, :, :]
            d -= np.round(d)
            if np.min(np.linalg.norm(d @ lattice, axis=2)) < min_dist:
                return False
        frac += list(orb)
        numbers += [z] * len(orb)
        return True
    for z, p in enumerate(GENERAL_POINTS, 1):
        add(orbit(p, rots, trans, lattice), z)
    n_general = len(orbit(GENERAL_POINTS[0], rots, trans, lattice))
    z = len(GENERAL_POINTS) + 1
    for p in SPECIAL_POINTS:
        if z > n_special + len(GENERAL_POINTS):
            break
        orb = orbit(p, rots, trans, lattice)
        if len(orb) < n_general and add(orb, z):                  # 一般位置より小さい軌道だけ
            z += 1
    return lattice, np.array(frac), numbers


def decomposition(result):
    return dict(Counter(m.label for m in result.modes))


def phonopy_decomposition(result, seed=0):
    """phonopy (別の実装; 指標表は Bilbao 由来) で同じ primitive cell の既約分解を求める。

    対称化した乱数の力定数で動力学行列を作り、固有ベクトルの縮退組ごとに phonopy が付ける
    ラベルを数える。phonopy が入っていなければ None。"""
    try:
        from phonopy import Phonopy
        from phonopy.structure.atoms import PhonopyAtoms
    except ImportError:
        return None
    from invariants import mech_rep_from_ops
    lattice = np.asarray(result.lattice, float)
    frac = np.linalg.solve(lattice.T, np.asarray(result.cart).T).T
    numbers = [int(z) for z in result.numbers]
    sym = spglib.get_symmetry((lattice, frac, numbers), symprec=1e-3)
    D, _, _ = mech_rep_from_ops(lattice, frac, numbers, sym["rotations"], sym["translations"])
    na = len(numbers)
    rng = np.random.default_rng(seed)
    phi0 = rng.standard_normal((3 * na, 3 * na))
    phi0 = phi0 + phi0.T
    phi = np.mean([d.T @ phi0 @ d for d in D], axis=0) + 30.0 * np.eye(3 * na)
    # 仮想原子 (番号 ≥ 1000) は phonopy に渡せないので、種類を保って実在元素に置き換える
    kinds = sorted(set(numbers))
    z_std = [kinds.index(z) + 1 for z in numbers]
    cell = PhonopyAtoms(numbers=z_std, cell=lattice, scaled_positions=frac)
    ph = Phonopy(cell, supercell_matrix=np.eye(3, dtype=int), primitive_matrix=np.eye(3),
                 log_level=0)
    ph.force_constants = phi.reshape(na, 3, na, 3).transpose(0, 2, 1, 3)
    try:
        ph.set_irreps([0, 0, 0], degeneracy_tolerance=1e-3)
    except RuntimeError as e:                 # phonopy が扱えない向き (P312 系など "Rotation symbols not found")
        return f"unlabelled ({e})"
    labels = ph.irreps._ir_labels
    if labels is None or any(lbl is None for lbl in labels):
        return "unlabelled"
    return dict(Counter(labels))


def run_hall(hall, bilbao=True, phonopy=False):
    """1 つの Hall 設定を検査して結果 dict を返す。"""
    t0 = time.time()
    t = spglib.get_spacegroup_type(hall)
    out = dict(hall=hall, number=t.number, setting=t.international_full, choice=t.choice,
               ok=True, problems=[], notes=[])
    try:
        lat, frac, num = synthetic_structure(hall)
        out["natoms"] = len(num)
        r = analyze_crystal(lat, frac, num)
        out["decomposition"] = decomposition(r)
        out["sg_detected"] = r.space_group
        if not r.space_group.endswith(f"(No. {t.number})"):
            out["problems"].append(f"detected {r.space_group}, expected No. {t.number}")
        probs, info = check_crystal(r, lat, frac, num, bilbao=bilbao)
        out["notes"] += info.get("bilbao_notes", [])
        # 設定の向きが Bilbao の点群表と違う群 (-4m2 など) では、B1/B2 の成分の違いは慣例の差として
        # 別に記録する (Mulliken の規則: 添字 1 = 2 回軸に対称)
        if any("not along the Bilbao" in n for n in out["notes"]):
            out["convention"] = [p for p in probs if "Bilbao label" in p or "vs Bilbao" in p]
            probs = [p for p in probs if p not in out["convention"]]
        out["problems"] += probs
        out["bilbao_labels"] = info.get("bilbao_labels", {})
        out["wyckoffs"] = sorted(set(r.meta.get("wyckoffs", [])))
        # 標準設定 (idealized) に直した構造でも同じ分解になるか
        ds = spglib.get_symmetry_dataset((lat, frac, num), symprec=1e-3)
        r2 = analyze_crystal(ds.std_lattice, ds.std_positions, ds.std_types)
        if decomposition(r2) != out["decomposition"]:
            out["problems"].append(f"standard-setting result differs: {decomposition(r2)}")
        out["std_setting"] = bool(np.allclose(ds.transformation_matrix, np.eye(3), atol=1e-6))
        if phonopy:
            # phonopy は独立検査ではなく別実装との照合なので、不一致は合否とは別に記録する
            # (既知の例: P6_3/mcm で phonopy は B1/B2 の 2 回軸の向きを逆にとる)
            out["phonopy"] = phonopy_decomposition(r)
    except Exception as e:                                   # 失敗も記録して続ける
        out["problems"].append(f"EXCEPTION {type(e).__name__}: {e}")
        out["trace"] = traceback.format_exc()
    out["ok"] = not out["problems"]
    out["seconds"] = time.time() - t0
    return out


def fmt_decomp(d):
    return " + ".join(f"{n if n > 1 else ''}{k}" for k, n in sorted(d.items()))


def write_markdown(results, path):
    by_sg = {}
    for r in results:
        by_sg.setdefault(r["number"], []).append(r)
    L = ["# 230 空間群の網羅検査 (合成構造)", "",
         "`tests/sweep_space_groups.py` の結果。各空間群について spglib の全 Hall 設定で合成構造"
         " (一般位置 + 特殊位置) を解析し、独立検査 (`tests/invariants.py`: 既約性・活性・"
         "ラマンテンソル成分・相関法の重複度・Bilbao 点群表のラベル) と標準設定との一致を確かめた。", "",
         "| No. | 設定 (Hall 設定数) | 既約分解 (標準設定・合成構造) | 結果 |", "|---|---|---|---|"]
    n_ok = 0
    for sg in sorted(by_sg):
        rs = by_sg[sg]
        bad = [r for r in rs if not r["ok"]]
        d = fmt_decomp(rs[0].get("decomposition", {}))
        status = "OK" if not bad else "; ".join(
            f"{r['setting']}: " + "; ".join(r["problems"][:3]) for r in bad)
        conv = [r for r in rs if r.get("convention")]
        if conv:
            status += " (Bilbao 点群表と向きが異なる設定: " + "; ".join(sorted(set(
                p for r in conv for p in r["convention"]))) + ")"
        ph = rs[0].get("phonopy")
        if ph is not None:
            status += "; phonopy " + ("一致" if isinstance(ph, dict) and ph == rs[0].get("decomposition")
                                    else "照合不可" if not isinstance(ph, dict)
                                    else f"不一致 ({fmt_decomp(ph)})")
        n_ok += not bad
        L.append(f"| {sg} | {rs[0]['setting']} ({len(rs)}) | {d} | {status} |")
    L.insert(4, f"合格: {n_ok} / {len(by_sg)} 空間群, {sum(r['ok'] for r in results)} / {len(results)} 設定")
    L.insert(5, "")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--hall", default="1-530", help="範囲 (例 1-530, 90, 290-300)")
    ap.add_argument("--md", default=None, help="マークダウンの出力先")
    ap.add_argument("--no-bilbao", action="store_true")
    ap.add_argument("--phonopy", action="store_true", help="phonopy の irreps とも照合する")
    ap.add_argument("--json", default=None)
    ap.add_argument("--from-json", default=None, help="保存した結果からマークダウンだけ作り直す")
    a = ap.parse_args(argv)
    if a.from_json:
        import json
        with open(a.from_json, encoding="utf-8") as f:
            results = json.load(f)
        for r in results:                                     # 旧形式: phonopy の不一致が problems に入っている
            r["problems"] = [p for p in r["problems"] if not p.startswith("phonopy decomposition")]
            r["ok"] = not r["problems"]
        write_markdown(results, a.md)
        return 0
    lo, _, hi = a.hall.partition("-")
    halls = range(int(lo), int(hi or lo) + 1)
    results = []
    for h in halls:
        r = run_hall(h, bilbao=not a.no_bilbao, phonopy=a.phonopy)
        results.append(r)
        flag = "ok " if r["ok"] else "NG "
        print(f"{flag}{h:3d} No.{r['number']:3d} {r['setting']:24s} "
              f"{fmt_decomp(r.get('decomposition', {}))}"
              + ("" if r["ok"] else "\n      " + "\n      ".join(r["problems"][:6])), flush=True)
    if a.md:
        write_markdown(results, a.md)
    if a.json:
        import json
        with open(a.json, "w", encoding="utf-8") as f:
            json.dump(results, f, ensure_ascii=False, indent=1, default=str)
    bad = [r for r in results if not r["ok"]]
    print(f"\n{len(results) - len(bad)} / {len(results)} settings passed")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
