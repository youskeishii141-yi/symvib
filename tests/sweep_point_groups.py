"""収録している 45 点群すべてについて、一般位置の軌道から作った分子で分子解析を検査する (数分)。

    .venv\\Scripts\\python.exe tests\\sweep_point_groups.py

各点群について
  1. 生成元から理想的な群を作り、一般位置 2 軌道 (+ 中心原子) の分子を作る (原子が近すぎる点は避ける)
  2. 解析して tests/invariants.check_molecule の独立検査 (既約性・直交性・自由度・活性・運動量・
     角運動量・対称操作での閉性) を通す
  3. 任意の回転・原子の並べ替え・0.02 Å の歪みを加えても、点群と既約分解が変わらず検査に通る
"""
from __future__ import annotations

import math
import os
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from invariants import check_molecule                         # noqa: E402
from structures import _closure, _mol, _orbit, distorted       # noqa: E402
from test_tables import generators                            # noqa: E402
from symvib import molecule as M                  # noqa: E402
from symvib.molecule_tables import CHAR_TABLES, COMPLEX_TABLES   # noqa: E402


def generic_molecule(pg, seed=0):
    """一般位置 2 軌道 (C, H) + 中心原子 N の分子。原子間距離 0.6 Å 以上で、点群が pg になるものを選ぶ
    (乱数の位置によっては偶然高い対称性 (正六角形など) になるので、何度か引き直す)。"""
    group = _closure(generators(pg))
    scale = max(1.0, math.sqrt(len(group) / 8))                # 大きな群ほど軌道を広げる
    rng = np.random.default_rng(seed)
    for _ in range(60):
        p1 = scale * (np.array([1.3, 0.4, 0.6]) + rng.uniform(-0.4, 0.4, 3))
        p2 = scale * (np.array([2.1, 1.3, 0.2]) + rng.uniform(-0.4, 0.4, 3))
        atoms = _orbit("C", p1, group) + _orbit("H", p2, group)
        if len(group) > 1:
            atoms = [("N", (0, 0, 0))] + atoms
        if len(group) <= 2:                                      # C1, Ci, Cs, C2: 2 軌道だけでは
            p3 = scale * (np.array([0.3, -1.1, 0.4]) + rng.uniform(-0.4, 0.4, 3))   # 平面・直線になる
            atoms += _orbit("O", p3, group)
        if len(group) == 1:                                      # C1: 4 原子目で平面も避ける
            atoms += [("F", tuple(scale * (np.array([-0.9, 0.5, -1.2]) + rng.uniform(-0.4, 0.4, 3))))]
        cart = np.array([a[1] for a in atoms], float)
        d = np.linalg.norm(cart[:, None] - cart[None], axis=2)
        d[np.eye(len(d), dtype=bool)] = 9
        if d.min() < 0.6:
            continue
        mol = _mol(atoms)
        try:
            if M.analyze_molecule(*mol).point_group == pg:
                return mol
        except Exception:                                        # noqa: BLE001 (ほぼ直線など) → 引き直す
            continue
    raise RuntimeError(pg)


def main():
    rng = np.random.default_rng(3)
    bad = 0
    for pg in sorted(set(CHAR_TABLES) | set(COMPLEX_TABLES)):
        cart, Z, sym = generic_molecule(pg)
        probs = []
        try:
            r = M.analyze_molecule(cart, Z, sym)
            if r.point_group != pg:
                probs.append(f"point group {r.point_group}")
            probs += check_molecule(r)
            base = Counter(m.label for m in r.modes)
            q, _ = np.linalg.qr(rng.standard_normal((3, 3)))
            perm = rng.permutation(len(Z))
            r2 = M.analyze_molecule(distorted((cart @ q.T)[perm], 0.02, seed=5),
                                    [Z[i] for i in perm], [sym[i] for i in perm])
            if r2.point_group != r.point_group or Counter(m.label for m in r2.modes) != base:
                probs.append(f"rotated/permuted/distorted: {r2.point_group} "
                             f"{dict(Counter(m.label for m in r2.modes))}")
            probs += [f"(distorted) {p}" for p in check_molecule(r2)]
        except Exception as e:                                  # noqa: BLE001
            probs.append(f"EXCEPTION {type(e).__name__}: {e}")
        bad += bool(probs)
        print(f"{'NG' if probs else 'ok'} {pg:4s} {len(Z):3d} atoms", flush=True)
        for p in probs:
            print("      ", p)
    print(f"\n{45 - bad} / 45 point groups passed")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
