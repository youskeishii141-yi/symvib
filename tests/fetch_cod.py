"""Crystallography Open Database (COD) から 230 空間群それぞれの代表 CIF を 1 つずつ取得する。

出典: COD https://www.crystallography.net/cod/ (パブリックドメイン)。
検索 API:  https://www.crystallography.net/cod/result?format=json&space_group_number=N&strictmax=3
CIF:       https://www.crystallography.net/cod/<id>.cif

選び方: 原子座標があり、元素数 ≤ 3、秩序化 (化学式に小数の係数がない)、
単位胞の原子数 (cellformula) が少ないものを優先する。取得した CIF は tests/data/cod/<id>.cif に
保存し、一覧を tests/data/cod/index.json に書く (space group → id, 物質名など)。

CIF は配布物に含めない (一覧 index.json だけを含める)。一覧にある空間群は検索し直さず、
記録した ID の CIF をそのまま取り直すので、検査 (tests/check_cod.py, tests/test_cod_regression.py)
を同じ入力で再現できる。

使い方:  python tests/fetch_cod.py [--sg 1-230] [--max-bytes 40000]
(ネットワークが必要。)
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
import urllib.request
from datetime import date

API = "https://www.crystallography.net/cod/result?format=json&space_group_number={sg}&strictmax={nel}"
CIF = "https://www.crystallography.net/cod/{id}.cif"
OUTDIR = os.path.join(os.path.dirname(__file__), "data", "cod")
INDEX = os.path.join(OUTDIR, "index.json")


def _get(url, timeout=120, tries=4):
    """URL を取得する。COD は短時間に続けて取りに行くと接続を拒否することがあるので、
    間隔を広げながら数回やり直す。"""
    req = urllib.request.Request(url, headers={"User-Agent": "symvib-tests"})
    for k in range(tries):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as f:
                return f.read()
        except OSError:
            if k == tries - 1:
                raise
            time.sleep(10 * 2 ** k)


def _atoms_in_cell(entry):
    """cellformula '- Ca4 O12 Ti4 -' → 原子数 20 (読めなければ大きな値)。"""
    s = entry.get("cellformula") or ""
    n, ok = 0, False
    for el, cnt in re.findall(r"([A-Z][a-z]?)(\d*\.?\d*)", s):
        ok = True
        if cnt and "." in cnt:
            return 10 ** 6                              # 部分占有 (小数の係数)
        n += int(cnt) if cnt else 1
    return n if ok else 10 ** 6


def _ordered(entry):
    """化学式 (formula / cellformula) に小数の係数がなければ秩序化した構造とみなす。"""
    return all("." not in (entry.get(k) or "") for k in ("formula", "cellformula"))


def candidates(sg, nel=3):
    data = json.loads(_get(API.format(sg=sg, nel=nel)).decode("utf-8"))
    rows = [e for e in data if "has coordinates" in (e.get("flags") or "") and _ordered(e)]
    rows.sort(key=lambda e: (_atoms_in_cell(e), int(e.get("nel") or 9), e["file"]))
    return rows


def fetch_one(sg, max_bytes, tried=12):
    """空間群 sg の代表 CIF を保存して index 用の情報を返す (見つからなければ None)。"""
    for nel in (2, 3, 4, 5, 6):
        rows = candidates(sg, nel)
        for e in rows[:tried]:
            cid = e["file"]
            try:
                text = _get(CIF.format(id=cid))
            except Exception as ex:                      # noqa: BLE001
                print(f"  {cid}: download failed ({ex})")
                continue
            if len(text) > max_bytes:
                continue
            path = os.path.join(OUTDIR, f"{cid}.cif")
            with open(path, "wb") as f:
                f.write(text)
            return dict(id=cid, url=CIF.format(id=cid), sg_number=sg, sg=e.get("sg"),
                        formula=(e.get("formula") or "").strip("- "), mineral=e.get("mineral"),
                        chemname=e.get("chemname") or e.get("commonname"),
                        atoms_in_cell=_atoms_in_cell(e), bytes=len(text))
        if rows:
            break                                       # 候補はあったが大きすぎた → 元素数を増やしても同じ
    return None


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--sg", default="1-230")
    ap.add_argument("--max-bytes", type=int, default=120000)
    ap.add_argument("--search-missing", action="store_true",
                    help="一覧で「入手不可」(null) の空間群を検索し直す (一覧が変わるので再現用には使わない)")
    a = ap.parse_args(argv)
    lo, _, hi = a.sg.partition("-")
    os.makedirs(OUTDIR, exist_ok=True)
    index = {}
    if os.path.exists(INDEX):
        with open(INDEX, encoding="utf-8") as f:
            index = json.load(f).get("entries", {})
    failed = []
    for sg in range(int(lo), int(hi or lo) + 1):
        if str(sg) in index and index[str(sg)] is None and not a.search_missing:
            continue                                     # 記録どおり「入手不可」のまま
        if str(sg) in index and index[str(sg)]:
            # 一覧にある ID はそのまま取り直す (検査の再現用。CIF は配布物に含めていない)
            e = index[str(sg)]
            path = os.path.join(OUTDIR, f"{e['id']}.cif")
            if not os.path.exists(path):
                try:
                    text = _get(e["url"])
                except Exception as ex:                  # noqa: BLE001
                    print(f"sg {sg}: download of {e['id']} failed ({ex})")
                    failed.append(e["id"])
                    continue
                with open(path, "wb") as f:
                    f.write(text)
                if len(text) != e.get("bytes", len(text)):
                    print(f"sg {sg}: {e['id']} differs in size from the recorded one "
                          f"({len(text)} vs {e['bytes']} bytes; COD may have revised the entry)")
                print(f"sg {sg:3d}: {e['id']} downloaded", flush=True)
                time.sleep(1.0)
            continue
        try:
            info = fetch_one(sg, a.max_bytes)
        except Exception as ex:                          # noqa: BLE001
            print(f"sg {sg}: query failed ({ex})")
            info = None
        index[str(sg)] = info
        print(f"sg {sg:3d}: " + (f"{info['id']} {info['formula']} ({info['atoms_in_cell']} atoms)"
                                if info else "not available"), flush=True)
        with open(INDEX, "w", encoding="utf-8") as f:
            json.dump(dict(source="Crystallography Open Database (https://www.crystallography.net/cod/)",
                           retrieved=str(date.today()), entries=index), f, indent=1, ensure_ascii=False)
        time.sleep(0.5)
    if failed:
        print(f"{len(failed)} downloads failed: {', '.join(failed)} (run again to retry)")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
