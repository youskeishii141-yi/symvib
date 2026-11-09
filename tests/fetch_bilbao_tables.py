"""Bilbao Crystallographic Server の点群指標表 (32 点群) を取得して JSON に保存する。

出典: Bilbao Crystallographic Server, "Point Group Tables"
      https://www.cryst.ehu.eus/cgi-bin/rep/programs/sam/point.py?sg=<空間群番号>
      (表の原典: Bradley & Cracknell (1972); Altmann & Herzig (1994))

保存先: tests/data/bilbao/point_tables.json
    { "<点群記号 (Hermann-Mauguin)>": {
          "schoenflies": "C4h", "source_sg": 83, "url": "...",
          "classes": [{"symbol": "4+", "mult": 1}, ...],
          "irreps":  [{"label": "Ag", "chars": [[re, im], ...], "functions": "x2+y2,z2,Jz"}, ...]
      }, ... }
複素共役ペア (1Eg / 2Eg) は別々の行として保存する (chars は複素数)。

使い方:  python tests/fetch_bilbao_tables.py
(ネットワークが必要。取得した表は配布物に含めない。未取得の場合、Bilbao との照合は省略される。)
"""
import cmath
import json
import os
import re
import sys
import urllib.request
from datetime import date
from html import unescape

URL = "https://www.cryst.ehu.eus/cgi-bin/rep/programs/sam/point.py?sg={sg}&num=1"
OUT = os.path.join(os.path.dirname(__file__), "data", "bilbao", "point_tables.json")

# 32 点群の代表空間群 (表は点群ごとなので、どの空間群を指定しても同じ表が返る)
REPRESENTATIVE_SG = {
    "1": 1, "-1": 2, "2": 3, "m": 6, "2/m": 10, "222": 16, "mm2": 25, "mmm": 47,
    "4": 75, "-4": 81, "4/m": 83, "422": 89, "4mm": 99, "-42m": 111, "4/mmm": 123,
    "3": 143, "-3": 147, "32": 149, "3m": 156, "-3m": 162,
    "6": 168, "-6": 174, "6/m": 175, "622": 177, "6mm": 183, "-6m2": 187, "6/mmm": 191,
    "23": 195, "m-3": 200, "432": 207, "-43m": 215, "m-3m": 221,
}

_W = cmath.exp(2j * cmath.pi / 3)


def _complex(token):
    """表の指標 ('1', '-1j', 'w', '-w2' など) を複素数に。"""
    s = token.strip().replace("^", "")
    table = {"w": _W, "w2": _W ** 2, "-w": -_W, "-w2": -_W ** 2}
    if s in table:
        return table[s]
    s = s.replace("i", "j")
    return complex(s)


def _cells(row_html):
    cells = re.findall(r"<td[^>]*>(.*?)</td>", row_html, flags=re.S)
    out = []
    for c in cells:
        c = c.replace("<br>", "\n")
        c = re.sub(r"<sup>(.*?)</sup>", r"^\1", c)
        c = re.sub(r"<sub>(.*?)</sub>", r"_\1", c)
        c = re.sub(r"<[^>]+>", "", c)
        out.append(unescape(c).strip())
    return out


def _label(raw):
    """'A_g' → 'Ag', "A'_1" → "A1'", '^1E_g' → '1Eg'。"""
    s = raw.replace("^", "").replace("_", "")
    # Bilbao はプライムを添字の前に書く (A'1) → 慣用の A1'
    m = re.fullmatch(r"([A-Z])('+)(\d?)([gu]?)", s)
    if m:
        s = m.group(1) + m.group(3) + m.group(2) + m.group(4)
    return s


def parse(html, sg):
    i = html.find("Character Table")
    j = html.find("</table>", i)
    seg = html[i:j]
    m = re.search(r"group\s*<nobr>(.*?)</nobr>", seg, flags=re.S)
    schoen = re.sub(r"<[^>]+>", "", m.group(1)) if m else ""
    schoen = unescape(schoen).replace("_", "")
    rows = [_cells(r) for r in re.findall(r"<tr>(.*?)</tr>", seg, flags=re.S)]
    header = rows[0]
    has_fn = header[-1] == "functions"           # C1 の表には functions 列がない
    end = -1 if has_fn else len(header)
    classes = [{"symbol": c.replace("^", "").replace("_", "")} for c in header[2:end]]
    body = [r for r in rows[1:] if len(r) >= 3]  # 'w = exp(2iπ/3)' のような注記行を除く
    if body and body[0][0] == "Mult.":
        for k, c in enumerate(body[0][2:end]):
            classes[k]["mult"] = int(c)
        body = body[1:]
    else:
        for c in classes:
            c["mult"] = 1
    irreps = []
    for r in body:
        labels = r[0].split("\n")
        n = len(labels)
        chars = [c.split("\n") for c in r[2:end]]
        fn = r[-1] if has_fn and r[-1] != "·" else ""
        for k in range(n):
            z = [_complex(col[k]) for col in chars]
            irreps.append({"label": _label(labels[k]), "gamma": r[1].split("\n")[k],
                           "chars": [[v.real, v.imag] for v in z], "functions": fn})
    return dict(schoenflies=re.sub(r"\(.*\)", "", schoen).strip(), source_sg=sg,
                url=URL.format(sg=sg), classes=classes, irreps=irreps)


def main():
    tables = {}
    for hm, sg in REPRESENTATIVE_SG.items():
        with urllib.request.urlopen(URL.format(sg=sg), timeout=60) as f:
            html = f.read().decode("latin-1")
        tables[hm] = parse(html, sg)
        print(hm, tables[hm]["schoenflies"], [c["symbol"] for c in tables[hm]["classes"]],
              [ir["label"] for ir in tables[hm]["irreps"]])
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(dict(source="Bilbao Crystallographic Server, Point Group Tables "
                              "(https://www.cryst.ehu.eus/rep/point.html)",
                       retrieved=str(date.today()), tables=tables), f, indent=1,
                  ensure_ascii=False)
    print("saved", OUT)


if __name__ == "__main__":
    sys.exit(main())
