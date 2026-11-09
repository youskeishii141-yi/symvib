"""CIF の読み込み (実在の CIF の「くせ」に強い自前の展開処理)。

pymatgen の CIF パーサは字句解析 (CifFile) だけに使い、原子の展開は自前で行う。
対応している主なケース:
  - 複数のデータブロックを含む CIF        → 原子座標を持つ最初のブロックを使う
  - 不確かさ付きの数値 "4.583(2)"          → 括弧を外して読む
  - ラベルが元素記号でない "Ti1", "O2a"    → _atom_site_type_symbol、なければラベルから元素を推定
  - 酸化数付きの型 "O-2", "Si4+"           → 元素記号だけを取り出す
  - 部分占有 (同じ座標に複数のラベル)      → 1 つの仮想原子にまとめ、名前はラベルの連結 (例 "Ti1V2")
  - 実在しない元素記号 ("X", "Wat" など)   → 質量 1 amu の仮想原子として扱い、注記を残す
  - 対称操作の記載がない                    → 空間群名・番号から操作を補う (なければ P1)
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

import numpy as np

from .common import AnalysisError, register_pseudo_species
from .i18n import Msg

# 同じ位置とみなす距離 [Å] (対称操作で重なる像の統合・部分占有の検出)
DUPLICATE_TOL = 0.02
# 同位体など、元素記号として扱いたい別名 → (元素, 質量)
_ISOTOPES = {"D": ("H", 2.014), "T": ("H", 3.016)}

_SYMOP_KEYS = ("_symmetry_equiv_pos_as_xyz", "_space_group_symop_operation_xyz",
               "_space_group_symop.operation_xyz", "_symmetry_equiv.pos_as_xyz")
_HM_KEYS = ("_symmetry_space_group_name_h-m", "_space_group_name_h-m_alt",
            "_space_group.name_h-m_alt", "_symmetry.space_group_name_h-m")
_ITN_KEYS = ("_space_group_it_number", "_symmetry_int_tables_number",
             "_space_group.it_number", "_symmetry.int_tables_number")


@dataclass
class CifStructure:
    lattice: np.ndarray            # 格子ベクトル (行ベクトル, Å)
    frac: np.ndarray               # 分率座標 (N, 3), [0, 1) に収めたもの
    numbers: list[int]             # 原子の種類 (原子番号; 仮想原子は 1000 以上)
    names: list[str]               # 表示名 (元素記号、部分占有ならラベルの連結)
    labels: list[str]              # CIF のサイトラベル (部分占有なら連結)
    block: str = ""                # 使ったデータブロック名
    notes: list[str] = field(default_factory=list)


def _number(value):
    """CIF の数値 ("0.4006(3)", "?", ".") を float に (読めなければ None)。"""
    if value is None:
        return None
    s = str(value).strip()
    if s in ("?", ".", ""):
        return None
    s = re.sub(r"\(.*?\)", "", s)
    try:
        if "/" in s:                                  # "1/2" のような分数
            num, den = s.split("/", 1)
            return float(num) / float(den)
        return float(s)
    except (ValueError, ZeroDivisionError):
        return None


# 小数で丸めて書かれやすい分数 (1/3, 1/6, 1/12 の倍数)。1/2, 1/4, 1/8 は小数で正確に書ける
_FRACTION_DENOMINATORS = (3, 6, 12)


def _coordinate(value):
    """原子の分率座標を読む。返り値 (値 or None, 分数に直したか)。

    "1/3" のような分数はそのまま正確に読む。"0.333", "0.6667" のように 1/3, 2/3, 1/6 …
    を小数で丸めた (または切り捨てた) 書き方は、その分数そのものとみなす。
    古い CIF に多い 3 桁の "0.333" のままだと、許容差 (既定 0.001 Å) を超えて対称性が
    下がってしまうため (例: 黒鉛が P6_3/mmc でなく Cmcm になる)。
    不確かさ付き ("0.3330(2)") は精密化された自由パラメータなので、そのまま使う。"""
    v = _number(value)
    s = str(value).strip() if value is not None else ""
    if v is None or "(" in s or "/" in s:
        return v, False
    m = re.fullmatch(r"[+-]?\d*\.(\d+)", s)
    if not m or len(m.group(1)) < 2:
        return v, False
    digits = len(m.group(1))
    scale = 10 ** digits
    # 小数 2 桁 ("0.33") は粗いので 1/3 の倍数だけ。3 桁以上なら 1/6, 1/12 の倍数も
    for n in (_FRACTION_DENOMINATORS if digits >= 3 else (3,)):
        k = round(v * n)
        exact = k / n
        if abs(v - exact) < 1e-9:
            return v, False                     # 小数で正確に書けている (0.25, 0.5 など)
        if abs(v - exact) > 1.0 / scale:
            continue
        a = abs(exact) * scale
        written = abs(v) * scale
        # 丸め (四捨五入) か切り捨ての結果とちょうど一致するときだけ
        if round(written) in (round(a), int(a)) and abs(written - round(written)) < 1e-6:
            return exact, True
    return v, False


def _snap_translation(t, tol=1e-3):
    """対称操作の並進 (小数で書かれた 0.3333 など) を最も近い k/12 にそろえる。"""
    t = np.asarray(t, float)
    snapped = np.round(t * 12) / 12
    return np.where(np.abs(t - snapped) < tol, snapped, t)


def _element_from(token):
    """"O-2", "Si4+", "Ti1", "SI1", "Fe2a" などから元素記号と質量を推定する。

    返り値 (元素記号 or None, 質量 or None)。実在元素でなければ (None, None)。"""
    from pymatgen.core import Element
    m = re.match(r"\s*([A-Za-z]+)", str(token))
    if not m:
        return None, None
    letters = m.group(1)
    # "Wat", "Uiso" のように小文字が 2 文字以上続く語は元素記号ではない
    if len(letters) >= 3 and letters[1:].islower():
        return None, None
    cands = []
    # 2 文字の元素は 2 文字目が小文字 ("Si1") か、全部大文字 ("SI1") のときだけ
    if len(letters) >= 2 and (letters[1].islower() or letters[:2].isupper()):
        cands.append(letters[:2].capitalize())
    cands.append(letters[:1].upper())
    for cand in cands:
        if cand in _ISOTOPES:
            return _ISOTOPES[cand]
        try:
            el = Element(cand)
        except ValueError:
            continue
        if el.Z <= 103:                               # 超重元素 (Nh など) はラベルの誤読
            return el.symbol, float(el.atomic_mass)
    return None, None


def _ci(block_data):
    """キーを小文字にした辞書 (CIF のタグは大文字小文字を区別しない)。"""
    return {k.lower(): v for k, v in block_data.items()}


def _first(d, keys):
    for k in keys:
        if k in d:
            v = d[k]
            return v[0] if isinstance(v, list) else v
    return None


_HALL_KEYS = ("_symmetry_space_group_name_hall", "_space_group_name_hall",
              "_space_group.name_hall")


def _norm_sg(name):
    """空間群記号の比較用の正規化 (空白・下線・引用符を除き小文字に)。"""
    return re.sub(r"[\s_'\"]", "", str(name)).lower()


def _setting_names(t):
    """spglib の空間群型 t (1 つの設定) の呼び名の集合 (正規化済み)。"""
    names = {_norm_sg(t.international_short), _norm_sg(t.international_full),
             _norm_sg(t.international)}
    tokens = t.international_full.split()
    if 3 <= t.number <= 15:                           # 単斜晶: "P 1 2_1/n 1" → P2_1/n
        names.add(_norm_sg(tokens[0] + "".join(x for x in tokens[1:] if x != "1")))
    elif 16 <= t.number <= 74:                        # 直方晶: "P 2_1/b 2_1/n 2_1/m" → Pbnm
        names.add(_norm_sg(tokens[0] + "".join(x.split("/")[-1] for x in tokens[1:])))
    return names


def _symmetry_candidates(d, lattice_params):
    """対称操作の候補 [(説明, [(R, t)])] と、記載から決まったかどうか。

    対称操作の記載があればそれだけを使う。なければ Hall 記号 → H-M 記号 → 番号の順に
    spglib の空間群データベース (530 設定) から候補を集める。原点の選び方
    (origin 1/2) や菱面体の軸 (H/R) が記号から決まらなければ、すべての候補を返し、
    呼び出し側が展開結果で選ぶ (原子が重なりすぎない・数が少ないもの)。"""
    import spglib
    from pymatgen.core.operations import SymmOp
    for k in _SYMOP_KEYS:
        if k in d:
            vals = d[k] if isinstance(d[k], list) else [d[k]]
            ops = []
            for s in vals:
                op = SymmOp.from_xyz_str(str(s).replace(" ", "").lower())
                ops.append((op.rotation_matrix, _snap_translation(op.translation_vector)))
            return [(None, ops)]
    types = [spglib.get_spacegroup_type(h) for h in range(1, 531)]

    def ops_of(t):
        sym = spglib.get_symmetry_from_database(t.hall_number)
        return [(np.array(R, float), np.array(tr, float))
                for R, tr in zip(sym["rotations"], sym["translations"])]

    hall = _first(d, _HALL_KEYS)
    if hall not in (None, "?", "."):
        hits = [t for t in types if _norm_sg(t.hall_symbol) == _norm_sg(hall)]
        if hits:
            return [(f"{hits[0].international_full} (Hall {hits[0].hall_symbol})", ops_of(hits[0]))]
    cands = []
    hm = _first(d, _HM_KEYS)
    if hm not in (None, "?", "."):
        text = str(hm).strip()
        choice = None
        m = re.search(r":\s*(\w+)\s*$", text)            # "F d -3 m :2", "R -3 c :H"
        if m:
            choice, text = m.group(1).upper(), text[:m.start()]
        else:
            m = re.search(r"\s([ZSHR])\s*$", text)       # ICSD の "F d -3 m Z" (Z = origin 2)
            if m:
                choice, text = m.group(1), text[:m.start()]
        choice = {"Z": "2", "S": "1"}.get(choice, choice)
        cands = [t for t in types if _norm_sg(text) in _setting_names(t)]
        if choice:
            cands = [t for t in cands if t.choice.upper().endswith(choice)] or cands
    if not cands:
        itn = _number(_first(d, _ITN_KEYS))
        if itn is not None:
            cands = [t for t in types if t.number == int(itn)]
    a, b, c, al, be, ga = lattice_params
    rhombo = (abs(a - b) < 1e-3 and abs(b - c) < 1e-3 and abs(al - 90) > 0.1
              and abs(al - be) < 1e-3 and abs(be - ga) < 1e-3)
    cands = [t for t in cands if t.choice not in ("H", "R") or t.choice == ("R" if rhombo else "H")]
    # 単斜晶: 90° でない角が 1 つだけなら、その軸 (α → a, β → b, γ → c) を unique axis とする
    # 設定に絞る (記号 "P 21/c" だけでは b 軸か c 軸かが決まらない)
    oblique = [abs(x - 90) > 0.1 for x in (al, be, ga)]
    if sum(oblique) == 1 and any(3 <= t.number <= 15 for t in cands):
        axis = "abc"[oblique.index(True)]
        same_axis = [t for t in cands if t.number < 3 or t.number > 15 or axis in t.choice]
        cands = same_axis or cands
    return [(f"{t.international_full}" + (f" ({t.choice})" if t.choice else ""), ops_of(t))
            for t in cands]


def _find_block(path):
    from pymatgen.io.cif import CifFile
    try:
        with open(path, encoding="utf-8-sig", errors="replace") as f:   # BOM 付きも読む
            cif = CifFile.from_str(f.read())
    except Exception as e:
        raise AnalysisError(Msg("err.cif_read", path=path, err=e)) from e
    blocks = [(name, _ci(b.data)) for name, b in cif.data.items()]
    usable = [(n, d) for n, d in blocks if "_atom_site_fract_x" in d and "_cell_length_a" in d]
    if not usable:
        raise AnalysisError(Msg("err.cif_no_atoms"))
    notes = []
    if len(usable) > 1:
        notes.append(Msg("cif.multi", n=len(usable), name=usable[0][0]))
    return usable[0][0], usable[0][1], notes


def read_cif(path) -> CifStructure:
    """CIF を読み、対称操作で単位胞内の全原子に展開して返す。"""
    from pymatgen.core import Element, Lattice

    block, d, notes = _find_block(path)
    cell = [_number(_first(d, [k])) for k in ("_cell_length_a", "_cell_length_b",
                                              "_cell_length_c", "_cell_angle_alpha",
                                              "_cell_angle_beta", "_cell_angle_gamma")]
    if any(v is None for v in cell):
        raise AnalysisError(Msg("err.cif_cell"))
    lattice = Lattice.from_parameters(*cell).matrix

    def column(key, default=None):
        v = d.get(key)
        if v is None:
            return [default] * len(d["_atom_site_fract_x"])
        return v if isinstance(v, list) else [v]

    xs, ys, zs = (column(f"_atom_site_fract_{c}") for c in "xyz")
    labels = column("_atom_site_label")
    types = column("_atom_site_type_symbol")
    occs = column("_atom_site_occupancy", "1")

    # ── 非対称単位のサイト ──
    asym = []
    unknown = []
    snapped = []                                 # 分数に直した座標の値
    for i, (x, y, z) in enumerate(zip(xs, ys, zs)):
        parsed = [_coordinate(x), _coordinate(y), _coordinate(z)]
        pos = [p for p, _ in parsed]
        snapped += [f"{p:g}" for (p, s), raw in zip(parsed, (x, y, z)) if s]
        if any(p is None for p in pos):
            raise AnalysisError(Msg("err.cif_coords", label=labels[i] or i + 1))
        label = str(labels[i]) if labels[i] not in (None, "?", ".") else f"X{i + 1}"
        elem, mass = _element_from(types[i]) if types[i] not in (None, "?", ".") else (None, None)
        if elem is None:
            elem, mass = _element_from(label)
        if elem is None:
            unknown.append(label)
        occ = _number(occs[i])
        asym.append(dict(label=label, elem=elem, mass=mass,
                         occ=1.0 if occ is None else occ, frac=np.array(pos)))
    if snapped:
        notes.append(Msg("cif.fraction", n=len(snapped)))
    if unknown:
        notes.append(Msg("cif.unknown", sites=", ".join(unknown)))

    # ── 対称操作で展開し、重なる位置をまとめる ──
    def expand(ops):
        sites = []                               # dict(frac, members=[asym index])
        for k, a in enumerate(asym):
            for R, t in ops:
                f = (R @ a["frac"] + t) % 1.0
                f[np.isclose(f, 1.0, atol=1e-6)] = 0.0
                hit = None
                for s in sites:
                    df = f - s["frac"]
                    df -= np.round(df)
                    if np.linalg.norm(df @ lattice) < DUPLICATE_TOL:
                        hit = s
                        break
                if hit is None:
                    sites.append(dict(frac=f, members=[k]))
                elif k not in hit["members"]:
                    hit["members"].append(k)       # 別ラベルが同じ位置 → 部分占有
        return sites

    def min_distance(sites):
        """展開した原子どうしの最短距離 [Å] (周期像を含む)。"""
        f = np.array([s["frac"] for s in sites])
        df = f[:, None, :] - f[None, :, :]
        df -= np.round(df)
        dist = np.linalg.norm(df @ lattice, axis=2)
        dist[np.eye(len(f), dtype=bool)] = np.inf
        return float(dist.min()) if len(f) > 1 else np.inf

    candidates = _symmetry_candidates(d, cell)
    if not candidates:
        notes.append(Msg("cif.p1"))
        candidates = [(None, [(np.eye(3), np.zeros(3))])]
    # 記号から設定 (原点・軸) が決まらないときは、原子が近づきすぎない (0.5 Å 以上)・
    # 原子数が最も少ない展開を採る (誤った原点では原子が重なったり増えたりする)
    expanded = [(name, expand(ops)) for name, ops in candidates]
    name, sites = min(expanded, key=lambda e: (min_distance(e[1]) < 0.5, len(e[1])))
    if name is not None:
        notes.append(Msg("cif.symops_from_sg", sg=name))

    numbers, names, site_labels = [], [], []
    partial, merged = set(), set()
    for s in sites:
        members = [asym[k] for k in s["members"]]
        # 同じ元素の別ラベルが同じ位置 (対称操作で重なる重複記載など) は 1 つの原子
        # (占有率は合計し、1 を上限とする)
        if len(members) > 1 and len({m["elem"] for m in members}) == 1                 and members[0]["elem"] is not None:
            top = max(members, key=lambda m: m["occ"])
            members = [dict(top, occ=min(sum(m["occ"] for m in members), 1.0))]
        label = "".join(m["label"] for m in members)
        site_labels.append(label)
        if len(members) > 1:
            merged.add(label)
        if len(members) == 1 and members[0]["elem"] is not None:
            m = members[0]
            if m["occ"] < 0.999:
                partial.add(f"{m['label']} ({m['occ']:g})")
            if m["elem"] == "H" and m["mass"] > 1.5:   # 重水素など
                z = register_pseudo_species(m["label"][:1], m["mass"], 1,
                                            key=("iso", m["mass"]))
                numbers.append(z)
                names.append(m["label"][:1])
            else:
                numbers.append(Element(m["elem"]).Z)
                names.append(m["elem"])
            continue
        # 部分占有 (複数ラベル) または未知の元素 → 仮想原子
        weights = [m["occ"] for m in members]
        masses = [m["mass"] if m["mass"] is not None else 1.0 for m in members]
        mass = float(np.dot(weights, masses) / max(sum(weights), 1e-12))
        main = max(members, key=lambda m: (m["elem"] is not None, m["occ"]))
        zdisp = Element(main["elem"]).Z if main["elem"] else 0
        key = (label, tuple(sorted((m["elem"] or m["label"], round(m["occ"], 3))
                                   for m in members)))
        numbers.append(register_pseudo_species(label, mass, zdisp, key=key))
        names.append(label)
    if merged:
        notes.append(Msg("cif.merged", sites=", ".join(sorted(merged))))
    if partial:
        notes.append(Msg("cif.partial", sites=", ".join(sorted(partial))))

    frac = np.array([s["frac"] for s in sites])
    return CifStructure(lattice=lattice, frac=frac, numbers=numbers, names=names,
                        labels=site_labels, block=block, notes=notes)

