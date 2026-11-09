"""力ベース Γ 点ブロック凍結フォノン (Quantum ESPRESSO 連携)。

因子群解析の対称性適応変位基底を使い、pw.x の SCF を数回行うだけで、
各既約表現の Γ 点振動数と固有ベクトルを調和近似・有限変位の精度で求める。
極性結晶の方向依存の LO–TO 分裂 (非解析的補正) は含まない。

────────────────────────────────────────────────────────────────────────
物理 (対称性によるブロック分解)
────────────────────────────────────────────────────────────────────────
Γ 点 (q=0) では全 primitive cell が同位相で振動するのでスーパーセルは不要。
力定数 Φ は結晶の全対称操作と可換なので、各既約表現の部分空間を保つ
(Schur の補題)。したがって既約表現ごとに、その部分空間 (次元 m·d) に限定
して有限変位すれば、全行列と等価なブロック固有値問題が得られる。
力・有限差分・電子状態計算の近似誤差は別に残る。

  Φ b_α ≈ −[F(+δ b_α) − F(−δ b_α)] / (2δ)      (中心差分; 残留力を相殺)
  K_lα  = b_l · (Φ b_α)                          力定数ブロック [eV/Å²]
  M_lα  = Σ_i M_i (b_l,i · b_α,i)                質量計量      [amu]
  K v = ω² M v                                   一般化固有値問題
  ω[cm⁻¹] = 521.47 · √(固有値[eV/Å²/amu])

使い方は qe_frozen_phonon.py (CLI) を参照。
"""
from __future__ import annotations

import json
import os
import re

import numpy as np

from .common import (AnalysisError, atomic_mass, element_symbol, ensure_outdir,
                     generalized_eigh, safe_label, write_xsf)
from .crystal import graphite_2h, load_crystal, mech_rep, symmetry_adapted_blocks
from .i18n import Msg, pad, t

# ── 物理定数 (CODATA 2018) ──
RY_TO_EV = 13.605693122994
BOHR_TO_ANG = 0.529177210903
FORCE_RYBOHR_TO_EVANG = RY_TO_EV / BOHR_TO_ANG      # Ry/Bohr → eV/Å (≈25.711)
SQRT_EVA2AMU_TO_CM1 = 521.47083                     # √(eV/Å²/amu) → cm⁻¹


# ════════════════════════════════════════════════════════════════════════
# 1. 対称性適応変位基底
# ════════════════════════════════════════════════════════════════════════
def prepare_symmetry_modes(lattice, positions, numbers, symprec=1e-3):
    """各既約表現の部分空間全体 (m·d 本) の Cartesian 正規直交変位基底を作る。

    この部分空間は対称性だけで Φ 不変なので、ブロック対角化で得る振動数は
    音響和則を仮定せずに厳密。音響 (並進) は analyze 側で識別して除く。"""
    from .common import is_pseudo
    if any(is_pseudo(z) for z in numbers):
        raise AnalysisError(Msg("err.qe_pseudo"))
    sa = symmetry_adapted_blocks(lattice, positions, numbers, symprec)
    na = len(sa["numbers"])
    modes = []
    for blk in sa["blocks"]:
        rows = np.array([row for copy in blk["copies"] for row in copy])   # (m·d, 3N)
        q, r = np.linalg.qr(rows.T)
        rank = int(np.sum(np.abs(np.diag(r)) > 1e-6))
        basis = q[:, :rank].T
        modes.append(dict(label=blk["label"], d=int(blk["d"]), m=len(blk["copies"]),
                          activity=blk["activity"], tensor=list(blk["tensor"]),
                          basis=[b.reshape(na, 3) for b in basis]))
    ds = sa["dataset"]
    numbers = [int(z) for z in sa["numbers"]]
    return dict(lattice=sa["lattice"], frac=sa["positions"],
                cart=sa["positions"] @ sa["lattice"], numbers=numbers,
                species=[element_symbol(z) for z in numbers],
                masses=[atomic_mass(z) for z in numbers],
                spacegroup=ds.international, pointgroup=ds.pointgroup,
                centered=na < sa["conv_natoms"], na=na, modes=modes,
                rots=sa["rots"], trans=sa["trans"])


# ════════════════════════════════════════════════════════════════════════
# 2. pw.x 入力生成
# ════════════════════════════════════════════════════════════════════════
def pw_input(cfg, lattice, cart, species, prefix):
    """pw.x の SCF 入力 (力を出力する設定) を文字列で返す。"""
    kinds = list(dict.fromkeys(species))                  # 出現順に重複を除く
    absent = [k for k in ("pseudo_dir", "ecutwfc", "pseudopotentials") if k not in cfg]
    if absent:
        raise AnalysisError(Msg("err.qe_cfg_key", keys=", ".join(absent)))
    missing = [s for s in kinds if s not in cfg.get("pseudopotentials", {})]
    if missing:
        raise AnalysisError(Msg("err.qe_pp_missing", el=", ".join(missing)))
    L = ["&control", "   calculation = 'scf'", f"   prefix = '{prefix}'",
         f"   outdir = '{cfg.get('outdir', './tmp')}'",
         f"   pseudo_dir = '{cfg['pseudo_dir']}'",
         "   tprnfor = .true.", "   tstress = .false.", "/",
         "&system", "   ibrav = 0", f"   nat = {len(species)}",
         f"   ntyp = {len(kinds)}", f"   ecutwfc = {cfg['ecutwfc']}"]
    if cfg.get("ecutrho"):
        L.append(f"   ecutrho = {cfg['ecutrho']}")
    if cfg.get("smearing"):
        L += ["   occupations = 'smearing'", f"   smearing = '{cfg['smearing']}'",
              f"   degauss = {cfg.get('degauss', 0.01)}"]
    for k, v in cfg.get("system_extra", {}).items():
        L.append(f"   {k} = {v}")
    L += ["/", "&electrons", f"   conv_thr = {cfg.get('conv_thr', '1.0d-10')}",
          f"   mixing_beta = {cfg.get('mixing_beta', 0.7)}", "/", "ATOMIC_SPECIES"]
    from pymatgen.core import Element
    for s in kinds:
        L.append(f"  {s}  {atomic_mass(Element(s).Z):9.4f}  {cfg['pseudopotentials'][s]}")
    L.append("CELL_PARAMETERS (angstrom)")
    for v in lattice:
        L.append(f"  {v[0]:18.12f} {v[1]:18.12f} {v[2]:18.12f}")
    L.append("ATOMIC_POSITIONS (angstrom)")
    for s, r in zip(species, cart):
        L.append(f"  {s:3s} {r[0]:18.12f} {r[1]:18.12f} {r[2]:18.12f}")
    k = cfg.get("kpoints", [2, 2, 2])
    L += ["K_POINTS automatic", f"  {k[0]} {k[1]} {k[2]} 0 0 0"]
    return "\n".join(L) + "\n"


def _write_text(path, text):
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)


def displacement_amplitudes(basis, delta):
    """各基底ベクトル b (na, 3; 3N で規格化) に掛ける振幅 s [Å]: s·max_i|b_i| = delta。

    基底は 3N 次元で長さ 1 なので、そのまま delta 倍すると原子 1 個あたりの変位は
    delta/√(動く原子数) に減り、大きなセルでは力の数値ノイズに埋もれる。
    最大の原子変位が常に delta になるように規格化する。"""
    delta = _validated_amplitudes(delta, len(basis))
    return [float(s / np.max(np.linalg.norm(np.asarray(b), axis=1))) for s, b in zip(delta, basis)]


def _validated_amplitudes(delta, count):
    """中心差分の分母に使う振幅を、有限・正であることを確かめて配列にする。"""
    try:
        amps = np.broadcast_to(np.asarray(delta, float), (count,))
    except (TypeError, ValueError) as e:
        raise AnalysisError(Msg("err.qe_delta")) from e
    if not np.all(np.isfinite(amps) & (amps > 0)):
        raise AnalysisError(Msg("err.qe_delta"))
    return amps


def generate(cif=None, config=None, delta=0.01, labels=None, outdir="qe_fp",
             include="raman", symprec=1e-3, log=print):
    """変位構造の pw.x 入力群と manifest.json を生成し、manifest (dict) を返す。

    include: "raman" (ラマン活性のみ) / "optical_all" (全既約表現)。labels 指定時はそれを優先。"""
    _validated_amplitudes(delta, 1)
    lat, pos, num, _ = load_crystal(cif) if cif else graphite_2h()
    data = prepare_symmetry_modes(lat, pos, num, symprec)
    ensure_outdir(outdir)

    selected = [mi for mi, mo in enumerate(data["modes"])
                if (mo["label"] in labels if labels is not None
                    else include != "raman" or "Raman" in mo["activity"])]

    man = dict(delta=delta, lattice=data["lattice"].tolist(), cart=data["cart"].tolist(),
               species=data["species"], numbers=data["numbers"], masses=data["masses"],
               spacegroup=data["spacegroup"], pointgroup=data["pointgroup"],
               centered=data["centered"], na=data["na"], blocks=[])
    n_files = 0
    if config:
        _write_text(os.path.join(outdir, "scf_ref.in"),
                    pw_input(config, data["lattice"], data["cart"], data["species"], "ref"))
        man["ref"] = "scf_ref.out"
    for mi in selected:
        mo = data["modes"][mi]
        amps = displacement_amplitudes(mo["basis"], delta)
        binfo = dict(label=mo["label"], d=mo["d"], m=mo["m"], activity=mo["activity"],
                     tensor=mo["tensor"], basis=[b.tolist() for b in mo["basis"]],
                     amplitudes=amps, files=[])
        for a, b in enumerate(mo["basis"]):
            for sign, tag in ((+1, "p"), (-1, "m")):
                name = f"scf_{safe_label(mo['label'])}_{mi}_{a}_{tag}"
                if config:
                    _write_text(os.path.join(outdir, name + ".in"),
                                pw_input(config, data["lattice"],
                                         data["cart"] + sign * amps[a] * b, data["species"], name))
                binfo["files"].append(dict(a=a, sign=sign, out=name + ".out"))
                n_files += 1
        man["blocks"].append(binfo)
    with open(os.path.join(outdir, "manifest.json"), "w", encoding="utf-8") as f:
        json.dump(man, f, indent=1)

    log(t("qe.gen_title", sg=data["spacegroup"], pg=data["pointgroup"]))
    if data["centered"]:
        log(t("qe.centered", n=data["na"]))
    target = t("qe.target.selected") if labels is not None else (
        t(f"qe.target.{include}") if include in ("raman", "optical_all") else include)
    log(t("qe.target", target=target, n=len(selected)))
    for mi in selected:
        mo = data["modes"][mi]
        log(t("qe.block", label=f"{mo['label']:5s}", d=mo["d"], m=mo["m"], act=mo["activity"],
              n=len(mo["basis"])))
    log(t("qe.inputs", n=n_files, delta=delta) if config else t("qe.no_config"))
    log(t("qe.outdir", d=outdir))
    return man


# ════════════════════════════════════════════════════════════════════════
# 3. 力の回収とブロック対角化
# ════════════════════════════════════════════════════════════════════════
_FORCE_LINE = re.compile(
    r"atom\s+(\d+)\s+type\s+\d+\s+force\s*=\s*"
    r"(-?\d+\.\d+(?:[EeDd][+-]?\d+)?)\s+"
    r"(-?\d+\.\d+(?:[EeDd][+-]?\d+)?)\s+"
    r"(-?\d+\.\d+(?:[EeDd][+-]?\d+)?)")


def parse_forces(outfile, na):
    """pw.x 出力の最後の 'Forces acting on atoms' ブロックから力を読む [eV/Å] (na, 3)。

    verbosity='high' では総力のあとに非局所項・イオン項などの内訳ブロックが続き、
    同じ書式の行が並ぶ。総力だけを使うよう、各原子について最初に現れた行を採る。"""
    if not os.path.exists(outfile):
        raise AnalysisError(Msg("err.qe_no_output", f=outfile))
    with open(outfile, encoding="utf-8", errors="replace") as f:
        txt = f.read()
    marker = "Forces acting on atoms"
    seg = txt[txt.rfind(marker):] if marker in txt else txt
    end = seg.find("Total force")
    if end > 0:
        seg = seg[:end]
    found = {}
    for m in _FORCE_LINE.finditer(seg):
        i = int(m.group(1)) - 1
        if i not in found:
            found[i] = [float(m.group(k).replace("D", "E").replace("d", "e")) for k in (2, 3, 4)]
    if len(found) < na:
        raise AnalysisError(Msg("err.qe_forces", f=outfile, k=len(found), n=na))
    return np.array([found[i] for i in range(na)]) * FORCE_RYBOHR_TO_EVANG


def solve_block(basis, f_plus, f_minus, masses, delta):
    """1 つの既約表現ブロックの一般化固有値問題を解く。

    basis: list of (na,3); f_plus/f_minus: ±δ 変位での力 [eV/Å]; masses [amu]。
    delta: 変位振幅 [Å] (基底ベクトルごとのリスト、または共通のスカラー)。
    返り値: (ω[cm⁻¹] 昇順, 固有ベクトル list of (na,3) 正規化済み)。虚振動数は負で返す。"""
    nb = len(basis)
    na = basis[0].shape[0]
    amps = _validated_amplitudes(delta, nb)
    bf = np.array([b.reshape(-1) for b in basis])                          # (nb, 3N)
    phi_b = np.array([-(f_plus[a] - f_minus[a]).reshape(-1) / (2 * amps[a]) for a in range(nb)])
    k = bf @ phi_b.T                                                        # eV/Å²
    mass3 = np.repeat(np.asarray(masses, float), 3)
    m = bf @ (mass3[:, None] * bf.T)                                        # amu
    lam, v = generalized_eigh(k, m)
    omega = np.sign(lam) * np.sqrt(np.abs(lam)) * SQRT_EVA2AMU_TO_CM1
    eig = []
    for j in range(nb):
        u = v[:, j] @ bf
        n = np.linalg.norm(u)
        eig.append((u / n if n > 1e-12 else u).reshape(na, 3))
    return omega, eig


def analyze(outdir="qe_fp", manifest_path=None, write_modes=True, log=print):
    """pw.x 出力から力を回収し、振動数と固有ベクトルを求める。光学モードのリストを返す。"""
    mpath = manifest_path or os.path.join(outdir, "manifest.json")
    if not os.path.exists(mpath):
        raise AnalysisError(Msg("err.qe_no_manifest", p=mpath))
    with open(mpath, encoding="utf-8-sig") as f:
        man = json.load(f)
    delta, masses, na = man["delta"], man["masses"], man["na"]
    lattice, cart, numbers = np.array(man["lattice"]), np.array(man["cart"]), man["numbers"]

    mass3 = np.repeat(np.asarray(masses, float), 3)
    sqrt_mass3 = np.sqrt(mass3)
    tb = np.zeros((3, 3 * na))                     # 質量重み座標での剛体並進
    for ax in range(3):
        tb[ax, ax::3] = np.sqrt(masses) / np.sqrt(np.sum(masses))

    def trans_fraction(u):
        # K u = ω² M u の光学解は Σ_i M_i u_i = 0。
        # Cartesian 内積では質量差の大きい光学モードを並進と誤認する。
        uf = u.reshape(-1) * sqrt_mass3
        n2 = uf @ uf
        return float(np.sum((tb @ uf) ** 2) / n2) if n2 > 1e-12 else 0.0

    results = []
    for blk in man["blocks"]:
        basis = [np.array(b) for b in blk["basis"]]
        nb = len(basis)
        f_plus, f_minus = [None] * nb, [None] * nb
        for fr in blk["files"]:
            forces = parse_forces(os.path.join(outdir, fr["out"]), na)
            (f_plus if fr["sign"] > 0 else f_minus)[fr["a"]] = forces
        # 古い manifest (振幅の記録なし) では基底ベクトルに delta をそのまま掛けていた
        amps = blk.get("amplitudes", [delta] * nb)
        omega, eig = solve_block(basis, f_plus, f_minus, masses, amps)
        # 対称性より固有値は d 重に縮退して出るので、昇順に d 本ずつ束ねる
        d = blk["d"]
        idx = list(np.argsort(omega))
        groups = ([idx[k:k + d] for k in range(0, nb, d)] if d >= 1 and nb % d == 0
                  else [[j] for j in idx])
        for gp in groups:
            # 多重項全体の質量重み並進成分で識別 (成分の選び方によらない)。
            ft = float(np.mean([trans_fraction(eig[j]) for j in gp]))
            results.append(dict(label=blk["label"], omega=float(np.mean(omega[gp])),
                                deg=len(gp), kind="acoustic" if ft > 0.5 else "optical", tfrac=ft,
                                activity=blk["activity"], tensor=blk["tensor"],
                                eig=eig[gp[0]], eigs=[eig[j] for j in gp]))

    opt = sorted([r for r in results if r["kind"] == "optical"], key=lambda r: r["omega"])
    n_ac = sum(r["deg"] for r in results if r["kind"] == "acoustic")
    log("=" * 72)
    log(t("qe.an_title"))
    log("=" * 72)
    log("")
    log(t("qe.an_head", sg=man["spacegroup"], pg=man["pointgroup"], n=na, delta=delta))
    log("")
    log(t("qe.an_count", n=len(opt), dof=sum(r["deg"] for r in opt), ac=n_ac))
    log("")
    widths = (5, 10, 15, 9, 11)
    cols = t("qe.an_cols").split("|")
    log("  " + "".join(pad(c, w) for c, w in zip(cols, widths)) + cols[-1])
    log("  " + "-" * 70)
    for n, r in enumerate(opt, 1):
        w = f"{r['omega']:10.2f}" if r["omega"] >= 0 else f"{-r['omega']:9.2f}i"
        tc = ", ".join(r["tensor"]) if "Raman" in r["activity"] else "—"
        cells = (n, r["label"], w, t("qe.deg", d=r["deg"]), r["activity"])
        log("  " + "".join(pad(c, wd) for c, wd in zip(cells, widths)) + tc)
    soft = [r for r in opt if r["omega"] < 0]
    if soft:
        log("")
        log(t("qe.soft", n=len(soft)))

    if write_modes:
        mdir = ensure_outdir(os.path.join(outdir, "modes"))
        for n, r in enumerate(opt, 1):
            hdr = [f"Mode {n}: {r['label']}  {r['omega']:.2f} cm^-1",
                   f"Activity: {r['activity']}",
                   f"Raman tensor: {', '.join(r['tensor'])}" if r["tensor"] else "Raman: -"]
            for component, u in enumerate(r["eigs"]):
                suffix = f"_{chr(97 + component)}" if r["deg"] > 1 else ""
                name = f"mode_{n:02d}_{safe_label(r['label'])}{suffix}.xsf"
                write_xsf(os.path.join(mdir, name), lattice, cart, numbers, u,
                          hdr + [f"Component {component + 1}/{r['deg']}"], scale=0.4)
        log("")
        log(t("qe.xsf", d=mdir, n=sum(r["deg"] for r in opt)))
    return opt


# ════════════════════════════════════════════════════════════════════════
# 4. 解析機構の自己検証 (QE なしで数学を検証)
#    対称性を満たす合成 Φ を作り、全 3N 対角化とブロック法が一致するか確かめる
# ════════════════════════════════════════════════════════════════════════
def validate(cif=None, seed=0, symprec=1e-3, log=print):
    lat, pos, num, _ = load_crystal(cif) if cif else graphite_2h()
    data = prepare_symmetry_modes(lat, pos, num, symprec)
    na = data["na"]
    n3 = 3 * na
    rep = mech_rep(data["lattice"], data["frac"], data["numbers"],
                   data["rots"], data["trans"], symprec).real

    # 対称性を満たす合成力定数 Φs = mean_g D(g)^T Φ0 D(g) (+ 正定値シフト)
    rng = np.random.default_rng(seed)
    phi0 = rng.standard_normal((n3, n3))
    phi0 = 0.5 * (phi0 + phi0.T)
    phis = np.mean([D.T @ phi0 @ D for D in rep], axis=0) + 50.0 * np.eye(n3)
    phis = 0.5 * (phis + phis.T)

    # (a) 全 3N 質量重み対角化 (真値)
    mass3 = np.repeat(np.asarray(data["masses"], float), 3)
    lam = np.linalg.eigvalsh(phis / np.sqrt(np.outer(mass3, mass3)))
    w_full = np.sort(np.sign(lam) * np.sqrt(np.abs(lam)) * SQRT_EVA2AMU_TO_CM1)

    # (b) ブロック法: 合成「力」F(±δb) = ∓δ Φs b
    delta, w_block, total = 0.01, [], 0
    for mo in data["modes"]:
        amps = displacement_amplitudes(mo["basis"], delta)
        fp = [(-s * phis @ b.reshape(-1)).reshape(na, 3) for s, b in zip(amps, mo["basis"])]
        fm = [(+s * phis @ b.reshape(-1)).reshape(na, 3) for s, b in zip(amps, mo["basis"])]
        omega, _ = solve_block(mo["basis"], fp, fm, data["masses"], amps)
        w_block.extend(omega.tolist())
        total += len(mo["basis"])
    w_block = np.sort(w_block)

    dim_ok = total == n3
    diff = float(np.max(np.abs(w_full - w_block))) if dim_ok else float("inf")
    log(f"  {data['spacegroup']:10} ({data['pointgroup']:6}) "
        f"Σm·d={total}/{n3} {'✓' if dim_ok else '✗'}  "
        f"{t('qe.validate')}={diff:.2e} cm⁻¹ {'✓' if diff < 1e-4 else '✗'}")
    return dim_ok and diff < 1e-4
