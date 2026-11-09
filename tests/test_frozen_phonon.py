"""QE 凍結フォノン: 合成した pw.x 出力で generate → analyze を通しで検証する。"""
import numpy as np
import pytest

from structures import crystals
from symvib import frozen_phonon as fp
from symvib.crystal import graphite_2h, mech_rep

CONFIG = {"pseudo_dir": "./pseudo", "ecutwfc": 40, "kpoints": [2, 2, 2],
          "pseudopotentials": {"Ti": "Ti.upf", "O": "O.upf", "C": "C.upf"}}
QUIET = dict(log=lambda *a: None)


@pytest.mark.parametrize("name", ["graphite", "rutile", "calcite"])
def test_validate(tmp_path, name):
    cif = None
    if name != "graphite":
        cif = str(tmp_path / f"{name}.cif")
        crystals()[name].to(filename=cif)
    assert fp.validate(cif, **QUIET)


def _fake_pw_output(path, forces_ev):
    """pw.x 風の出力 (力は Ry/Bohr)。verbosity='high' の内訳ブロックも付ける。"""
    f_ry = forces_ev / fp.FORCE_RYBOHR_TO_EVANG
    atom = "     atom {:4d} type  1   force = {:14.8f}{:14.8f}{:14.8f}"
    lines = ["     Forces acting on atoms (cartesian axes, Ry/au):", ""]
    lines += [atom.format(i + 1, *f) for i, f in enumerate(f_ry)]
    lines += ["     The non-local contrib.  to forces"]          # 読んではいけない内訳
    lines += [atom.format(i + 1, 9.0, 9.0, 9.0) for i in range(len(f_ry))]
    lines += ["", "     Total force =     0.000000     Total SCF correction =     0.000000"]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_generate_and_analyze_with_synthetic_forces(tmp_path):
    """対称な合成力定数 Φ から ±δ 変位の力を作り、解析した振動数が全対角化と一致する。"""
    man = fp.generate(cif=None, config=CONFIG, delta=0.01, outdir=str(tmp_path),
                      include="optical_all", **QUIET)
    assert (tmp_path / "scf_ref.in").exists()
    data = fp.prepare_symmetry_modes(*graphite_2h()[:3])
    na = data["na"]
    rep = mech_rep(data["lattice"], data["frac"], data["numbers"],
                   data["rots"], data["trans"]).real
    rng = np.random.default_rng(3)
    phi0 = rng.standard_normal((3 * na, 3 * na))
    phi = np.mean([D.T @ (phi0 + phi0.T) @ D for D in rep], axis=0) + 40 * np.eye(3 * na)
    # 音響和則: 剛体並進で力が出ないよう並進成分を除く (→ 音響 3 本は振動数 0)
    t = np.zeros((3, 3 * na))
    for ax in range(3):
        t[ax, ax::3] = 1 / np.sqrt(na)
    p_opt = np.eye(3 * na) - t.T @ t
    phi = p_opt @ phi @ p_opt
    for blk in man["blocks"]:
        for fr in blk["files"]:
            # 生成された pw.x 入力の変位 (最大の原子変位が δ) をそのまま使って力を合成する
            text = (tmp_path / fr["out"].replace(".out", ".in")).read_text(encoding="utf-8")
            pos = np.array([[float(x) for x in ln.split()[1:4]] for ln in
                            text.split("ATOMIC_POSITIONS (angstrom)")[1].strip().splitlines()[:na]])
            disp = pos - np.array(man["cart"])
            assert np.isclose(np.linalg.norm(disp, axis=1).max(), man["delta"])
            assert np.allclose(disp, fr["sign"] * blk["amplitudes"][fr["a"]]
                               * np.array(blk["basis"][fr["a"]]))
            _fake_pw_output(tmp_path / fr["out"], -(phi @ disp.reshape(-1)).reshape(na, 3))
    opt = fp.analyze(outdir=str(tmp_path), **QUIET)

    mass3 = np.repeat(man["masses"], 3)
    lam = np.linalg.eigvalsh(phi / np.sqrt(np.outer(mass3, mass3)))
    expected = np.sort(np.sqrt(np.abs(lam)) * fp.SQRT_EVA2AMU_TO_CM1)[3:]   # 音響 3 本を除く
    got = np.sort([r["omega"] for r in opt for _ in range(r["deg"])])
    assert len(got) == 3 * na - 3
    assert np.allclose(got, expected, atol=1e-3)
    assert len(list((tmp_path / "modes").glob("*.xsf"))) == sum(r["deg"] for r in opt)
    for mode in opt:
        assert len(mode["eigs"]) == mode["deg"]
        for u in mode["eigs"]:
            # Every exported partner must satisfy the full generalized eigenproblem.
            lam_mode = (mode["omega"] / fp.SQRT_EVA2AMU_TO_CM1) ** 2
            assert np.allclose(phi @ u.ravel(), lam_mode * mass3 * u.ravel(), atol=1e-5)


@pytest.mark.parametrize("delta", [0., -0.01, float("nan"), float("inf")])
def test_invalid_displacement_is_rejected_before_writing(tmp_path, delta):
    with pytest.raises(fp.AnalysisError):
        fp.generate(delta=delta, outdir=str(tmp_path), **QUIET)
    assert not list(tmp_path.iterdir())


def test_solve_block_rejects_invalid_amplitudes():
    basis = [np.array([[1., 0, 0]])]
    with pytest.raises(fp.AnalysisError):
        fp.solve_block(basis, basis, basis, [1.], [0.])


def test_missing_output_is_reported(tmp_path):
    fp.generate(cif=None, config=CONFIG, outdir=str(tmp_path), **QUIET)
    with pytest.raises(fp.AnalysisError, match="出力がありません"):
        fp.analyze(outdir=str(tmp_path), **QUIET)


def test_config_missing_pseudopotential(tmp_path):
    bad = dict(CONFIG, pseudopotentials={"O": "O.upf"})
    with pytest.raises(fp.AnalysisError, match="pseudopotentials"):
        fp.generate(cif=None, config=bad, outdir=str(tmp_path), **QUIET)
