"""Independent reproductions of issues found in the prerelease audit."""
import json

import numpy as np
import pytest

from symvib.crystal import format_tensor
from symvib.molecule import translation_rotation_basis, vibrational_space
from symvib import frozen_phonon as fp


def test_mechanical_rep_requires_a_bijective_atom_mapping():
    from symvib.common import AnalysisError
    from symvib.crystal import mech_rep
    lattice = np.diag([4., 5., 6.])
    # Inversion plus a half-cell translation exchanges the two atoms.
    rep = mech_rep(lattice, np.array([[0., 0, 0], [.5, 0, 0]]), [6, 6],
                   np.array([-np.eye(3)]), np.array([[.5, 0, 0]]))
    assert np.allclose(rep[0], np.block([[np.zeros((3, 3)), -np.eye(3)],
                                        [-np.eye(3), np.zeros((3, 3))]]))
    # Overlapping sites must not silently produce a singular representation.
    with pytest.raises(AnalysisError):
        mech_rep(lattice, np.array([[0., 0, 0], [1e-6, 0, 0]]), [6, 6],
                 np.array([np.eye(3)]), np.zeros((1, 3)))


@pytest.mark.parametrize("isotropic", [False, True])
@pytest.mark.parametrize("periodic", [False, True])
def test_spring_hessian_matches_pair_energy(isotropic, periodic):
    from symvib.common import network_hessian
    cart = np.array([[0., 0, 0], [1., 0, 0]])
    lattice = np.diag([2., 4., 5.]) if periodic else None
    phi = network_hessian(cart, lattice=lattice, cutoff=1.01, k0=3., isotropic=isotropic)
    # E = (k/2) (u1-u0)^T B (u1-u0); periodic cell has two bonds per pair.
    block = (6. if periodic else 3.) * (np.eye(3) if isotropic else np.diag([1., 0, 0]))
    assert np.allclose(phi, np.block([[block, -block], [-block, block]]))
    assert np.allclose(phi @ np.tile(np.eye(3), (2, 1)), 0)


def test_network_cutoff_preserves_degenerate_irrep_copies():
    from symvib.common import refine_in_subspace
    from symvib.crystal import isotypic_basis
    from symvib.molecule import mechanical_matrices
    angles = np.arange(3) * 2 * np.pi / 3
    cart = np.column_stack([np.cos(angles), np.sin(angles), np.zeros(3)]) * 4 / np.sqrt(3)
    ops, chi = [], []
    for angle in angles:
        c, s = np.cos(angle), np.sin(angle)
        rz = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.]])
        for mirror in (False, True):
            r = rz @ (np.diag([1, -1, 1]) if mirror else np.eye(3))
            perm = np.argmin(np.linalg.norm((cart @ r.T)[:, None] - cart[None], axis=2), axis=1)
            ops.append((r, perm))
            chi.append(0. if mirror else 2 * c)
    rep = mechanical_matrices(ops, 3)
    basis = isotypic_basis(rep, np.array(chi), 2)
    rows = refine_in_subspace(basis, cart, np.ones(3), cutoff=4., sym_mats=rep)
    # Roundoff at d == cutoff must not select only some symmetry-equivalent bonds.
    for pair in rows.reshape(-1, 2, 9):
        q = np.linalg.qr(pair.T)[0]
        for r in rep:
            assert np.allclose(q @ q.T @ r @ q, r @ q, atol=1e-7)


@pytest.mark.parametrize("axis", [[1, 1, 0], [1, 2, 0], [0, 1, 1], [1, 1, 1]])
def test_linear_rotations_span_all_rigid_rotations(axis):
    axis = np.asarray(axis, float)
    axis /= np.linalg.norm(axis)
    cart = np.outer([-1.2, 0, 1.2], axis)
    masses = np.array([16., 12., 16.])
    trans, rot = translation_rotation_basis(cart, masses)
    rigid = np.array([(np.cross(e, cart) * np.sqrt(masses[:, None])).ravel()
                      for e in np.eye(3)])
    assert rot.shape == (2, 9)
    assert np.allclose(rigid, rigid @ rot.T @ rot, atol=1e-12)
    vib = vibrational_space(cart, masses)
    assert vib.shape == (4, 9)
    assert np.allclose(vib @ rigid.T, 0, atol=1e-12)
    assert np.allclose(vib @ trans.T, 0, atol=1e-12)


def test_tensor_format_converts_orthonormal_shear_coordinates():
    # T = outer((x+y)/sqrt(2), (x+y)/sqrt(2)); Txx = Tyy = Txy.
    assert format_tensor([[.5, .5, 0, 1 / np.sqrt(2), 0, 0]]) == (
        "[[a,a,·],[a,a,·],[·,·,·]]")


def test_frozen_phonon_optical_classification_uses_mass_metric(tmp_path):
    # An ASR-satisfying Hessian with one heavy and two light atoms.
    # Some true optical modes have > 50% *Euclidean* translation overlap.
    na = 3
    basis = np.eye(9).reshape(9, na, 3)
    masses = [200., 1., 1.]
    trans = np.tile(np.eye(3), (1, na)) / np.sqrt(na)
    phi = np.eye(9) - trans.T @ trans
    files = []
    delta = .01
    for a, b in enumerate(basis):
        for sign in (1, -1):
            name = f"force_{a}_{sign}.out"
            forces = (-sign * delta * phi @ b.ravel()).reshape(na, 3)
            lines = ["Forces acting on atoms (cartesian axes, Ry/au):"]
            for i, f in enumerate(forces / fp.FORCE_RYBOHR_TO_EVANG, 1):
                lines.append(f"atom {i} type 1 force = {f[0]:.14f} {f[1]:.14f} {f[2]:.14f}")
            (tmp_path / name).write_text("\n".join(lines), encoding="utf-8")
            files.append(dict(a=a, sign=sign, out=name))
    manifest = dict(delta=delta, na=na, masses=masses, lattice=(10*np.eye(3)).tolist(),
                    cart=[[0, 0, 0], [1, 2, 3], [3, 2, 1]], numbers=[80, 1, 1],
                    spacegroup="P1", pointgroup="1",
                    blocks=[dict(label="A", d=1, m=9, activity="Raman+IR", tensor=[],
                                 basis=basis.tolist(), files=files)])
    (tmp_path / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    optical = fp.analyze(str(tmp_path), write_modes=False, log=lambda *a: None)
    assert sum(r["deg"] for r in optical) == 6
    mass3 = np.repeat(masses, 3)
    lam = np.linalg.eigvalsh(phi / np.sqrt(np.outer(mass3, mass3)))[3:]
    assert np.allclose(sorted(r["omega"] for r in optical),
                       np.sqrt(lam) * fp.SQRT_EVA2AMU_TO_CM1, atol=1e-5)
