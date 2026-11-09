English | [日本語](README_frozen_phonon.ja.md)

# Force-Based Γ-Point Block Frozen Phonon (Quantum ESPRESSO Interface)

Using the symmetry-adapted displacement basis obtained with `factor_group_analysis.py` (factor-group analysis),
this tool obtains, from the forces of pw.x SCF calculations, the **Γ-point frequencies of the Raman-active modes** and
the **eigenspaces in which the mixing between copies of the same symmetry species is resolved**.
The accuracy depends on the electronic-structure method adopted, the convergence of the forces, and the finite-displacement amplitude. The orientation of the basis
within a degenerate eigenspace and the overall sign remain arbitrary even after the QE forces are incorporated.
The computational cost and the relative merits compared with `ph.x` (DFPT) depend on the target structure, the selected modes, and the calculation conditions.

## Requirements
```
pip install numpy spglib spgrep pymatgen
```
Place it in the same folder as `factor_group_analysis.py`.

## Usage (two phases)

### Phase 1: Generate pw.x inputs for the displaced structures
```
python qe_frozen_phonon.py generate mycrystal.cif --config qe_config.json
```
- `scf_*.in` (pw.x inputs with displaced atoms) and `manifest.json` are written to `qe_fp/`.
- By default, only Raman-active modes are included. For all optical modes, use `--include optical_all`.
- For specific irreducible representations only, use `--labels Eg A1g`.
- The displacement amplitude is set with `--delta 0.01` (Å, default). δ is "the displacement of the atom that moves the most"
  (because the symmetry-adapted basis is normalized in 3N dimensions, multiplying it directly by δ would reduce the displacement per atom
  to 1/√(number of moving atoms), burying it in noise for large cells. From v2.1 on, the basis is normalized so that the maximum atomic displacement = δ,
  and the amplitudes used are recorded in `amplitudes` in `manifest.json`).

Run these through pw.x on the HPC system (each `scf_XXX.in` → `scf_XXX.out`).

### Phase 2: Collect the forces and analyze the frequencies
```
python qe_frozen_phonon.py analyze --dir qe_fp
```
- Reads the forces from each pw.x output and diagonalizes the force-constant blocks with the mass metric included.
- Displays the list of Raman peaks (irreducible representation, frequency, degeneracy, tensor components).
- Writes the eigenvector XSF files (for VESTA) to `qe_fp/modes/`. For degenerate modes, all components
  are saved, e.g. `mode_02_Eg_a.xsf`, `mode_02_Eg_b.xsf`.
  In the Python return value, `eigs` holds all components, and the legacy `eig` is the first component.

### Self-validation (no QE required)
```
python qe_frozen_phonon.py validate mycrystal.cif
```
Checks, using synthetic symmetric force constants, whether the block method agrees with full 3N diagonalization.

## Important caveats
- This implementation obtains the analytic part of the force constants at the Γ point under periodic boundary conditions. It
  **does not include the direction-dependent LO–TO splitting** of polar insulators. That calculation additionally requires Born effective charges, the high-frequency dielectric constant, and
  the non-analytic correction ([Phonopy formulation](https://phonopy.github.io/phonopy/formulation.html)).
- For the block decomposition to be correct, the structure and the mass configuration, as well as the computed force constants, must preserve the adopted
  symmetry. The acoustic/optical distinction is an approximate criterion based on the mass-metric overlap with translations,
  and needs to be checked when the force errors are large or when there is mixing with zero-frequency optical modes.
- **Always input a relaxed structure (vc-relax).** For a structure that is not at a stationary point,
  the "frequencies" have little physical meaning (mathematically, central differences still yield a curvature).
- For force accuracy, set **conv_thr to about 1e-10**.
- δ of 0.01–0.02 Å is recommended (too small gives numerical noise; too large introduces anharmonic contributions).
- To extend this to anharmonic scans (e.g., double wells of soft modes), vary δ over a wide range
  and fit V(Q) with a polynomial (a separate tool, equivalent to ModeMap / CRYSTAL).

## Estimated computational cost
Only ±δ displacements for the total number of Raman branches (= 2 × the number of branches SCF runs) plus a single reference.
Example: calcite Raman A1g+4Eg → 9 branches → 18 displaced SCFs + 1 reference.
The blocking itself is equivalent to diagonalizing the full matrix as long as the symmetry is preserved, but it does not remove the
approximation errors of finite differences or DFT.
