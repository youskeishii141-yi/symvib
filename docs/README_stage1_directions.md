English | [日本語](README_stage1_directions.ja.md)

# Stage 1: Visualizing the "Direction" of Vibrational Modes (No DFT Required)

Given just a CIF, `factor_group_analysis.py` outputs, for each vibrational mode,
"which atoms move in which direction" as XSF (for VESTA).
**There are two kinds of certainty for the directions**, and they are explicitly distinguished in the output.

## The difference between "symmetry-exact" and "network guess"

### Symmetry-exact (exact)
For a mode whose irreducible representation occurs **only once** (m=1), the displacement direction is
**uniquely determined by group theory alone**. Neither DFT nor any model is needed; the output is directly
the direction of the true eigenvector. The T2g mode of diamond and the Raman modes of many
high-symmetry crystals fall into this category.

### Network guess (only when m>1; the true orientation is determined in Stage 2)
For modes whose irreducible representation occurs **more than once** (m>1), symmetry alone does not
determine the individual orientations (the freedom of linear combination within the 2m-dimensional space).
Therefore a force-constant matrix is built with a **central-force spring-network model with 1/d² decay**, diagonalized
within the subspace of that irreducible representation, and the resulting "guessed direction" is presented.

**Important caveats:**
- This is a **geometric guess** and is not necessarily the true eigenvector.
- The spring constants (1/d²) are not the actual bond strengths (which are determined by the electronic structure).
- If you need the true orientation, **Stage 2 (qe_frozen_phonon.py)** uses the QE forces and
  block diagonalization to determine it rigorously.
- However, since the model does reflect the geometric notion that "shorter bond = stronger bond," it is practically useful
  for getting a **rough idea** of the orientation (e.g., it correctly separates the two E2g modes of graphite into the
  interlayer shear mode and the in-plane G band).

## Reading the output

"Direction" column of the mode table:
- `symmetry-exact` … unique by group theory (exact)
- `guess (network)` … guess from the spring-network model because m>1
- `—`          … acoustic mode (translation; the notion of direction does not apply)

The same distinction is also recorded in the header of the XSF file (`# Direction:`):
- `symmetry-determined (exact)`
- `NETWORK-MODEL GUESS (m>1; 1/d^2 springs) - true eigenvector requires DFT`

## Usage
```
python factor_group_analysis.py mycrystal.cif
```
→ An XSF file for each mode is written to `modes_mycrystal/`.
Open it in VESTA and display the arrows via Properties → Vectors → Show vectors.

## Connecting to Stage 2
When you want the true orientations and frequencies of the m>1 modes:
```
python qe_frozen_phonon.py generate mycrystal.cif --config qe_config.json
# (run pw.x on the HPC system)
python qe_frozen_phonon.py analyze --dir qe_fp
```
By incorporating the QE forces, the symmetry-adapted basis of Stage 1 is elevated to the true eigenvectors.
