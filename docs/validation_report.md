English | [日本語](validation_report.ja.md)

# Pre-release validation report (v2.1, 2026-11-08)

These are the results of checking the symmetry-based vibrational mode analysis (`symvib/`) against
external references (Bilbao Crystallographic Server, Crystallography Open Database, literature values, and phonopy,
an independent implementation). All scripts used for the checks are included in the repository and can be rerun.
The external data (the Bilbao character tables, the COD CIFs) are not included in the distribution. Retrieve them from
the respective official sites with `python tests/fetch_bilbao_tables.py` and `python tests/fetch_cod.py` (for the COD,
the script retrieves the same CIFs as the IDs recorded in `tests/data/cod/index.json`). If they have not been
retrieved, the tests that use them are skipped.

| Check | Content | How to run |
|---|---|---|
| Independent checks (invariants) | Irreducibility, orthogonality, degrees of freedom, 3 acoustic modes, fixed centre of mass, activity, Raman tensor components (projection operators), multiplicities from the correlation method, labels from the Bilbao point-group tables | `tests/invariants.py` (`pytest tests/test_invariants.py`) |
| 2a Sweep of all space groups | For every one of spglib's 530 Hall settings, a synthetic structure (2 orbits on general positions + special positions) is built and analysed. Agreement with the structure transformed to the standard setting, agreement with phonopy | `python tests/sweep_space_groups.py --phonopy` → [validation_space_groups.md](validation_space_groups.md) |
| 2b Real CIFs from the COD | Analysis of representative CIFs for 227 of the 230 space groups (IDs in `tests/data/cod/index.json`). Notes on reading, agreement with a synthetic structure on the same Wyckoff positions | `python tests/check_cod.py` → [validation_cod.md](validation_cod.md) |
| 2c Literature values | Numbers of Raman / IR / inactive modes for 33 substances (7 crystal systems) | `pytest tests/test_literature.py` |
| Molecules | Group-theoretical verification of the character tables of all 45 point groups, literature values for representative molecules, invariance under arbitrary rotation, permutation of the atom order and small distortions | `pytest tests/test_tables.py tests/test_molecule.py tests/test_point_groups.py`, `python tests/sweep_point_groups.py` (45/45) |

## 1. External references

- **Bilbao Crystallographic Server, Point Group Tables** (Bradley & Cracknell 1972; Altmann & Herzig 1994):
  `https://www.cryst.ehu.eus/cgi-bin/rep/programs/sam/point.py?sg=<space-group number>`. The character tables of the
  32 point groups (operation symbols and axis directions of the classes, Mulliken labels, assignment of x, y, z and of
  the quadratic forms) are retrieved with `tests/fetch_bilbao_tables.py`, saved to `tests/data/bilbao/point_tables.json`
  and used from there (not included in the distribution). For each irreducible representation of a crystal, the
  operations of the standard setting were mapped to the Bilbao classes by their axis directions, and the label,
  activity and tensor components were compared with those of the row whose characters match.
  - Bilbao's **SAM itself** (`nph-sam`: irreducible decomposition per Wyckoff position) requires human verification
    through Cloudflare Turnstile and could not be retrieved automatically. Instead, the multiplicities per Wyckoff
    position were verified independently with the correlation method (Bilbao characters × number of atoms left fixed
    by each operation) and with phonopy.
- **COD** (`https://www.crystallography.net/cod/`, public domain): using the search API
  (`result?format=json&space_group_number=N`), one entry was selected for each space group that is ordered (no
  fractional coefficients in the chemical formula), has ≤ 6 elements and a file size ≤ 120 kB (relaxed to 182 kB for
  No. 48 only), taking the one with the fewest atoms (`tests/fetch_cod.py`).
  No. 89 (P422), 106 (P4₂bc) and 168 (P6) could not be found under these conditions and are marked "not available".
- **phonopy 4.10** (`pip install phonopy`; an additional dependency for the tests only): a dynamical matrix is built
  from symmetrized random force constants, and the irreducible decomposition is compared with phonopy's irreps
  (Mulliken labels from character tables derived from Bilbao).
- Literature: the source for each substance is given in the comments of `tests/test_literature.py`.

## 2. Errors found and fixes

| # | Severity | Location | Description | Reproduction | Fix |
|---|---|---|---|---|---|
| 1 | Medium | `crystal.py` `_RAMAN_FORMS` | The Raman tensor forms of the E-type irreps (xz, yz) / (x²−y², xy) of 4/m, -4, 4, 6/m, 6, -6 and of E of 3, -3 were displayed in the same one-parameter form as for 4/mmm etc. Correctly they have two parameters (four for 3, -3) (Loudon, Adv. Phys. 13, 423 (1964)). The sets of components and the activities were correct | `tests/test_literature.py` / `test_invariants.py::test_raman_tensor_forms_match_loudon` (Eg of P4/m, etc.) | Removed the table; the forms are now computed exactly from the displacements of each mode and the row projection operator P = (1/h) Σ ⟨u\|D(g)\|u⟩ T(g) (`raman_tensor_forms`). The set of components is also determined from the diagonal elements of the projection operator |
| 2 | Medium | `frozen_phonon.py` `generate` | The displacement δ was multiplied onto basis vectors normalized in 3N dimensions, so the displacement per atom became δ/√(number of moving atoms), and in large cells the forces were buried in numerical noise | `tests/test_frozen_phonon.py::test_generate_and_analyze_with_synthetic_forces` (checks that the maximum atomic displacement = δ) | Uses amplitudes for which the maximum atomic displacement is δ (`displacement_amplitudes`) and records them in the manifest. The analysis side takes the central differences with the recorded amplitudes |
| 3 | Low | `molecule.py` `analyze_molecule` | Symmetrization of a distorted input stopped "once the number of operations no longer changed", so the operations did not form an exact group; residuals of ~1e-7 Å remained, and the subspaces of degenerate modes were off by about 1e-6 | Check applying a 0.02 Å distortion to a general molecule of every point group (`tests/invariants.check_molecule`) | The iteration that re-determines the operations with a strict tolerance and averages is continued until convergence (1e-12 Å) |
| 4 | Low | `cif_reader.py` `_symmetry_candidates` | When no symmetry operations were listed and only a short symbol was given ("P 21/c", "C 2/c"), the setting was always expanded with the b-axis setting, without looking at the cell angles (unique axis) | `tests/test_cif_reader.py::test_monoclinic_unique_axis_follows_the_cell_angles` | The unique axis is determined from the angle that is not 90° (α→a, β→b, γ→c) to narrow down the setting |
| 5 | Low (latent) | `crystal.py` `mulliken_labels` | The principal axis was taken as "the first of the highest-order proper rotation axes", so in D2d / S4 a C2' could be taken as the principal axis (no actual harm with spglib's operation order) | `tests/test_invariants.py::test_mulliken_labels_do_not_depend_on_operation_order` | If a rotoreflection axis (−4, −6) has a higher order, it is taken as the principal axis |
| 6 | Medium | `crystal.py` `symmetry_adapted_blocks` | For inputs whose lattice parameters are not exactly consistent with the detected symmetry (e.g. AsI₃, COD 9008098, R-3, rhombohedral angles 51.683/51.683/51.685°), the Cartesian rotation matrices were not orthogonal (8e-5), the force constants of the spring-network model did not commute with the symmetry operations, and the "direction guess" for multiplicity m>1 mixed different copies of the same irreducible representation (closure error 1e-5). The irreducible decomposition and the labels were correct; only the display of the degenerate components was affected | `tests/test_cod_regression.py` (9008098) | Before the analysis, the structure is idealized by averaging the lattice metric over the group and the coordinates over the orbits (`idealize_cell`). When the change exceeds 1e-5, a note is added to the report |
| 7 | Medium | `common.py` `network_hessian` | The spring-network model counted only the ±1 images of neighbouring cells, so in skewed primitive cells (body-centred, face-centred, rhombohedral) neighbouring atoms within the cutoff were missed and the force constants did not commute with the symmetry operations (COD 1010030, I-42m: closure error 0.06). As with 6, only the direction of the degenerate components for m>1 is affected | `tests/test_cod_regression.py` (1010030) | The required range of images is determined from the cutoff and the interplanar spacings of the lattice |
| 8 | Low | `molecule.py` `find_point_group_operations` | When the reference atom for the search for symmetry operations was chosen only by "fewest atoms of the same species", an atom right next to the centre of mass could be chosen; with a distortion of about 0.02 Å the estimate of the rotation then went far off and the symmetry was lost (a C3 molecule was identified as C1) | `tests/test_molecule.py::test_distorted_molecule_with_atom_near_centre_of_mass` | The reference is chosen from the atoms of the outer shell (distance from the centre of mass ≥ 40% of the maximum) |
| 9 | Medium (GUI) | `gui.py` `_build_scene` | For crystals consisting only of element pairs that are not in the VESTA bond table (AsI₃, etc.) there was not a single bond; passing an empty Line3DCollection to `add_collection3d` made matplotlib raise an exception, and the analysis results were not shown on screen (found while checking that the exe starts) | `tests/test_gui.py::test_crystal_without_table_bonds_builds_a_scene` | No line collection is created when there are no bonds |

Not errors, but conventions that were stated explicitly:
- **-4m2 family (No. 115–120)**: the Bilbao point-group table is written in the -42m orientation (2-fold axis ∥ [100]).
  In -4m2 the 2-fold axis is along [110], so the code's B1 (symmetric with respect to the 2-fold axis: Mulliken's rule)
  is xy and B2 (z) is x²−y², which looks swapped relative to B1 = x²−y² in the Bilbao table. Group-theoretically,
  B2 ∋ z is the same in both settings; this is a difference of convention.
  It has not been confirmed directly with SAM (see above).
- The inactive mode of **h-BN** (P6₃/mmc, B 2c, N 2d) is written as B1g in the literature (Geick 1966), but in the
  Bilbao convention (subscript 1 = symmetric with respect to the 2-fold axis along [120]) it is B2g, the same as in
  graphite. This was confirmed by computing the characters directly
  (`tests/test_literature.py::test_b2g_of_graphite_and_hbn_follow_bilbao_convention`).
- The Bilbao point-group tables for **monoclinic** groups take z as the 2-fold axis, but the standard setting for
  crystals is unique axis b. The code uses the standard setting (Ag ∋ xz, Bg ∋ xy, yz) (the same as the convention in
  the literature).
- The Raman tensor components of molecules are expressed in the standard orientation of the character table
  (z = principal axis); a note in the report states that they do not depend on the orientation of the input
  coordinates.

## 3. Results

### 3a. Synthetic structures for the 230 space groups (530 Hall settings)
`tests/sweep_space_groups.py --phonopy` (run with the final code, about 1 hour). The results for each setting are in
[validation_space_groups.md](validation_space_groups.md).

- **Independent checks and agreement with the standard setting: 530 / 530 settings and 230 / 230 space groups pass.**
  For each setting, we verified 3N degrees of freedom, 3 acoustic modes, a fixed centre of mass, closure of each block
  under D(g), |χ|²/h ∈ {1, 2}, consistency of the labels (orthogonality), activity (inner products with the vector
  representation and the symmetric rank-2 tensor representation), the Raman tensor components (projection operators),
  the multiplicities per Wyckoff position from the correlation method, the labels of the Bilbao point-group tables and
  the assignment of functions (x, y, z, quadratic forms), as well as agreement with the result of analysing the same
  structure after transforming it to the standard setting with spglib.
- Settings whose **orientation differs** from the Bilbao tables (recorded as label swaps, not errors): the 6 settings
  of the -4m2 family (No. 115–120; see Section 2).
- **Comparison with phonopy: the irreducible decomposition agrees for 451 settings; 78 settings could not be compared
  because phonopy could not assign labels ("Rotation symbols not found" for the P3₁2 types, groups in which A/B or
  one-dimensional complex representations cannot be distinguished, etc.); in 1 setting (No. 193 P6₃/mcm)
  B1g/B1u ↔ B2g/B2u were swapped.** For No. 193, the code's labels were confirmed to be correct by matching the
  characters directly against the Bilbao point-group table (subscript 1 = symmetric with respect to the 2-fold axis
  along [100]), and the discrepancy was judged to be due to a difference in how phonopy chooses the orientation of the
  2-fold axes (for No. 194 P6₃/mmc, phonopy could not be compared).
- For orthorhombic groups (P222, Pmmm types, etc.), B1/B2/B3 follow the axis order determined by spglib's
  standardization (a < b < c). If the axes of the input are chosen differently, the labels may be swapped relative to
  the literature (when the input is not in the standard setting, the report adds the note "labels and tensor
  components are expressed in the axes of the standard setting").

### 3b. Real CIFs from the COD
`tests/check_cod.py` (run with the final code, about 2.5 hours). The results for each CIF are in [validation_cod.md](validation_cod.md).

- **227 / 227 entries (227 of the 230 space groups) pass, with 0 problems.** Each CIF was taken through reading →
  analysis → independent checks (the same as in 3a); in addition, a synthetic structure with random free coordinates
  on the same Wyckoff positions of the same space group (a different element for each orbit) was analysed, and the
  multiplicities of the irreducible decomposition were confirmed to agree (a check that the CIF reading and the
  handling of Wyckoff positions are correct). The number of atoms ranges from 1 to 1344 (unit cell).
- 73 entries received notes on reading: 50 reinterpretations of rounded fractional coordinates (0.333 → 1/3, etc.),
  partially occupied sites treated as fully occupied, and sites sharing the same position merged into pseudo-atoms.
  The notes appear in the "Notes on reading" column of the table and in the GUI / CLI reports.
- In 49 entries, a **higher symmetry** than the space group given in the CIF was detected (the rows of the table where
  "Detected space group" differs from the COD number). All of them give the same result even with symprec = 1e-6 Å;
  this is because the coordinates lie exactly on high-symmetry positions (idealized structures, hypothetical ordered
  structures of hydrides and intermetallic compounds, disordered structures with partially occupied sites treated as
  fully occupied (e.g. Bi₂O₃ COD 1528292 Fm-3 → Pm-3, cell a/2), etc.). No CIF ended up with a **lower** symmetry than
  declared (`tests/test_cod_regression.py` also confirms for 135 entries that the order of the point group does not
  decrease).
- The 135 small entries (12 atoms or fewer in the unit cell) also run as ordinary tests (`tests/test_cod_regression.py`).

### 3c. Literature values (33 substances)
All items of `tests/test_literature.py` agree. Substances: diamond, NaCl, CsCl, CaF₂, ZnS, SrTiO₃,
spinel MgAl₂O₄, YAG, pyrite, Cu₂O (cubic); rutile, anatase, zircon, scheelite, chalcopyrite, PbTiO₃
(tetragonal); α-quartz, corundum, calcite, LiNbO₃, Bi₂Se₃, ilmenite (trigonal); graphite, h-BN, 2H-MoS₂, ZnO
(hexagonal); forsterite (Pbnm → Pnma reinterpretation), KTiOPO₄ (orthorhombic); baddeleyite ZrO₂, β-Ga₂O₃,
diopside (monoclinic); kyanite (triclinic).

Invariance with respect to how the input is given: non-standard settings (Pbnm / Pnma, R / H axes, origin 1 / 2,
unique axis a/b/c), primitive cell / superlattice, permutation of the atom order and the way coordinates are written
(1/3, 0.333, 0.3333) are checked by `tests/sweep_space_groups.py` (each of the 530 settings compared with the standard
setting), `tests/test_crystal.py` and `tests/test_cif_reader.py`.

### 3d. Molecules
- Character tables of the 45 point groups (real + complex): verified that they give projection operators on the
  regular representation, and that the IR/Raman assignments and components agree with the decompositions of the
  vector representation and the symmetric rank-2 tensor representation (`tests/test_tables.py`).
- Literature values: water, formaldehyde, NH₃, CH₄, CO₂, HCN, C₂H₂, ethylene, benzene, SF₆, BF₃, XeF₄, ethane,
  allene, ferrocene (D5h / D5d), B(OH)₃, B₁₂H₁₂²⁻, C₆₀, PCl₅, XeF₂, N₂O, adamantane, Ni(CO)₄,
  Mo(CO)₆, tropylium, IF₇, spiropentane, formic acid, HNCS and others (`tests/test_molecule.py`,
  `tests/test_point_groups.py`).
- For molecules with atoms on general positions in all 45 point groups, verified that the point group and the
  decomposition do not change under arbitrary rotation, permutation of the atom order and a 0.02 Å distortion, and
  that the independent checks pass (the check script was run when this report was prepared; the exact symmetrization
  of the distorted structure was improved by fix 3).

## 4. Refactoring
- `crystal.py`: removed the table of Raman tensor forms (`_RAMAN_FORMS`, 13 lambda entries) and replaced it with a
  computation by projection operators (`sym2_rep`, `raman_tensor_forms`, `format_tensor`). Moved the acoustic /
  optical separation into `split_acoustic_optical`. Made `tensor_components` use exact projection operators. Removed
  duplicate imports.
- `frozen_phonon.py`: consolidated the displacement amplitudes in `displacement_amplitudes`.
- i18n: added the new display strings (notes on the Raman tensors, the form of each component, the orientation of
  molecular tensors) in English and Japanese.
- Tests: consolidated the independent checks in `tests/invariants.py`, which each test reuses.

## 5. Remaining concerns and unverified items
- A direct comparison with Bilbao SAM (decomposition per Wyckoff position) has not been carried out (Turnstile). The
  correlation method, phonopy and literature values were used instead. The B1/B2 convention of the -4m2 family has
  also not been confirmed with SAM.
- The default `symprec = 1e-3 Å`: among the 227 COD entries, all 49 in which a higher symmetry than declared was
  detected have coordinates exactly on high-symmetry positions (the same result with 1e-6 Å), and there was no case in
  which the tolerance either raised or lowered the symmetry. 1e-3 Å is reasonable for the precision of the coordinates
  of experimental structures (around 1e-4), but old CIFs with low precision may require `--symprec 1e-2`. Conversely,
  slightly distorted low-symmetry phases (pseudo-cubic perovskites, etc.) can be rounded up to a higher symmetry even
  at 1e-3 Å, so the detected space group is shown at the top of the report.
- `cif_reader._coordinate`: three-digit decimals without an uncertainty (0.333, 0.667, 0.167 …) are reinterpreted as
  fractions (1/3 …). If a refined free parameter happens to be written as 0.333, this moves it by about 0.0003 (below
  symprec, so it does not affect the symmetry; we judged that a loss of symmetry in old CIFs with many values rounded
  from 1/3 does more actual harm). A note appears in the report.
- The spring-network model (direction guess for m>1) and the choice of degenerate components are display rules, not
  physical predictions, and the report states this explicitly.
- Tensor components of molecules are expressed in the standard orientation of the point group (no transformation to
  the input coordinate system is performed).
- Since the GUI contains no mathematical processing, this validation only checked that it starts and that it
  displays the results.
