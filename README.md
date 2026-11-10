English | [日本語](README.ja.md)

# SymVib

**Symmetry-based vibrational analysis for crystals and molecules**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24085930.svg)](https://doi.org/10.5281/zenodo.24085930)

A set of Python tools for assigning and visualizing the vibrational modes of
crystals and molecules on the basis of symmetry (group theory). It is specialized
for understanding "which vibrational modes exist and in which direction each atom
moves" without performing a DFT phonon calculation. It is intended as an aid to mode
assignment in Raman / infrared spectra and for teaching group theory and vibrational
spectroscopy.

It does not compute absolute vibrational frequencies (that is the role of quantum-chemical
and first-principles calculations such as Gaussian / ORCA / Quantum ESPRESSO). This tool
deals with what symmetry determines: "the types of modes, their degeneracy, their activity,
and their displacement patterns".

- Author: Yosuke Ishii (石井 陽祐), Nagoya Institute of Technology
- License: [MIT License + Commons Clause](LICENSE) — free to use, modify and redistribute,
  including for corporate R&D, but **selling the software (or modified versions) is prohibited**
  (see "[License](#license)" for details).
- Contact: [GitHub Issues](https://github.com/youskeishii141-yi/symvib/issues)
  or e-mail (see "[Contact and feature requests](#contact-and-feature-requests)").

## Design philosophy: "symmetry-exact" and "network guess"

For each vibrational mode, if its symmetry species (irreducible representation) occurs only
once in the space (m=1), the displacement subspace is uniquely determined by group theory
alone (**exact**). The orientation of the degenerate components and the overall sign are
arbitrary. If the same symmetry species occurs more than once (m>1), symmetry alone does not
fix the direction, and a freedom of linear combination remains. In this case, an
**estimated** direction from a 1/d² spring-network model is presented, explicitly labelled as
such. If the true eigenvectors are needed, one proceeds to a quantum-chemical calculation;
this is the division of roles.

A distinctive feature is that this "difference between symmetry-adapted modes and true
eigenvectors" can be explored within a consistent framework for both molecules and crystals.

### Note: what a quantum-chemical calculation determines, and what remains

Normal modes are the eigenvectors of the mass-weighted force-constant matrix (Hessian).
Group theory determines only the **form**, namely that "this matrix splits into blocks, one
for each irreducible representation"; it says nothing about the **values** of the matrix
elements.

- **m=1**: The block of a symmetry species with degeneracy g (1 for A and B, 2 for E, 3 for T)
  is g×g, and a symmetry-preserving harmonic dynamical matrix is eigenvalue × identity matrix
  within this space. The eigenspace is determined independently of the values of the force
  constants.
- **m>1**: The full block is (m·g)×(m·g). For irreducible representations of real type,
  aligning the bases of the partners reduces this to an m×m mixing problem (complex-conjugate
  pairs are treated in real space). How the m copies are mixed is determined by the actual
  values of the matrix elements. Quantum-chemical and first-principles calculations supply
  these values, so within the harmonic approximation the mixing freedom disappears and the
  eigenvectors are determined. The 1/d² spring-network model substitutes a crude approximation
  (central forces only, no bending stiffness) for these values, and is therefore treated as
  an "estimate".

Example: in H₂O with Γvib = 2A₁ + B₂, B₂ (antisymmetric stretch) has m=1, and its displacement
pattern is determined by symmetry alone. A₁ has m=2, and the two SALCs "symmetric stretch" and
"bend" mix according to the force constants. How much they mix cannot be known without a
calculation.

However, even after a calculation, the following points remain.

1. **Orientation within a degenerate symmetry species (E, T, etc.)**: Within a degenerate
   space, any choice of orthogonal basis gives equally correct eigenvectors. This is genuinely
   arbitrary physically and is not removed by any calculation (the output of a calculation
   program is simply one of them, picked by chance). This is a different kind of freedom from
   the "mixing between different copies" for m>1.
2. **The overall sign (±)** is always arbitrary.
3. **Dependence on the level of calculation**: The mixing coefficients change with the choice
   of functional, basis set, etc. When two modes of the same symmetry species have close
   frequencies, the way they mix becomes sensitive to small differences in the force constants
   (if they happen to be accidentally degenerate, it again becomes arbitrary between those two).
4. **Beyond the harmonic approximation**: Once anharmonicity, Fermi resonance, and crystal-field
   or solvent effects are included, the very picture of "eigenvectors" is an approximation.
   Moreover, eigenvectors cannot be measured fully and directly in experiments; they can only be
   verified partially, e.g. from isotope shifts or neutron-scattering intensities.

In this repository, for crystals, the QE frozen phonon (`qe_frozen_phonon.py`) corresponds to
this "proceed to a quantum-chemical calculation" route: it obtains the force constants from
forces and determines even the mixing for m>1.

| What determines it | Target | Treatment |
|---|---|---|
| Symmetry | Mode types, degeneracy, activity; displacement patterns for m=1 | Exact in this tool |
| Values of the force constants | Mixing of SALCs for m>1 | Estimated by the spring-network model → determined by first-principles calculation if needed |
| (Nothing determines it) | Orientation within a degenerate space; overall sign | Arbitrary with either method |

## Included tools

| Entry point | Target | Overview |
|---|---|---|
| `vibration_gui.py` (distribution version `SymVib.exe`) | Molecules and crystals | GUI. Mode list and 3D view (displacement arrows, animation), saving XSF/reports. |
| `molecule_vibration_analysis.py` | Isolated molecules | Point-group-based molecular vibration analysis (CLI). Supports 45 point groups. |
| `factor_group_analysis.py` | Crystals (periodic systems) | Factor-group analysis of space groups (CLI). Decomposition into irreducible representations at the Γ point. |
| `qe_frozen_phonon.py` | Crystals (QE interface) | Force-based Γ-point frozen phonon (CLI; Quantum ESPRESSO). |

The core of the analysis is collected in the package `symvib/`, and the CLI and GUI
call the same functions (`analyze_molecule` / `analyze_crystal`).

### 1. Molecular vibration analysis

Given a CIF or xyz input, it automatically determines the molecular point group and labels the
internal vibrational modes by irreducible representation (Mulliken symbols). It outputs the
Raman / infrared activity of each mode, the polarizability-tensor components, and the
displacement vector of each atom in XSF format (for VESTA).

- **Supports 45 point groups**: nearly complete coverage of the Cn / Cnv / Cnh / Sn / Dn / Dnd /
  Dnh series for n = 2–8, the cubic groups T / Th / Td / O / Oh, the icosahedral groups I / Ih,
  and the linear molecules C∞v / D∞h.
- Linear molecules (infinite groups) are decomposed exactly with true continuous-group projection
  operators, obtained by numerically integrating over the continuous family of rotation angles
  (no approximation by a finite subgroup).
- Groups with complex irreducible representations (Cn, Cnh, Sn, T, Th) are handled correctly
  with complex characters.
- Small distortions in experimental or computed structures (up to 0.1 Å by default) are
  symmetrized exactly to the detected point group before analysis (adjustable with `--tol`; the
  symmetrization is stated explicitly in the report).

For details, see [docs/README_molecule.md](docs/README_molecule.md).

### 2. Crystal factor-group analysis

Given a CIF input, it determines the space group and decomposes the lattice vibrations at the
Γ point (k=0) into irreducible representations. It supports automatic assignment of Mulliken
labels, polarized Raman selection rules, and XSF output of the displacement pattern of each
mode. For modes with m>1, estimated directions from a network model are presented
([docs/README_stage1_directions.md](docs/README_stage1_directions.md)).

### 3. QE frozen phonon

In combination with Quantum ESPRESSO, it obtains Γ-point phonons based on force calculations.
It provides a workflow of generating displaced inputs (`generate`), analysis from forces
(`analyze`), and symmetry validation (`validate`). See
[examples/qe_config_example.json](examples/qe_config_example.json) for a configuration example
and [docs/README_frozen_phonon.md](docs/README_frozen_phonon.md) for details.

## GUI version

### For users (distribution version)

Download `SymVib.zip` from
[Releases](https://github.com/youskeishii141-yi/symvib/releases), extract it, and
double-click `SymVib.exe` in the folder. No Python installation is required. The zip
also contains this tool's license (`LICENSE`), the documentation (`README.md` / `README.ja.md`),
and the license notices of the bundled libraries (`THIRD_PARTY_LICENSES.txt`).

- Requirements: Windows 10 / 11 (64-bit). No installation or administrator rights are needed.
- If "Windows protected your PC" (SmartScreen) appears the first time you start it, choose
  "More info" → "Run anyway" (this is shown for apps without a digital signature).
- If it does not start, see "[Troubleshooting](#troubleshooting)" below.

1. Choose a `.xyz` or `.cif` file with "Open…" (molecule/crystal is selected automatically from
   the extension; for a molecular CIF, switch to "Molecule" and click "▶ Analyze"). You can also
   try water and graphite from "Demo".
2. Selecting a mode in the list on the left shows the atoms, bonds and displacement arrows in 3D
   on the right. **Left-drag to rotate, wheel to zoom, right-drag to pan**. The orientation of
   the crystal axes (a, b, c; x, y, z for molecules) is shown at the lower left. For crystals,
   all atoms inside the unit cell of the input file (including those on the boundary) are shown,
   and image atoms move in the same way as the original atoms. For degenerate modes, switch
   between a/b/c with "Component:", and make the mode vibrate with "Animation".
   As with the VESTA defaults, bonds are determined from a table of maximum bond lengths for each
   pair of elements (the VESTA defaults bundled with pymatgen). For pairs not in the table
   (metal–metal pairs such as Al–Al), no bond is drawn. For molecules, pairs not in the table
   (H–H, etc.) are also treated as bonded if within 1.25 × the sum of the covalent radii.
   The contributing sites in the list (the "sites" part of the "Raman tensor / sites" column) are
   the sites of the atoms that can move in that mode (element + Wyckoff letter, e.g. `Al(c)`).
   This is exact information determined by symmetry alone.
3. With "Save XSF/VESTA…", choose the format ("VESTA only" / "XSF only" / "Both") and write out
   all modes; with "Save report…", save a text report. `.vesta` files can be opened directly in
   VESTA and already include thin gray arrows (the midpoint of each arrow at the atom's center,
   ±0.67 Å from the center at the maximum displacement, ±0.53 Å for molecules) and the bond
   settings.
4. The "Language / 言語" control at the upper right switches between English and Japanese (the
   default is English; the selected language is also used at the next start-up). Both the screen
   and the reports are switched.

For the components a, b, c of degenerate modes (E, T, etc.), symmetry alone does not determine
the orientation (any re-choice within the degenerate space is equally correct). The components
shown are chosen by the following common rules, independent of the material (there are no
material-specific settings or special cases).

1. In the order of axes 1, 2, 3 (a, b, c of the input cell for crystals; x, y, z for molecules),
   take the motion within the degenerate space that maximizes "the fraction of the displacement
   of all atoms along that axis". If it is nearly the same direction as an already chosen
   component, choose again among the remainder orthogonal to it.
2. When the axes alone do not decide (e.g. when all atoms move only along c, as in the Eg mode of
   rutile), the remaining components are determined by ordering the atoms by their coordinates.
3. The sign is fixed so that "the sum of the displacements along the axis that determined the
   component is positive".

These rules do not depend on the order of the atoms (how the CIF is written) and always give the
same components for the same structure. For modes that can move along the axes, such as E1/E1u
of graphite and ZnO or T1u of cubic crystals, component a is exactly parallel to axis a and
component b to axis b. For other modes (such as E of quartz), it is "the motion best aligned with
axis a", and the arrows are not necessarily parallel to the axis. The components are not
necessarily orthogonal to each other (a and b of a hexagonal crystal are at 120°). The initial
orientation of the 3D view is also chosen by computation, independent of the crystal system, so
that all three cell axes are easy to see.

The colors in the list indicate activity (red = Raman, blue = IR, purple = both, gray = inactive),
and the "Direction" column distinguishes "symmetry-exact / guess (network)".

### Troubleshooting

If SymVib does not start and shows a message such as "SymVib could not start" or "DLL load failed",
try the following in order. When the startup fails, SymVib automatically checks the extracted
files for missing or damaged ones and names them in the dialog.

1. **Extract the zip again, completely**: delete the SymVib folder, right-click the downloaded
   `SymVib.zip` → "Extract All", wait until the extraction has finished, and then start
   `SymVib.exe`. Do not run the exe directly from the window that shows the contents of the zip.
   A short location such as `C:\SymVib` is recommended.
   Before extracting, right-clicking `SymVib.zip` → "Properties" → ticking "Unblock" at the bottom
   → "OK" makes SmartScreen warnings and loading restrictions less likely.
   **On some PCs, "Extract All" stops partway** (seen on a managed PC; if folders such as
   `numpy.libs` or `_tcl_data` are missing from `_internal` after extraction, it stopped partway).
   In that case open "PowerShell" from the Start menu and run the following lines one by one to
   extract and start it (example for a zip saved in `Downloads`, extracted to `SymVib` in your user
   folder). If something fails, the error is shown.
   ```powershell
   Unblock-File "$HOME\Downloads\SymVib.zip"
   Expand-Archive -Path "$HOME\Downloads\SymVib.zip" -DestinationPath "$HOME\SymVib" -Force
   & "$HOME\SymVib\SymVib\SymVib.exe"
   ```
   You can check the extraction with the self-test below (`RESULT: OK` at the end means it is fine).
   ```powershell
   Start-Process "$HOME\SymVib\SymVib\SymVib.exe" -ArgumentList "--self-test", "$HOME\symvib_check.txt" -Wait
   Get-Content "$HOME\symvib_check.txt"
   ```
2. **Check whether security software quarantined files**: in Windows Security → "Virus & threat
   protection" → "Protection history", look for quarantined SymVib files (such as `.dll` files).
   If you use other security software, check its quarantine. On a PC managed by a university or
   company, running apps may be restricted; ask the administrator.
3. **Install the Microsoft Visual C++ Redistributable (x64)**: the required files are bundled with
   SymVib, but installing it solves the problem on some PCs. Get the "X64" version from Microsoft's
   official page (https://learn.microsoft.com/cpp/windows/latest-supported-vc-redist).
4. **Look at the details of the error**: they are saved in `%USERPROFILE%\symvib_error.log`
   (for example `C:\Users\<user name>\symvib_error.log`). If the problem persists, please report it
   on [Issues](https://github.com/youskeishii141-yi/symvib/issues) with the content of this file and
   your Windows version.
5. **Use it from Python**: even where the exe cannot be used, you can start the same GUI with
   `python vibration_gui.py` after installing Python, following "Installation (for use from
   Python)" below.

### For developers (building the exe)

In Windows PowerShell:

```powershell
.\build_exe.ps1
```

This creates `.venv`, installs the dependencies, runs the tests, builds with PyInstaller, and
produces `dist\SymVib\` (the folder containing the exe) and `dist\SymVib.zip`
(for distribution). Because pymatgen is bundled, the extracted size is several hundred MB.

On GitHub, pushing a tag with the same name as `__version__` in `symvib/__init__.py`
(e.g. `v2.1.1`) makes GitHub Actions (`.github/workflows/release.yml`) run the tests, build, and
check start-up on Windows, and attach `SymVib.zip` to the release for that tag (the
release text is `.github/release_notes/<tag>.md`). The zip contains `THIRD_PARTY_LICENSES.txt`
(the full license texts of Python itself and the bundled libraries) generated by
`tools/third_party_notices.py`, `LICENSE`, and the READMEs.

## Installation (for use from Python)

Python 3.9 or later is recommended.

```bash
git clone https://github.com/youskeishii141-yi/symvib.git
cd symvib
pip install -r requirements.txt
```

Dependencies: `numpy`, `pymatgen` (structure reading, point-group symbols), `spglib`, `spgrep`
(crystal symmetry, irreducible representations), `matplotlib` (3D view in the GUI).

## Usage (CLI)

```bash
# GUI
python vibration_gui.py

# Molecular vibration analysis: water demo without arguments
python molecule_vibration_analysis.py
python molecule_vibration_analysis.py mymolecule.xyz --outdir out --tol 0.05

# Crystal factor-group analysis: graphite demo without arguments
python factor_group_analysis.py mycrystal.cif --symprec 1e-3

# Output (report) is in English by default. For Japanese, use --lang ja
python factor_group_analysis.py mycrystal.cif --lang ja

# Formats to write: both (default) / xsf / vesta
python factor_group_analysis.py mycrystal.cif --format vesta
```

The `.vesta` and `.xsf` files for each mode are written to `modes_<name>/`. A `.vesta` file opens
directly with arrows. Open a `.xsf` file in VESTA and display the arrows with Properties →
Vectors → Show vectors. Example output for graphite is in
[examples/modes_graphite/](examples/modes_graphite/). For VESTA, see the next section.

To use it directly from Python:

```python
from symvib.molecule import analyze_molecule, load_molecule
result = analyze_molecule(*load_molecule("benzene.xyz"))
print(result.decomposition)          # Γvib = 2A1g + A2g + ...
for mode in result.modes:
    print(mode.label, mode.activity, mode.displacements[0])   # (N, 3) displacements
```

## About VESTA (crystal-structure visualization software)

The vibrational-mode files written by this tool (`.vesta`, `.xsf`) are meant to be opened in the
crystal-structure visualization software **VESTA** (Koichi Momma and Fujio Izumi) to display the
atomic displacements as arrows in 3D. If the 3D view of the GUI version is sufficient, VESTA is
not needed.

- Obtaining it: it can be downloaded free of charge from the official website
  <https://jp-minerals.org/vesta/> (Windows / macOS / Linux). See the official website for the
  terms of use. VESTA is not bundled with this tool. This tool is also independent and
  unaffiliated with the developers of VESTA.
- `.vesta`: when opened directly, it is displayed with the gray arrows (the displacement of each
  atom) and the bond settings applied.
- `.xsf`: after opening, display the arrows with Properties → Vectors → Show vectors.
- If you use figures made with VESTA in papers or elsewhere, you are asked to cite the VESTA
  reference (K. Momma and F. Izumi, *J. Appl. Crystallogr.* **44**, 1272–1276 (2011)).
- The bond determination in the GUI's 3D view uses VESTA's default bond-length table (the one
  bundled with pymatgen).

## Reading CIF files

The following cases common in real-world CIF files are supported (how they are handled on reading
is noted in the report).

- A CIF containing multiple structures → the first data block with atomic coordinates is used
- Labels that are not element symbols (`Ti1`, `O2a`), with oxidation states (`O-2`, `Si+4`), and
  numbers with uncertainties (`4.583(2)`)
- Partial occupancy: multiple sites occupying the same coordinates are merged into one virtual
  (pseudo) atom whose name is the concatenation of the labels (e.g. `Ti1` and `V2` → `Ti1V2`; the
  mass is averaged by occupancy). Partial occupancy of a single site is treated as full occupancy.
  Duplicate entries of the same element (atoms written with different labels at positions that
  coincide under symmetry operations) are merged into one atom
- Nonexistent element symbols (`Wat`, `X`, etc.) → treated as virtual (pseudo) atoms with a mass
  of 1 amu
- Fractional coordinates (`1/2`), files with a BOM
- No symmetry operations listed → supplied from the Hall symbol, the H-M symbol (`P 21/n`,
  `F d -3 m :2`, ICSD's `Z`/`S`, `R` for rhombohedral axes), or the number, in that order. When the
  choice of origin, etc. is not determined by the symbol, the setting in which atoms do not overlap
  too much (are at least 0.5 Å apart) and the number of atoms is smallest is chosen (P1 if there is
  none)
- The input cell may be conventional, primitive (e.g. rhombohedral axes), or a supercell. Irrep
  labels and Raman-tensor components are expressed in the standard setting of the space group (the
  same orientation as Bilbao)

Structures containing virtual (pseudo) atoms cannot be used for the QE frozen phonon (input an
ordered structure).

## Choice of axes (C2v / D2 / D2h molecules)

The labels B1/B2 (C2v) and B1/B2/B3 (D2, D2h) depend on the choice of coordinate axes. Following
Mulliken's recommendation, this tool takes the normal to the molecular plane as the x axis for
planar molecules (for C2v, z = the C2 axis; for D2/D2h, z = the C2 axis passing through the most
atoms). With this choice, water is 2A₁ + B₂ and formaldehyde is 3A₁ + B₁ + 2B₂. The adopted axes
are shown in the report.

## Validation

Regression tests are in `tests/` (`python -m pytest tests` after
`pip install -r requirements-dev.txt`). In addition to agreement with literature values,
independent checks that do not rely on the code's internal paths (`tests/invariants.py`) confirm
for every mode that "it is an irreducible-representation subspace closed under the symmetry
operations", "it contains no translations or rotations", "the activity and Raman-tensor components
agree with the projection operators", and "the Mulliken labels agree with the point-group tables
of the Bilbao Crystallographic Server".

- Molecules: water 2A₁ + B₂, formaldehyde 3A₁ + B₁ + 2B₂, ethylene, methane, CO₂, HCN, C₂H₂,
  benzene, SF₆, B₁₂H₁₂²⁻, C₆₀ (2Ag + 3T1g + 4T2g + 6Gg + 8Hg + Au + 4T1u + 5T2u + 6Gu + 7Hu), PCl₅,
  ferrocene (D5h / D5d), ethane, allene, B(OH)₃, and others. The character tables of all 45 point
  groups are verified group-theoretically (`tests/test_tables.py`).
- Crystals: the numbers of Raman / IR / inactive modes of 33 materials (7 crystal systems) are
  checked against literature values (`tests/test_literature.py`; the source is given for each
  entry). Synthetic structures for all 230 space groups and 530 settings are covered exhaustively
  (`tests/sweep_space_groups.py`; results in
  [docs/validation_space_groups.md](docs/validation_space_groups.md)). One real CIF from COD is
  analyzed for each space group (`tests/check_cod.py`,
  [docs/validation_cod.md](docs/validation_cod.md)). Also agrees with the irreducible
  representations of phonopy.
- QE frozen phonon: synthetic pw.x outputs are passed through generate → analyze and agree with
  full diagonalization.

External reference data (Bilbao point-group character tables, COD CIF files) are not included in
the distribution. The tests pass without fetching them; only the checks that use them are skipped.
To reproduce these checks as well, fetch the data from each official website and then run the tests
(for COD, fetch the CIF files with the same IDs as recorded in `tests/data/cod/index.json`). The
same can be run manually with the "Full validation" workflow of GitHub Actions
(`.github/workflows/full-validation.yml`).

```bash
python tests/fetch_bilbao_tables.py
python tests/fetch_cod.py
python -m pytest tests
python tests/check_cod.py          # all 227 COD entries (takes a long time)
python tests/sweep_space_groups.py # exhaustive coverage of 230 space groups and 530 settings
```

For details of the checks and a list of the errors found, see
[docs/validation_report.md](docs/validation_report.md).

## Notes and limitations

- This tool **does not compute absolute vibrational frequencies**. It is meant for determining the
  symmetry, degeneracy, activity, and displacement patterns of modes from symmetry.
- The spring-network model (central-force springs) has no bending stiffness, so for out-of-plane modes of
  planar molecules and the like, the estimated direction for m>1 becomes arbitrary. In such cases
  this is noted in the report (the identification of the irreducible representations itself is done
  by pure group theory and is unaffected).
- A molecular CIF is treated as one molecule consisting of the atoms in the unit cell as they are.
  Use a file in which the molecule is not cut at the cell boundary. Structures containing partially
  occupied sites cannot be handled.
- If a point group not covered is input, the number of vibrational modes is shown, and the program
  terminates safely without outputting symmetry-species labels.

## Changelog

### v2.2.3
- GUI: in the mode list on the left, the rows were too short and the tops and bottoms of the
  characters were cut off (most visibly on high-DPI screens). The row height now follows the actual
  height of the font, and the text size matches the analysis report below. The columns fit their
  headings and contents, and scrollbars appear only when the list does not fit vertically or
  horizontally. The border between the list and the report can be dragged up and down.

### v2.2.2
- The startup self-check now also covers failures while setting up the window (Tcl/Tk), not only
  failures to load modules (when an interrupted extraction left out `_tcl_data`, only "Can't find a
  usable init.tcl" was shown).
- For PCs on which "Extract All" stops partway, the README now explains how to extract with
  PowerShell and how to run the self-test (`--self-test`). The self-test result is written as UTF-8
  with a BOM so that Windows PowerShell displays it correctly.

### v2.2.1
- When the distributed exe fails to start, it now checks the extracted files for missing ones
  (an interrupted extraction or files quarantined by security software) and shows the missing
  files and how to fix the problem (previously only NumPy's "DLL load failed" error was shown).
- Added the system requirements and a "Troubleshooting" section to the README.
- Stronger checks when building the exe: a self-test that runs the analyses inside the exe, a check
  of the extracted zip, and a check that the diagnosis works when files are missing.

### v2.2.0 (public release)
- Renamed the software to **SymVib** (Python package `symvib`, Windows version `SymVib.exe`;
  formerly vibration-analysis / `vibration_analysis` / `VibrationAnalysis.exe`). The settings and
  error-log files are now `~/.symvib.json` and `~/symvib_error.log`.
- For the public release, set the license to **MIT + Commons Clause** (use, modification and
  redistribution are free; selling the software is prohibited). Stated a disclaimer, contact
  information, and the fact that development was AI-assisted in the README.
- Added English documentation (README.md and each document in docs/; the Japanese versions are
  `*.ja.md`).
- Added an introduction to VESTA and the sources of the dependency libraries and reference data.
- Included the license notices of the bundled libraries (`THIRD_PARTY_LICENSES.txt`), `LICENSE`, and
  the READMEs in the distribution zip.
- Removed the external data used for testing (Bilbao character tables, COD CIF files) from the
  distribution and made them reproducible with fetch scripts (if not fetched, the corresponding
  checks are skipped).

### v2.1.1
- Bug fixes (mathematical audit before release):
  - Raman tensor form: the factor 1/√2 was missing when converting the shear component
    (xy+yx)/√2 of the orthonormal basis into matrix elements, so the coefficient ratios were wrong
    for forms in which diagonal and off-diagonal components are combined into one letter (such as E
    of trigonal crystals).
  - QE frozen phonon: acoustic/optical classification is now based on the translational component
    in mass-weighted coordinates. Optical modes were sometimes misclassified as acoustic modes in
    crystals with large mass differences. Degenerate modes return all components, and all components
    are also saved to XSF. Invalid displacement amplitudes (zero, negative, NaN) are rejected before
    input generation.
  - Molecules: when the axis of a linear molecule was not along a coordinate axis, the rotational
    subspace could be misidentified.
  - Spring-network model: bonds exactly at the cutoff boundary were only partly counted due to rounding
    errors, so degenerate modes could fail to close under the symmetry operations. The projected
    model matrix is group-averaged to preserve symmetry.
  - When the atom mapping is not bijective (atoms overlap), an error is raised instead of silently
    constructing a singular representation.
  - GUI: the location of the error log was not shown in the error dialog.
- Added the note "What a quantum-chemical calculation determines, and what remains" to the README.
- Tests now run automatically on GitHub Actions on Windows / Linux, and pushing a tag automatically
  builds the distribution exe and attaches it to the release.

### v2.1 (pre-release checks)
- Added checks against external references: the 32 point-group tables of Bilbao, an exhaustive sweep
  of 230 space groups (530 settings), 227 real CIF files from COD, literature values for 33
  materials, and comparison with phonopy (see "Validation" above).
- Bug fixes:
  - Raman tensor form: for the E-type modes of 4/m, -4, 4, 6/m, 6, -6, 3, -3 there are two
    independent parameters (four for 3, -3), but a one-parameter form (the form for 4/mmm, etc.) was
    shown. Instead of looking up the form in a table, it is now computed exactly from the
    displacements of each mode with projection operators.
  - QE frozen phonon: the displacement δ was multiplied by the 3N-normalized basis vector, so the
    displacement per atom shrank by 1/√(number of moving atoms) (forces were buried in noise in
    large cells). Changed to maximum atomic displacement = δ, and the amplitude used is recorded in
    manifest.json.
  - Molecules: the symmetrization of distorted inputs did not converge to the exact group, and a
    residual (~1e-7 Å) remained in the displacement patterns.
  - CIF: for monoclinic structures without listed symmetry operations, with a short symbol
    ("P 21/c") the structure was expanded in the b-axis setting without looking at the lattice angles
    (unique axis).
  - In crystal labelling, the principal axis of D2d / S4 was taken as "the first C2 axis found",
    which could depend on the ordering of the operations → the rotoreflection axis is now taken as the
    principal axis.
  - Crystals: in CIF files whose lattice constants do not match the detected symmetry exactly (e.g.
    the rhombohedral angles 51.683/51.683/51.685° of AsI₃), the rotation matrices deviated from
    orthogonality, the force constants of the spring-network model no longer commuted with the symmetry
    operations, and the estimated orientation of degenerate components (m>1) mixed in different copies
    of the same irreducible representation → the lattice and coordinates are now idealized to the
    detected symmetry before analysis (noted in the report when the change is large).
  - Crystals: the spring-network model counted only ±1 images of neighboring cells, so neighbors within the
    cutoff could be missed in body-centered, face-centered, and rhombohedral primitive cells (COD
    1010030 I-42m) → the range of images is now determined from the cutoff.
  - Molecules: when an atom near the center of mass was chosen as the reference atom for the
    symmetry-operation search, the symmetry could be lost with a slight distortion → it is now chosen
    from the outer atoms.
  - GUI: for crystals consisting only of element pairs not in VESTA's bond table (AsI₃, etc.), when
    there was not a single bond, an exception occurred just before the 3D view and the results were
    not shown.
- Display: added a note that the Raman-tensor components of molecules are expressed in the standard
  orientation of the point group.

### v2.0 (GUI version)
- Added the GUI (`vibration_gui.py`) and the build of the Windows distribution exe (`build_exe.ps1`).
- Reorganized the analysis core into the package `vibration_analysis/`, separating computation and
  display (the CLI commands can be used as before). Added regression tests.
- Bug fixes:
  - On Japanese Windows, XSF output for linear molecules (CO₂, etc.) stopped with a character-encoding
    error. Display errors when the output was piped/redirected were also fixed.
  - In D₂h / D₂ / C₂ᵥ, C2(z/y/x) and σv(xz/yz) could not be distinguished, giving wrong decompositions
    for ethylene and others (e.g. "B1u with degeneracy 3"). In C₂ᵥ, B₂ never appeared.
  - For Th molecules, symmetry operations were missed (6 of 24), giving decompositions with an
    inconsistent number of degrees of freedom. For C₃h / C₅h / S6 / S8 as well, the classification of
    reflections and rotoreflections and their correspondence to the character table were wrong.
  - In slightly distorted structures (experimental or optimized structures), the correspondence
    between the point group and the atoms was off, giving wrong results → changed to symmetrize
    automatically.
  - In crystals, when the same irreducible representation occurred more than once, the individual
    copies were not closed under the symmetry operations, and the acoustic/optical separation and the
    estimated directions could be mixed. The treatment of masses in the spring-network-model estimate was not
    consistent with the molecular version (corrected to a generalized eigenvalue problem).
  - QE input generation stopped for CIF files with oxidation states (e.g. Ti4+).
  - In reading forces from QE output, the breakdown blocks of verbosity='high' could be mistaken for
    the total forces. Missing output files now also give a clear error.
  - The choice of degenerate components, and orientations not determined by the spring-network model (such as
    out-of-plane modes of planar molecules), changed with the order of atoms in the input file →
    changed to common rules determined by geometry alone.
  - Molecules: rebuilt the determination of the classes (conjugacy classes) of symmetry operations so
    that they are obtained from the group's multiplication table and matched to the character table.
    Previously the classes were determined by geometric guessing, which gave wrong decompositions for
    the cubic groups (Td, Oh), D5, D6d, D7h, D8h, Cs in some orientations, and D2d depending on the order
    of atoms (e.g. Mo(CO)₆, Ni(CO)₄, adamantane). Also fixed an error in the D8h character table (the
    S8 and S8³ columns). The consistency of the assignment and the character table is checked with the
    regular representation, and the program stops instead of giving a wrong result if they are
    inconsistent.
  - Molecules: nearly linear molecules (HNCS, etc.) and elongated molecules were misidentified as
    linear molecules and their atoms were moved onto the axis → now determined by whether the atoms lie
    within the tolerance of the molecular axis. Deuterium (D) in xyz files was treated as H. The way
    the C2v/D2/D2h axes were chosen could change with the orientation of the molecule.
  - Crystals: optical modes did not satisfy the condition "the center of mass does not move"
    (Σ M u = 0) (e.g. in T1u of NaCl, Na and Cl moved by the same amount).
  - Crystals: the program stopped for primitive-cell CIF files and CIF files with rhombohedral axes.
    B2/B3 of orthorhombic crystals could be swapped due to numerical errors. B1/B2 of mm2 and 6mm were
    the reverse of the Bilbao convention (wurtzite ZnO is 2B1). The Raman-tensor forms of Ag of
    orthorhombic crystals and B of tetragonal crystals (4/m, etc.) were wrong.
  - CIF: when no symmetry operations were listed, the structure could be expanded into a wrong one
    depending on how the H-M symbol was written (`P 21/n`, etc.) or on the choice of origin. `Wat` was
    read as tungsten and `NH4` as nihonium.
  - GUI: starting another analysis during an analysis could swap the results or the names of the save
    destinations. Save failures and unexpected errors were not shown in the exe (now shown in a dialog
    and the log `~/vibration_analysis_error.log`). The program did not start if the settings file was
    corrupted. The display was slow for large cells. The axis indicator at the lower left kept its old
    orientation right after loading.

## License

This tool is released under the **MIT License with the Commons Clause condition (no selling)**.
The authoritative text is the English [LICENSE](LICENSE). The following is a summary to aid
understanding; where it differs from the license text, the license text prevails.

- **What you may do**: use, modify and redistribute it free of charge for any purpose. Besides
  education, research and personal study, this includes **use by companies for for-profit research
  and development or internal work**.
- **What you may not do**: **sell** the software (including modified versions), that is, provide to
  third parties, for a fee or other consideration, a product or service whose value derives entirely
  or substantially from the functionality of the software. This includes paid hosting, consulting or
  support services related to the software.
- **Conditions for redistribution and distributing modified versions**: include the copyright notice
  in `LICENSE` (Yosuke Ishii, Nagoya Institute of Technology), the MIT license text and the Commons
  Clause condition.
- Citation in papers and notifying the author of your use are **not obligatory**. However, if the tool
  was helpful, a citation or acknowledgment would be encouraging (see "[Citation](#citation)" below).
- If you are unsure whether your use counts as selling (for example, when the software is part of
  your own product), please contact the author.

The license governs the use of this software's code. Copyright protects the expression of the code,
not group-theoretical methods or ideas, so no license can restrict an independent implementation of
the same methods that does not use this code.

Because selling is prohibited, this is not an "open source" license as defined by the OSI, but a
source-available license. Dependency libraries (numpy, pymatgen, spglib, spgrep, matplotlib, etc.)
are subject to their respective licenses. Those bundled in the distribution zip are collected in
`THIRD_PARTY_LICENSES.txt`.

## Disclaimer

This tool is provided free of charge "as is", and the author makes no warranty whatsoever as to the
correctness or completeness of the computed results or behavior, or fitness for a particular
purpose. To the extent permitted by law, the author accepts no liability for errors or defects in
this tool, or for any damage arising from its use. When using the results in research or education,
please verify them against literature values or other means of calculation. The author is also under
no obligation to fix bugs, maintain the tool, or respond to inquiries.

## Contact and feature requests

Bug reports and requests or suggestions for improvement are accepted via
[GitHub Issues](https://github.com/youskeishii141-yi/symvib/issues) (recommended; shared
with other users as well) or by e-mail. However, not every request can be addressed, and replies may
be considerably delayed or may not be possible. Thank you for your understanding.

- Yosuke Ishii (石井 陽祐), Nagoya Institute of Technology
- E-mail: ishii.yosuke [at] nitech.ac.jp (please replace [at] with @)

## Citation

Citation is optional. If you would like to cite this tool, please use the DOI below, which always
resolves to the latest version (each version also has its own DOI on Zenodo):

> Yosuke Ishii. SymVib [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.24085930

APA / BibTeX formats are also available from the repository's "Cite this repository" (generated by
GitHub from [CITATION.cff](CITATION.cff)) and from the [Zenodo record](https://doi.org/10.5281/zenodo.24085930).

## Development (use of AI)

The planning, the conception of features, and the validation strategy of this tool are the author's.
Most of the code writing and checking was done with AI coding assistance (Anthropic's Claude Code and
OpenAI's Codex), and the author has reviewed the content. The theoretical background (group theory,
the great orthogonality theorem, point-group character tables) is based on standard textbooks and
publicly available character tables.

## References: software and data used

This tool uses the following software.

- spglib (crystal symmetry): A. Togo, K. Shinohara, I. Tanaka, *Sci. Technol. Adv. Mater.: Methods* **4**,
  2384822 (2024).
- spgrep (irreducible representations of space groups): K. Shinohara, A. Togo, I. Tanaka, *J. Open Source Softw.* **8**, 5269 (2023).
- pymatgen (structure reading, point-group symbols, VESTA's default bond-length table): S. P. Ong *et al.*,
  *Comput. Mater. Sci.* **68**, 314–319 (2013).
- NumPy: C. R. Harris *et al.*, *Nature* **585**, 357–362 (2020).
- Matplotlib (3D view in the GUI): J. D. Hunter, *Comput. Sci. Eng.* **9**, 90–95 (2007).
- VESTA (viewing the output files; not bundled): K. Momma, F. Izumi, *J. Appl. Crystallogr.* **44**,
  1272–1276 (2011).
- Quantum ESPRESSO (interfaced for frozen phonons; not bundled): P. Giannozzi *et al.*,
  *J. Phys.: Condens. Matter* **21**, 395502 (2009).

The following external references were used for validation (the data are not included in the
distribution).

- Bilbao Crystallographic Server (point-group character tables): M. I. Aroyo *et al.*, *Z. Kristallogr.* **221**, 15–27 (2006).
- Crystallography Open Database (real CIF files): S. Gražulis *et al.*, *J. Appl. Crystallogr.* **42**,
  726–729 (2009).
- phonopy (comparison of irreducible representations): A. Togo, L. Chaput, T. Tadano, I. Tanaka, *J. Phys.: Condens. Matter* **35**,
  353001 (2023).
