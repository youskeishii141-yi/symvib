English | [日本語](README_molecule.ja.md)

# Molecular Vibrational Mode Analysis (molecule_vibration_analysis.py)

An educational tool that visualizes the vibrational modes of isolated molecules without any DFT calculation.
It lets you experience, with familiar molecules, the same "symmetry-exact vs. network guess"
philosophy as the crystal version (factor_group_analysis.py). A bridge from molecules to crystals.

## What it does
Given a CIF or xyz file, it:
- Automatically determines the molecular point group, 45 point groups in total (nearly all of the Cn/Cnv/Cnh/Sn/Dn/Dnd/Dnh series for n=2–8, the cubic groups T/Th/Td/O/Oh, the icosahedral groups I/Ih, and the linear-molecule groups C∞v/D∞h)
- Labels the internal vibrational modes by irreducible representation (Mulliken symbols; Σ/Π for linear molecules)
- Gives the Raman / IR activity and the polarizability tensor components of each mode
- Outputs "in which direction each atom moves" as XSF (for VESTA)
- For m=1 the direction is fixed by symmetry; for m>1 it is estimated with a 1/d² spring-network model (with a disclaimer)

It does not compute absolute vibrational frequencies (for that, turn to quantum-chemistry calculations such as Gaussian / ORCA).
This tool focuses on understanding "what modes exist and how they move."

## Usage
```
pip install numpy spglib pymatgen
python molecule_vibration_analysis.py mymolecule.xyz
python molecule_vibration_analysis.py mymolecule.cif
python molecule_vibration_analysis.py            # no arguments → demo with water
```
→ An XSF file for each mode is written to `modes_<name>/`.
Open it in VESTA and display the arrows via Properties → Vectors → Show vectors.

## Rigorous treatment of linear molecules (infinite groups)
Linear molecules such as CO2 and acetylene belong to infinite groups (C∞v, D∞h).
This tool decomposes them rigorously with a "true continuous-group projection operator"
P = (l/h) ∫ χ*(R) D(R), obtained by numerically integrating over the continuous family of rotation angles φ (it takes no shortcut via a finite subgroup).
Theoretical rigor is prioritized so that students reading the code can see the great orthogonality
theorem for infinite groups implemented as it stands.

Validation examples (in agreement with literature values):
- Water H2O (C2v): 2A1 + B2 (taking x along the normal to the molecular plane, as recommended by Mulliken)
- Ammonia NH3 (C3v): 2A1 + 2E
- Methane CH4 (Td): A1 + E + 2T2
- CO2 (D∞h): Σg⁺ + Σu⁺ + Πu
- HCN (C∞v): 2Σ⁺ + Π
- Acetylene C2H2 (D∞h): 2Σg⁺ + Σu⁺ + Πg + Πu
- SF6 (Oh): A1g + Eg + 2T1u + T2g + T2u
- C60 (Ih): 2Ag + 3T1g + 4T2g + 6Gg + 8Hg + Au + 4T1u + 5T2u + 6Gu + 7Hu
  (10 Raman-active, 4 IR-active — consistent with the standard assignment in fullerene spectroscopy)

## Certainty of the directions (same as the crystal version)
- **Symmetry-exact**: the irreducible representation occurs only once (m=1) → group theory determines the direction uniquely (exact)
- **Guess (network)**: the same irreducible representation occurs more than once (m>1) → estimate from the spring-network model
  (the true orientation requires a quantum-chemistry calculation). Example: the two A1 modes of water (symmetric stretch and bend)
  can mix, so they are network guesses, whereas B2 (antisymmetric stretch) is symmetry-exact.

## Accuracy of the irreducible decomposition
- The symmetry operations are searched for by the tool itself (because pymatgen's enumeration of operations misses some for Sn and Th); the conjugacy
  classes are obtained from the group multiplication table and mapped to the columns of the character table. For classes with the same magnitude and rotation angle, such as C2'/C2'' or σv/σd,
  the choice follows the Mulliken convention of designating the axes/planes that pass through atoms as C2' and σv (verified on benzene D6h).
- The consistency between the assignment and the character table is checked by confirming that, on the regular representation, the projection for each irreducible representation is a projection operator (P² = P, tr P = d²);
  if they are inconsistent, the program stops instead of producing an incorrect result.
- For groups with complex irreducible representations (Cn, Cnh, Sn, T, Th), each complex-conjugate pair is treated as a physical two-dimensional representation (E notation);
  the operations are classified by their signed rotation angle about the principal axis, and the multiplicities are obtained with the complex characters.
- For the icosahedral groups (I, Ih), the conjugacy classes are determined by the rotation angle alone, so they are matched by angle.
- The Raman tensor components (xx, xz, …) are expressed in the standard orientation of the point group (z = principal axis).
  Even if the input coordinates have a different orientation, the notation of the components does not change (the displacements are output in the input coordinate system).
- For the character tables of all 45 point groups, the projection-operator properties and the IR/Raman assignments
  (decomposition of the vector representation and of the symmetric rank-2 tensor representation) are verified on the ideal groups in `tests/test_tables.py`.

## Relation to the crystal version
You can experience the same "difference between symmetry-adapted modes and true eigenvectors" in both
molecules and crystals. The intended learning path is to grasp the concepts with molecules first and then
move on to crystals (periodic systems) with factor_group_analysis.py.
