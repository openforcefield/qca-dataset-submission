# OpenFF Protein Dipeptide 2-D TorsionDrive v1.1

## Description

Two-dimensional TorsionDrives on phi and psi for dipeptides of the most populated rotamer of the 20 canonical amino acids and common alternate tautomers or protomers (ASH, CYX, GLH, HID, HIP, LYN)

## Changelog

Note: v1.0 was submitted with `openff-qcsubmit` 0.2.4, but this dataset requires 0.3.0 to make use of the constraints defined in `additional_keywords` for each torsiondrive optimization.

- v1.0 Dipeptides of alanine and two rotamers each of proline and tryptophan with `openff-qcsubmit` 0.2.4
- v1.1 Dipeptides of alanine and two rotamers each of proline and tryptophan with `openff-qcsubmit` 0.3.0
- v2.0 Dipeptides of most populated rotamer of 20 canonical amino acids and ASH, CYX, GLH, HID, HIP, and LYN

## General Information

- Date: 2021-11-18
- Class: OpenFF TorsionDrive Dataset
- Purpose: 2-D scans of (phi, psi) for dipeptides
- Collection: TorsiondriveDataset
- Name: OpenFF Protein Dipeptide 2-D TorsionDrive v2.0
- Number of unique molecules: 26
- Number of filtered molecules: 0
- Number of torsion drives: 26
- Number of conformers min mean max: 408 569.54 576
- Dataset Submitter: Chapin Cavender
- Dataset Generator: Chapin Cavender
- Set of charges: [-1.0, 0.0, 1.0]
- Mean molecular weight: 200.12
- Max molecular weight: 350.46
- Enumerate stereoisomers: False
- Enumerate tautomers: False
- Enumerate protomers: False

## QCSubmit Generation procedure

- `openff_protein_dipeptide_2d_torsiondrive.ipynb`: This notebook shows how the dataset was prepared from the input files
    - Initial dipeptide conformations were taken from [https://github.com/openforcefield/amber-ff-porting/blob/master/AllDipeptides.tar.gz]
    - Dipeptide rotamers were generated from the MolProbity ultimate rotamer library and backbone dihedrals scanned using the OpenEye Toolkit
    - Conformers were minimized with the OpenFF 2.0.0 (Sage) force field using OpenMM

## QCSubmit Manifest

- `AllDipeptides.tar.gz` contains MOL2 files for 26 main chain and 23 N-terminal and C-terminal dipeptides from [https://github.com/openforcefield/amber-ff-porting/blob/master/AllDipeptides.tar.gz]
- `dataset-vX.Y.json.bz2`: Compressed dataset ready for submission
- `dataset.pdf`: Visualization of 2D graphs for dataset molecules with driven torsions highlighted
- `dataset.smi`: Smiles strings for dataset molecules
- `dipeptide_rotamers` contains SDF files with conformers scanning phi and psi for dipeptides both before and after minimization
- `openff_protein_dipeptide_2d_torsiondrive.ipynb`: Notebook describing dataset generation and submission

## Metadata

- elements: {'H', 'N', 'C', 'O', 'S'}
- unique molecules: 26
- scf properties:
    - dipole
    - quadrupole
    - wiberg_lowdin_indices
    - mayer_indices
- qc_spec
    - name: default
    - method: B3LYP-D3BJ
    - basis: DZVP

