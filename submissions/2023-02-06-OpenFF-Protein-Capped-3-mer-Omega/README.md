# OpenFF Protein Capped 3-mer Omega v1.0

## Description

TorsionDrives on omega for capped 3-mers of Ace-Ala-X-Ala-Nme with X = 27 canonical amino acids (Ash, Cyx, Glh, Hid, Hip, and Lyn)

## Changelog

- v1.0 Initial submission

## General Information

- Date: 2023-02-06
- Class: OpenFF TorsionDrive Dataset
- Purpose: Scans of omega for protein capped 3-mers
- Collection: TorsiondriveDataset
- Name: OpenFF Protein Capped 3-mer Omega v1.0
- Number of unique molecules: 26
- Number of filtered molecules: 0
- Number of torsion drives: 26
- Number of conformers min mean max: 24 24.00 24
- Dataset Submitter: Chapin Cavender
- Dataset Generator: Chapin Cavender
- Set of charges: [-1.0, 0.0, 1.0]
- Mean molecular weight: 342.28
- Max molecular weight: 492.61
- Enumerate stereoisomers: False
- Enumerate tautomers: False
- Enumerate protomers: False

## QCSubmit Generation procedure

- `openff_protein_capped_3-mer_omega.ipynb`: This notebook shows how the dataset was prepared from the input files
    - Initial conformations were generated from SMILES using the OpenEye Toolkit
    - Omega dihedral was scanned using the OpenEye Toolkit
    - Conformers were minimized with the OpenFF 2.0.0 (Sage) force field using OpenMM
    - Minimized conformers were screened by changes in connectivity and in tetrahedral geometry around tetravalent atoms

## QCSubmit Manifest

- `dataset.json.bz2`: Compressed dataset ready for submission
- `dataset.pdf`: Visualization of 2D graphs for dataset molecules with driven torsions highlighted
- `dataset.smi`: Smiles strings for dataset molecules
- `capped_3mer_conformations` contains SDF files with conformers scanning omega for capped 3-mers both before and after minimization
- `openff_protein_capped_3-mer_omega.ipynb`: Notebook describing dataset generation and submission

## Metadata

- elements: {'H', 'C', 'N', 'O', 'S'}
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

