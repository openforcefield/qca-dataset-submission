# OpenFF Protein Capped 1-mer Sidechains v1.0

## Description

TorsionDrives on phi and psi for capped 3-mers of Ace-Y-X-Y-Nme with Y = {Ala, Val} and X = 27 canonical amino acids (Ash, Cyx, Glh, Hid, Hip, and Lyn and Pro with omega both cis and trans)

## Changelog

- v1.0 Initial submission

## General Information

- Date: 2022-05-30
- Class: OpenFF TorsionDrive Dataset
- Purpose: 2-D scans of (phi, psi) for protein capped 3-mers
- Collection: TorsiondriveDataset
- Name: OpenFF Protein Capped 3-mer Backbones v1.0
- Number of unique molecules: 52
- Number of filtered molecules: 0
- Number of torsion drives: 54
- Number of conformers min mean max: 389 561.94 576
- Dataset Submitter: Chapin Cavender
- Dataset Generator: Chapin Cavender
- Set of charges: [-1.0, 0.0, 1.0]
- Mean molecular weight: 369.22
- Max molecular weight: 548.72
- Enumerate stereoisomers: False
- Enumerate tautomers: False
- Enumerate protomers: False

## QCSubmit Generation procedure

- `openff_protein_capped_3-mer_backbones.ipynb`: This notebook shows how the dataset was prepared from the input files
    - Initial conformations were generated from SMILES using the OpenEye Toolkit
    - Sidechain dihedrals scanned using the OpenEye Toolkit
    - Conformers were minimized with the OpenFF 2.0.0 (Sage) force field using OpenMM
    - Minimized conformers were screened by changes in connectivity and in tetrahedral geometry around tetravalent atoms

## QCSubmit Manifest

- `dataset.json.bz2`: Compressed dataset ready for submission
- `dataset.pdf`: Visualization of 2D graphs for dataset molecules with driven torsions highlighted
- `dataset.smi`: Smiles strings for dataset molecules
- `capped_3mer_conformations` contains SDF files with conformers scanning phi and psi for capped 3-mers both before and after minimization
- `openff_protein_capped_3-mer_backbones.ipynb`: Notebook describing dataset generation and submission

## Metadata

- elements: {'H', 'C', 'N', 'O', 'S'}
- unique molecules: 52
- scf properties:
    - dipole
    - quadrupole
    - wiberg_lowdin_indices
    - mayer_indices
- qc_spec
    - name: default
    - method: B3LYP-D3BJ
    - basis: DZVP

