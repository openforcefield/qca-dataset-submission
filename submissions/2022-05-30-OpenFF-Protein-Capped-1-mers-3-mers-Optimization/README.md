# OpenFF Protein Capped 1-mers 3-mers Optimization Dataset v1.0

## Description

Optimization dataset for protein capped 1-mers Ace-X-Nme and capped 3-mers Ace-Y-X-Y-Nme with Y = {Ala, Val} and X = 26 canonical amino acids with common protomers/tautomers (Ash, Cyx, Glh, Hid, Hip, and Lyn).

## Changelog

- v1.0 Initial submission

## General Information

- Date: 2022-05-30
- Class: OpenFF Optimization Dataset
- Purpose: Optimization dataset for protein capped 1-mers and capped 3-mers
- Collection: OptimizationDataset
- Name: OpenFF Protein Capped 1-mers 3-mers Optimization Dataset v1.0
- Number of unique molecules: 78
- Number of filtered molecules: 0
- Number of conformers: 759
- Number of conformers min mean max: 389 561.94 576
- Dataset Submitter: Chapin Cavender
- Dataset Generator: Chapin Cavender
- Set of charges: [-1.0, 0.0, 1.0]
- Mean molecular weight: 313.59
- Max molecular weight: 548.72
- Enumerate stereoisomers: False
- Enumerate tautomers: False
- Enumerate protomers: False

## QCSubmit Generation procedure

- `openff_protein_capped_1-mers_3-mers_optimization.ipynb`: This notebook shows how the dataset was prepared from the input files

## QCSubmit Manifest

- `dataset.json.bz2`: Compressed dataset ready for submission
- `dataset.pdf`: Visualization of 2D graphs for dataset molecules
- `capped_1-mers_3-mers.smi`: Smiles strings for dataset molecules
- `openff_protein_capped_1-mers_3-mers_optimization.ipynb`: Notebook describing dataset generation and submission

## Metadata

- elements: {'H', 'C', 'N', 'O', 'S'}
- unique molecules: 78
- scf properties:
    - dipole
    - quadrupole
    - wiberg_lowdin_indices
    - mayer_indices
- qc_spec
    - name: default
    - method: B3LYP-D3BJ
    - basis: DZVP

