# OpenFF Protein Dipeptide 2-D TorsionDrive v1.0

## Description

Two-dimensional TorsionDrives on phi and psi for dipeptides of alanine and two rotamers of proline and tryptophan

## General Information

- Date: 2021-11-18
- Class: OpenFF TorsionDrive Dataset
- Purpose: 2-D scans of (phi, psi) for dipeptides
- Collection: TorsiondriveDataset
- Name: OpenFF Protein Dipeptide 2-D TorsionDrive v1.0
- Number of unique molecules: 3
- Number of filtered molecules: 0
- Number of torsion drives: 5
- Number of conformers min mean max: 408 508.80 576
- Dataset Submitter: Chapin Cavender
- Dataset Generator: Chapin Cavender
- Set of charges: [0.0]
- Mean molecular weight: 200.64
- Max molecular weight: 259.30
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
- `dataset.json.bz2`: Compressed dataset ready for submission
- `dataset.pdf`: Visualization of 2D graphs for dataset molecules with driven torsions highlighted
- `dataset.smi`: Smiles strings for dataset molecules
- `dipeptide_rotamers` contains SDF files with conformers scanning phi and psi for dipeptides of alanine and two rotamers of proline and tryptophan both before and after minimization
- `openff_protein_dipeptide_2d_torsiondrive.ipynb`: Notebook describing dataset generation and submission

## Metadata

- elements: {'H', 'N', 'C', 'O'}
- unique molecules: 3
- scf properties:
    - dipole
    - quadrupole
    - wiberg_lowdin_indices
    - mayer_indices
- qc_spec
    - name: default
    - method: B3LYP-D3BJ
    - basis: DZVP

