# OpenFF Protein Capped 1-mer Sidechains v1.0

## Description

TorsionDrives on chi1 and chi2 for capped 1-mers of amino acids with a rotatable bond in the sidechain (not ALA, GLY, or PRO)

## Changelog

- v1.0 Initial submission
- v1.1 Screened initial conformers with bad geometry using qcelemental `guess_connectivity`
- v1.2 Additional screen of initial conformers based on tetrahedral geometry for tetravalent atoms

## General Information

- Date: 2022-02-10
- Class: OpenFF TorsionDrive Dataset
- Purpose: 2-D scans of (chi1, chi2) for protein capped 1-mers
- Collection: TorsiondriveDataset
- Name: OpenFF Protein Capped 1-mer Sidechains v1.2
- Number of unique molecules: 23
- Number of filtered molecules: 0
- Number of torsion drives: 46
- Number of conformers min mean max: 24 454.91 576
- Dataset Submitter: Chapin Cavender
- Dataset Generator: Chapin Cavender
- Set of charges: [-1.0, 0.0, 1.0]
- Mean molecular weight: 206.90
- Max molecular weight: 350.46
- Enumerate stereoisomers: False
- Enumerate tautomers: False
- Enumerate protomers: False

## QCSubmit Generation procedure

- `openff_protein_capped_1-mer_sidechains.ipynb`: This notebook shows how the dataset was prepared from the input files
    - Initial conformations of capped 1-mers were taken from [https://github.com/openforcefield/amber-ff-porting/blob/master/AllDipeptides.tar.gz]
    - Sidechain dihedrals scanned using the OpenEye Toolkit
    - Conformers were minimized with the OpenFF 2.0.0 (Sage) force field using OpenMM
    - Minimized conformers were screened by changes in connectivity and in tetrahedral geometry around tetravalent atoms

## QCSubmit Manifest

- `AllDipeptides.tar.gz` contains MOL2 files for 26 main chain and 23 N-terminal and C-terminal capped 1-mers from [https://github.com/openforcefield/amber-ff-porting/blob/master/AllDipeptides.tar.gz]
- `dataset.json.bz2`: Compressed dataset ready for submission
- `dataset.pdf`: Visualization of 2D graphs for dataset molecules with driven torsions highlighted
- `dataset.smi`: Smiles strings for dataset molecules
- `backbone_conformations` contains SDF files with conformers scanning chi1 and chi2 for capped 1-mers with two backbone conformations both before and after minimization
- `openff_protein_capped_1-mer_sidechains.ipynb`: Notebook describing dataset generation and submission

## Metadata

- elements: {'H', 'N', 'C', 'O', 'S'}
- unique molecules: 23
- scf properties:
    - dipole
    - quadrupole
    - wiberg_lowdin_indices
    - mayer_indices
- qc_spec
    - name: default
    - method: B3LYP-D3BJ
    - basis: DZVP

