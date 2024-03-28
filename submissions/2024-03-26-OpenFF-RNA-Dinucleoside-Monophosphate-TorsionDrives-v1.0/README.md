# OpenFF RNA Dinucleoside Monophosphate TorsionDrives v1.0

## Description

TorsionDrives for RNA dinucleoside monophosphates (DNMP), i.e. XpY 2-mers,
on epsilon, zeta, alpha, beta, gamma, chi, 2' OH, 3' OH, and 5' OH
with non-driven dihedrals constrained to A-form helix references values.

### Change log
* v1.0 Initial submission

## General information

* Date: 2024-03-26
* Class: OpenFF TorsionDrive Dataset
* Purpose: Torsion scans for RNA backbone, glycosidic, and hydroxyl dihedrals
* Name: OpenFF RNA Dinucleoside Monophosphate TorsionDrives v1.0
* Number of unique molecules: 16
* Number of conformers (min, mean, max): 24 24.0 24
* Molecular weight (min, mean, max): 547.39 579.91 627.44
* Charges: -1
* Dataset submitter: Chapin Cavender
* Dataset generator: Chapin Cavender


## QCSubmit generation pipeline

* `generate-dataset.ipynb`: This notebook shows how the dataset was prepared from the input files.

## QCSubmit Manifest

* `dataset.json.bz2`: Compressed dataset ready for submission
* `dataset.pdf`: Visualization of 2D graphs for dataset molecules with driven torsions highlighted
* `dataset.smi`: Smiles strings for dataset molecules
* `generate-dataset.ipynb`: Notebook describing dataset generation and submission
* `initial-conformers.tar.gz`: SDF files containing reference conformers with A-form helix dihedrals and initial conformers used to seed TorsionDrives
* `input-environment.yaml`: Environment file used to create Python environment for the notebook
* `full-environment.yaml`: Fully-resolved environment used to execute the notebook

## Metadata

* elements: {'H', 'C', 'N', 'O', 'P'}
* unique molecules: 16
* Spec: default
    * SCF properties:
        * dipole
        * quadrupole
        * wiberg_lowdin_indices
        * mayer_indices
    * QC Spec:
        * name: default
        * method: B3LYP-D3BJ
        * basis: DZVP
