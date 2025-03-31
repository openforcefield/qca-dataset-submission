# OpenFF Cresset Additional Coverage Hessian v4.0

## Description

Hessian single points for the final molecules in the [OpenFF Cresset Additional Coverage Optimizations v4.0 dataset](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2025-03-06-OpenFF-Cresset-Additional-Coverage-Optimizations-v4.0) at the B3LYP-D3BJ/DZVP level of theory.

## General information

* Date: 2024-11-12
* Class: OpenFF Basic Dataset
* Purpose: Hessian dataset generation for MSM parameters
* Name: OpenFF Cresset Additional Coverage Hessian v4.0
* Number of unique molecules: 70
* Number of conformers: 393
* Number of conformers (min, mean, max): 1.00, 5.61, 10.00
* Molecular weight (min, mean, max): 58.08, 144.98, 280.75
* Charges: [-1.0, 0.0, 1.0]
* Dataset generator: Matthew Habgood
* Dataset submitter: Lily Wang


## QCSubmit generation pipeline

* `generate-dataset.ipynb`: This notebook shows how the dataset was prepared from the input files.


## QCSubmit Manifest

* `dataset.json.bz2`: Compressed dataset ready for submission
* `dataset.pdf`: Visualization of dataaset molecules
* `dataset.smi`: Smiles strings for dataset molecules
* `generate-dataset.ipynb`: Notebook describing dataset generation and submission
* `input-env.yaml`: Environment file used to create Python environment for the notebook
* `full-env.yaml`: Fully-resolved environment used to execute the notebook.


## Metadata

* elements: {O, C, F, S, H, N, Br, Cl}
* unique molecules: 70
* Spec: B3LYP-D3BJ/DZVP
    * SCF properties:
        * dipole
        * quadrupole
        * wiberg_lowdin_indices
        * mayer_indices
    * QC Spec:
        * name: default
        * method: B3LYP-D3BJ
        * basis: DZVP
