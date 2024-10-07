# OpenFF Aniline Para Hessian v1.0

## Description

Hessian single points for the final molecules in the (OpenFF Aniline Para Opt v1.0 dataset)[https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2021-04-02-OpenFF-Aniline-Para-Opt-v1.0] at the B3LYP-D3BJ/DZVP level of theory.

## General information

* Date: 2024-10-07
* Class: OpenFF Basic Dataset
* Purpose: Hessian dataset generation for MSM parameters
* Name: OpenFF Aniline Para Hessian v1.0
* Number of unique molecules: 50
* Number of conformers: 223
* Number of conformers (min, mean, max): 1 4.46 10
* Molecular weight (min, mean, max): 107.15370300000001 150.423601985 343.842162
* Charges: -1.0, 0.0, 1.0
* Dataset submitter: Alexandra McIsaac
* Dataset generator: Alexandra McIsaac


## QCSubmit generation pipeline

* `generate-dataset.ipynb`: This notebook shows how the dataset was prepared from the input files.


## QCSubmit Manifest

* `dataset.json.bz2`: Compressed dataset ready for submission
* `dataset.pdf`: Visualization of dataaset molecules
* `dataset.smi`: Smiles strings for dataset molecules
* `generate-dataset.ipynb`: Notebook describing dataset generation and submission
* `input-environment.yaml`: Environment file used to create Python environment for the notebook
* `full-environment.yaml`: Fully-resolved environment used to execute the notebook.


## Metadata

* elements: {'O', 'Cl', 'S', 'Br', 'H', 'F', 'N', 'C'}
* unique molecules: 50
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
