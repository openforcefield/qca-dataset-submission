# OpenFF Gen2 Hessian Dataset Protomers v1.0

## Description

Hessian single points for the final molecules in the [OpenFF Gen2 Optimization Dataset Protomers v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2021-12-21-OpenFF-Gen2-Optimization-Set-Protomers) at the B3LYP-D3BJ/DZVP level of theory.

## General information

* Date: 2024-10-07
* Class: OpenFF Basic Dataset
* Purpose: Hessian dataset generation for MSM parameters
* Name: OpenFF Gen2 Hessian Dataset Protomers v1.0
* Number of unique molecules: 108
* Number of conformers: 597
* Number of conformers (min, mean, max): 1, 5.53, 10
* Molecular weight (min, mean, max): 82.06, 282.05, 542.59
* Charges: -3.0, -2.0, -1.0, 0.0, 1.0, 2.0
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

* elements: {'H', 'C', 'Cl', 'P', 'F', 'Br', 'O', 'N', 'S'}
* unique molecules: 108
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
