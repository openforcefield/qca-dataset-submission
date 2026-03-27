# OpenFF Iodine Chemistry Hessian Dataset v1.0 

## Description

Hessian single points for the final molecules in the [OpenFF Iodine Chemistry Optimization Dataset v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2022-07-27-OpenFF-iodine-optimization-set) at the B3LYP-D3BJ/DZVP level of theory.

## General information

* Date: 2024-11-11
* Class: OpenFF Basic Dataset
* Purpose: Hessian dataset generation for MSM parameters
* Name: OpenFF Iodine Chemistry Hessian Dataset v1.0 
* Number of unique molecules: 99
* Number of conformers: 327
* Number of conformers (min, mean, max): 1, 3.30, 22
* Molecular weight (min, mean, max): 221.95 318.54 533.92
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

* elements: I, F, Br, C, Cl, O, S, N, H
* unique molecules: 99
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
