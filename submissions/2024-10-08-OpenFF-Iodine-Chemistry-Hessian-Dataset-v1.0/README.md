# OpenFF Iodine Chemistry Hessian Dataset v1.0 

## Description

Hessian single points for the final molecules in the [OpenFF Iodine Chemistry Optimization Dataset v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2022-07-27-OpenFF-iodine-optimization-set) at the B3LYP-D3BJ/DZVP level of theory.

## General information

* Date: 2024-10-07
* Class: OpenFF Basic Dataset
* Purpose: Hessian dataset generation for MSM parameters
* Name: OpenFF Iodine Chemistry Hessian Dataset v1.0 
* Number of unique molecules: 
* Number of conformers: 
* Number of conformers (min, mean, max): 
* Molecular weight (min, mean, max): 
* Charges: 
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

* elements: 
* unique molecules: 
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
