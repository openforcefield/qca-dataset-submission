# OpenFF Optimization Hessians 2019-07 to 2025-03 v4.0

## Description

Hessian single points for the final molecules in the [OpenFF Cresset Additional Coverage Optimizations v4.0 dataset](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2025-03-06-OpenFF-Cresset-Additional-Coverage-Optimizations-v4.0) at the B3LYP-D3BJ/DZVP level of theory. These are used for calculating MSM starting points in force field fits. The molecules here include the F, N, H, O, Cl, S, Br, C elements and the charge states -1, 0, +1. They range from 58-281 Da (mean 145) and 4-19 heavy atoms.

## General information

* Date: 2025-04-14
* Class: OpenFF Basic Dataset
* Purpose: Hessian dataset generation for MSM parameters, from OpenFF datasets from 2019-07 to 2025-03
* Name: OpenFF Optimization Hessians 2019-07 to 2025-03 v4.0
* Number of unique molecules: 63446
* Number of conformers: 297934
* Number of conformers (min, mean, max): 1, 4.83, 275
* Molecular weight (min, mean, max): 16.04, 224.36, 1425.34
* Charges: [-1.0, 0.0, 1.0]
* Dataset generator: Lily Wang
* Dataset submitter: Lily Wang


## QCSubmit generation pipeline

* `generate-dataset.ipynb`: This notebook shows how the dataset was prepared from the input files.


## QCSubmit Manifest

* `optimizations.json.tar.gz`: Compressed combined input OptimizationResultCollection
* `dataset.json.bz2`: Compressed dataset ready for submission
* `dataset.pdf.tar.gz`: Compressed visualization of dataset molecules
* `dataset.smi`: Smiles strings for dataset molecules
* `generate-dataset.ipynb`: Notebook describing dataset generation and submission
* `input-env.yaml`: Environment file used to create Python environment for the notebook
* `full-env.yaml`: Fully-resolved environment used to execute the notebook.


## Metadata

* elements: {S, H, O, Br, F, N, P, Cl, I, C}
* unique molecules: 63446
* Spec: default
  * basis: DZVP
  * implicit_solvent: None
  * keywords: {}
  * maxiter: 200
  * method: B3LYP-D3BJ
  * program: psi4
  * SCF properties:
    * dipole
    * quadrupole
    * wiberg_lowdin_indices
    * mayer_indices

