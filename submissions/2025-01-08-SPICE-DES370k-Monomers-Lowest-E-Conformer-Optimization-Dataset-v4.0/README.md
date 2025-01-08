# SPICE DES370k Monomers Lowest E Conformer Optimization Dataset v4.0

## Description
A dataset containing the lowest energy conformer of all molecules from the [`SPICE DES Monomers Single Points Dataset v1.1`](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2021-11-15-QMDataset-DES-monomers-single-points) dataset, optimized at the OpenFF default level of theory (B3LYP-D3BJ/DZVP). Detailed description on how the original dataset is generated can be found at https://github.com/openmm/qmdataset/tree/main/des370k.

## General information
* Date: 2025-01-08
* Class: OpenFF Optimization Dataset
* Purpose: Conformer optimization
* Name: SPICE DES370k Monomers Lowest E Conformer Optimization Dataset v4.0
* Number of unique molecules: 374
* Number of conformers: 374
* Number of conformers (min, mean, max): 1.00, 1.00, 1.00
* Molecular weight (min, mean, max): 16.04, 95.89, 284.78
* Charges: -1.0 0.0 1.0
* Dataset submitter: Alexandra McIsaac
* Dataset generator: Alexandra McIsaac

## QCSubmit generation pipeline
* `des370k_monomers_minEconf.json`: Dataset containing the minimum energy conformers to use as a starting point for the optimization. Needed as an input file to `generate-dataset.ipynb`
* `generate-dataset.ipynb`: Notebook used to generate dataset

## QCSubmit Manifest
* `dataset.json.bz2`: Compressed dataset ready for submission
* `dataset.pdf`: Visualization of dataset molecules
* `dataset.smi`: Smiles strings for dataset molecules
* `des370k_monomers_minEconf.json`: Dataset containing the minimum energy conformers to use as a starting point for the optimization. Needed as an input file to `generate-dataset.ipynb`
* `generate-dataset.ipynb`: Notebook used to generate dataset
* `input_environment.yaml`: Environment file used to create Python environment for the notebook
* `input_environment_full.yaml`: Fully-resolved environment used to execute the notebook.

## Metadata
* Elements: {Br, S, O, N, H, Cl, I, P, C, F}
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