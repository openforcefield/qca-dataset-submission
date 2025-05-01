# SPICE Dipeptides Lowest E Conformer Optimization Dataset v4.0

## Description
A dataset containing the lowest energy conformer of all molecules from the `Dipeptides` subset of the SPICE dataset, optimized at the OpenFF default level of theory (B3LYP-D3BJ/DZVP). Detailed description on how the original dataset is generated can be found at https://github.com/openmm/spice-dataset/tree/main/dipeptides.

## General information
* Date: 2025-01-08
* Class: OpenFF Optimization Dataset
* Purpose: Conformer optimization
* Name: SPICE Dipeptides Lowest E Conformer Optimization Dataset v4.0
* Number of unique molecules: 677
* Number of conformers: 677
* Number of conformers (min, mean, max): 1.00, 1.00, 1.00
* Molecular weight (min, mean, max): 187.20, 313.73, 445.52
* Charges: -2.0 -1.0 0.0 1.0 2.0
* Dataset submitter: Alexandra McIsaac
* Dataset generator: Alexandra McIsaac

## QCSubmit generation pipeline
* `dipeptide_minEconf.json`: Dataset containing the minimum energy conformers to use as a starting point for the optimization. Needed as an input file to `generate-dataset.ipynb`
* `generate-dataset.ipynb`: Notebook used to generate dataset

## QCSubmit Manifest
* `dataset.json.bz2`: Compressed dataset ready for submission
* `dataset.pdf`: Visualization of dataset molecules
* `dataset.smi`: Smiles strings for dataset molecules
* `dipeptide_minEconf.json`: Dataset containing the minimum energy conformers to use as a starting point for the optimization. Needed as an input file to `generate-dataset.ipynb`
* `generate-dataset.ipynb`: Notebook used to generate dataset
* `input_environment.yaml`: Environment file used to create Python environment for the notebook
* `input_environment_full.yaml`: Fully-resolved environment used to execute the notebook.

## Metadata
* Elements: {H, S, C, O, N}
* Spec: default-mbis
  * basis: DZVP
  * implicit_solvent: None
  * keywords: {}
  * maxiter: 200
  * method: B3LYP-D3BJ
  * program: psi4
  * SCF properties:
    * dipole
    * quadrupole
    * mbis_charges
    * wiberg_lowdin_indices
    * mayer_indices
