# OpenFF CX3-CX4 singlepoints v4.0

## Description

A dataset of single-point calculations generated to train bond and angle parameters for a test case of torsion drives for the t17 and t18 torsions in Sage 2.2.1. Conformers were generated following the process in https://github.com/lilyminium/refit-t17-t18-torsions/tree/main/01_generate-singlepoints .

In short, for a particular molecular graph, a single conformer was generated. This molecule was assigned parameters from the Sage 2.2.1 force field. Additional torsion restraints of 1e5 kJ/mol were applied to restrain every torsion in the molecule. An MD simulation was then conducted at 500K and with 0.1 fs timestep. Frames were grouped by torsion similarity, and conformers were chosen from the biggest cluster.

This process was repeated for every molecule in which a t17 or t18 torsion has been driven. Some additional small molecules were also included. Molecules for which the process only generated a single conformer were excluded from the dataset.

This dataset is computed at the default OpenFF level of theory. The molecules here are all neutral and include the Br, C, Cl, F, H, I, N, O, S elements. They range from 30-307 Da (mean 162) and 2-21 heavy atoms. 

## General information

* Date: 2025-05-21
* Class: OpenFF Basic Dataset
* Purpose: Single point dataset generation for training bond and angle parameters
* Name: OpenFF CX3-CX4 singlepoints v4.0
* Number of unique molecules: 365
* Number of conformers: 2938
* Number of conformers (min, mean, max): 2, 8, 10
* Molecular weight (min, mean, max): 30, 162, 307
* Charges: [0.0]
* Dataset generator: Lily Wang
* Dataset submitter: Lily Wang


## QCSubmit generation pipeline

* `generate-dataset.ipynb`: This notebook shows how the dataset was prepared from the input files.


## QCSubmit Manifest

* `input-structures.tar.gz`: Compressed combined input files
* `dataset.json.bz2`: Compressed dataset ready for submission
* `dataset.pdf`: Visualization of dataset molecules
* `dataset.smi`: Smiles strings for dataset molecules
* `generate-dataset.ipynb`: Notebook describing dataset generation and submission
* `input-env.yaml`: Environment file used to create Python environment for the notebook
* `full-env.yaml`: Fully-resolved environment used to execute the notebook.


## Metadata

* elements: Br, C, Cl, F, H, I, N, O, S
* unique molecules: 365
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
