# TorsionNet500-Re-optimization-TorsionDrives-v1.0

### Description

TorsionNet500 re-optimized with OpenFF default spec

## General Information

 - Date: 2026-02-02
 - Class: OpenFF TorsionDrive Dataset
 - Purpose: Generate a dataset for benchmarking torsion profiles
 - Name: 2021-11-09-TorsionNet500-single-points
 - Number of unique molecules:        500
 - Number of filtered molecules:      0
 - Number of conformers:              12000
 - Number of conformers min mean max: 24, 24.00, 24
 - Set of charges: [0.0]
 - Dataset Submitter: Matt Thompson
 - Dataset Source: https://github.com/PfizerRD/TorsionNet/tree/main/data

## Changelog

Here any information regarding dataset changes are recorded.

## QCSubmit generation pipeline

 - `Dataset_Generation.ipynb`: This notebook shows how the dataset was prepared from the input files. 
 
## QCSubmit Manifest

- `Dataset_Generation.ipynb`: Dataset creation notebook with instructions for submission.
- `dataset.json.bz2`: The compressed constrained optimization dataset ready for submission.
- `TorsionNet500.smi`: The smiles file of the peptide molecules.
- `TorsionNet500.pdf`: A pdf file containing molecule 2D structures.
- `TorsionNet500_qm_opt_geometries.sdf.gz`: File which contains structures and is the main input for the dataset
 
## Metadata

* Number of unique molecules: 500
* Number of driven torsions: 500
* Number of filtered molecules: 0
* Number of conformers: 500
* Number of conformers per molecule (min, mean, max): 1, 1.00, 1
* Mean molecular weight: 183.52
* Min molecular weight: 70.13
* Max molecular weight: 268.74
* Charges: [0.0]
## Metadata
* Elements: {N, Cl, O, S, H, C, F}
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
