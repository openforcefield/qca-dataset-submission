# OpenFF Phosphate Torsion Drives v1.0

## Description

A torsion drive data set created to improve the coverage of lipid-like phosphate
parameters in Sage.

## General Information

* Date: 2024-07-17
* Class: OpenFF TorsionDrive Dataset
* Purpose: Improve phosphate coverage in Sage
* Name: OpenFF Phosphate Torsion Drives v1.0
* Number of unique molecules: 63
* Number of driven torsions: 318
* Number of filtered molecules: 0
* Number of conformers: 1526
* Number of conformers per molecule (min, mean, max): 1, 4.80, 5
* Mean molecular weight: 392.25
* Max molecular weight: 785.55
* Charges: [-5.0, -4.0, -3.0, -2.0, -1.0, 0.0]
* Dataset generator: Patrick Frankel
* Dataset submitter: Brent Westbrook

## QCSubmit Generation Pipeline

* `generate-dataset.py`: This script shows how the dataset was prepared from the
  input file `input.smi`.
* The list of SMILES in `input.smi` were curated by hand to correspond to
  phosphate motifs in lipids.

## QCSubmit Manifest

### Input files
* `input.smi`: Input SMILES strings for dataset molecules
* `generate-dataset.py`: Script describing dataset generation and submission
* `input-environment.yaml`: Environment file used to create Python environment for the script
* `full-environment.yaml`: Fully-resolved environment used to execute the script
* `td.toml`: Experimental input file for defining variables used throughout the QCA submission process

### Output Files
* `dataset.json.bz2`: Compressed dataset ready for submission
* `dataset.pdf`: Visualization of dataset molecules
* `output.smi`: SMILES strings for dataset molecules

## Metadata
* Elements: {C, S, N, H, O, P}
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
