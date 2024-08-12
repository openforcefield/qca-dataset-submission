# OpenFF Alkane Torsion Drives v1.0

## Description
A torsion drive data set created to improve the coverage of alkane and alkene
parameters in Sage.
## General Information

* Date: 2024-08-09
* Class: OpenFF TorsionDrive Dataset
* Purpose: Improve proper torsion coverage in Sage
* Name: OpenFF Alkane Torsion Drives v1.0
* Number of unique molecules: 50
* Number of driven torsions: 192
* Number of filtered molecules: 0
* Number of conformers: 250
* Number of conformers per molecule (min, mean, max): 1, 1.30, 3
* Mean molecular weight: 92.04
* Max molecular weight: 126.24
* Charges: [0.0]
* Dataset submitter: Brent Westbrook
* Dataset generator: Michael Shirts

## QCSubmit Generation Pipeline
* `generate-dataset.py`: This script shows how the dataset was prepared from the input file `input.smi`.
* The list of SMILES in `input.smi` were curated by hand to correspond to
alka/ene parameters missing coverage in Sage.

## QCSubmit Manifest
### Input
* `generate-dataset.py`: Script describing dataset generation and submission
* `input.smi`: Input SMILES
* `input-environment.yaml`: Environment file used to create Python environment for the script
* `full-environment.yaml`: Fully-resolved environment used to execute the script
* `td.toml`: Experimental input file for defining variables used throughout the QCA submission process
### Output
* `dataset.json.bz2`: Compressed dataset ready for submission
* `dataset.pdf`: Visualization of dataset molecules
* `dataset.smi`: Smiles strings for dataset molecules

## Metadata
* Elements: {C, H}
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
