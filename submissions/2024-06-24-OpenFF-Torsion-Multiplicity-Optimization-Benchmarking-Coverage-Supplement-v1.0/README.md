# OpenFF Torsion Multiplicity Optimization Benchmarking Coverage Supplement v1.0

## Description
An optimization data set created to improve the benchmarking coverage of both
existing Sage 2.2.0 proper torsion parameters and new parameters added through
the torsion multiplicity project. The molecules in this data set were partly
selected from the ChEMBL 33 database and partly generated manually to match
parameters not covered by ChEMBL.

## General Information

* Date: 2024-06-24
* Class: OpenFF Optimization Dataset
* Purpose: Improver proper torsion coverage in Sage
* Name: OpenFF Torsion Multiplicity Optimization Benchmarking Coverage Supplement v1.0
* Number of unique molecules: 73
* Number of filtered molecules: 0
* Number of conformers: 451
* Number of conformers per molecule (min, mean, max): 1, 6.18, 10
* Mean molecular weight: 214.06
* Max molecular weight: 489.48
* Charges: [-1.0, 0.0, 1.0, 2.0]
* Dataset submitter: Brent Westbrook
* Dataset generator: Brent Westbrook

## QCSubmit Generation Pipeline

* `generate-dataset.py`: This script shows how the dataset was prepared from the input files:
`bench.opt.smi` and `ff.offxml`.
* The list of proper torsion parameter ID and SMILES pairs in `bench.opt.smi` were
collected by searching the ChEMBL database for all of the molecules matching the
parameters of interest. The code used for all of these steps can be found
[here](https://github.com/ntBre/valence-fitting/tree/c1e98fb20e7a4c9622ff031d8b23fb0b1846be7d/02_curate-data/frag).

## QCSubmit Manifest

### Input Files
* `generate-dataset.py`: Script describing dataset generation and submission
* `input-environment.yaml`: Environment file used to create Python environment for the script
* `full-environment.yaml`: Fully-resolved environment used to execute the script
* `ff.offxml`: Draft force field with Sage 2.2.0 proper torsions split to ensure single multiplicities
* `opt.toml`: Experimental input file for defining variables used throughout the QCA submission process

### Output Files
* `dataset.json.bz2`: Compressed dataset ready for submission
* `dataset.pdf`: Visualization of dataset molecules
* `dataset.smi`: Smiles strings for dataset molecules

## Metadata
* Elements: {Cl, H, I, S, O, N, Br, C, P}
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
