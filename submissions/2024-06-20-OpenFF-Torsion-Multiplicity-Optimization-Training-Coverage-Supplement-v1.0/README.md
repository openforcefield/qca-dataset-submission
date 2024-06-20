# OpenFF Torsion Multiplicity Optimization Training Coverage Supplement v1.0

## Description
An optimization data set created to improve the coverage of both existing Sage
2.2.0 proper torsion parameters and new parameters added through the torsion
multiplicity project. The molecules in this data set were partly selected from
the ChEMBL 33 database and partly generated manually to match parameters not
covered by ChEMBL.

## General Information

* Date: 2024-06-20
* Class: OpenFF Optimization Dataset
* Purpose: Improve proper torsion coverage in Sage
* Name: OpenFF Torsion Multiplicity Optimization Training Coverage Supplement v1.0
* Number of unique molecules: 37
* Number of filtered molecules: 0
* Number of conformers: 259
* Number of conformers per molecule (min, mean, max): 1, 7.00, 10
* Mean molecular weight: 187.31
* Max molecular weight: 489.48
* Charges: [0.0, 1.0, 2.0]
* Dataset submitter: Brent Westbrook
* Dataset generator: Brent Westbrook

## QCSubmit Generation Pipeline
* `generate-dataset.ipynb`: This notebook shows how the dataset was prepared
from the input files: `train.opt.smi` and `ff.offxml`.

* The list of proper torsion parameter ID and SMILES pairs in `train.opt.smi`
were collected by searching the ChEMBL database for all of the molecules
matching the parameters of interest. The code used for all of these steps can be
found [here][frag].

## QCSubmit Manifest

### Input Files
* `generate-dataset.ipynb`: Notebook describing dataset generation and submission
* `input-environment.yaml`: Environment file used to create Python environment for the notebook
* `full-environment.yaml`: Fully-resolved environment used to execute the notebook
* `ff.offxml`: Draft force field with Sage 2.2.0 proper torsions split to ensure single multiplicities
* `test.toml`: Experimental input file for defining variables used throughout the QCA submission process

### Output Files
* `dataset.json.bz2`: Compressed dataset ready for submission
* `dataset.pdf`: Visualization of dataset molecules
* `dataset.smi`: Smiles strings for dataset molecules

## Metadata
* Elements: {C, Cl, S, O, H, P, N, Br}
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

[frag]: https://github.com/ntBre/valence-fitting/tree/c1e98fb20e7a4c9622ff031d8b23fb0b1846be7d/02_curate-data/frag
