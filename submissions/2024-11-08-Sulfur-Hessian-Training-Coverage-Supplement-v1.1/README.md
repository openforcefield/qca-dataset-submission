# OpenFF Sulfur Hessian Training Coverage Supplement v1.1

## Description

A basic data set created to improve the training coverage of sulfonic and
phosphonic acids, sulfone, sulfonate, sulfinyl, sulfoximine, sulfonamides,
thioether, and 1,3-thiazole groups. The structures in this data set are the
optimized geometries from `OpenFF Sulfur Optimization Training Coverage
Supplement v1.0`.

## General Information

* Date: 2024-11-08
* Class: OpenFF Optimization Dataset
* Purpose: Improve coverage in Sage
* Name: OpenFF Sulfur Hessian Training Coverage Supplement v1.1
* Number of unique molecules: 129
* Number of filtered molecules: 0
* Number of conformers: 899
* Number of conformers per molecule (min, mean, max): 1, 6.97, 10
* Mean molecular weight: 218.80
* Max molecular weight: 493.37
* Charges: [-2.0, -1.0, 0.0]
* Dataset submitter: Brent Westbrook
* Dataset generator: Brent Westbrook

## QCSubmit Generation Pipeline

* `generate.py`: This script shows how the dataset was prepared.


## QCSubmit Manifest

* `generate.py`: Script describing dataset generation and submission
* `input-environment.yaml`: Environment file used to create the Python environment for the script
* `full-environment.yaml`: Fully-resolved environment used to execute the script
* `opt.toml`: Experimental [qcaide](https://github.com/ntBre/qcaide) input file for defining
variables used throughout the QCA submission process
* `dataset.json.bz2`: Compressed dataset ready for submission
* `dataset.pdf`: Visualization of dataset molecules
* `dataset.smi`: SMILES strings for dataset molecules

## Metadata

* Elements: {O, S, C, Cl, P, N, F, Br, H}
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

## Changelog
v1.0 included a manual implementation of
`OptimizationResultCollection.create_basic_dataset` that failed to preserve
QCArchive molecule IDs between the optimization and single-point datasets.
Unfortunately, this issue would not have been avoided by that version of
`create_basic_dataset` either. The issue has been fixed in openff-qcsubmit
[version
0.54](https://github.com/openforcefield/openff-qcsubmit/releases/tag/0.54.0), so
the environment has been updated to use this release, and the `generate.py`
script has been updated to use `create_basic_dataset`.
