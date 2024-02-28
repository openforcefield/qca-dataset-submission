# OpenFF Torsion Coverage Supplement v1.0

## Description

A torsion drive data set created to improve the coverage of both existing Sage
2.1.0 proper torsion parameters and new parameters added through the [torsion
multiplicity][tm] project. The molecules in this data set were primarily
selected from the [ChEMBL 33][chembl] database, with three additional molecules
from the [NCI Open 2012-05-01][nci] database.

## General Information

* Date: 2024-01-31
* Class: OpenFF TorsionDrive Dataset
* Purpose: Improve proper torsion coverage in Sage
* Name: OpenFF Torsion Coverage Supplement v1.0
* Number of unique molecules: 43
* Number of filtered molecules: 0
* Number of conformers: 43
* Number of conformers (min, mean, max): 1, 1.88, 7
* Mean molecular weight: 202.95
* Max molecular weight: 431.57
* Charges: [-2, -1, 0, 1, 2]
* Dataset submitter: Brent Westbrook
* Dataset generator: Brent Westbrook

## QCSubmit generation pipeline

* `generate-dataset.ipynb`: This notebook shows how the dataset was prepared
  from the input files: `all.smiles` and `tm.v2.offxml`.

## QCSubmit Manifest

* `dataset.json.bz2`: Compressed dataset ready for submission
* `dataset.pdf`: Visualization of dataset molecules
* `dataset.smi`: Smiles strings for dataset molecules
* `generate-dataset.ipynb`: Notebook describing dataset generation and submission
* `input-environment.yaml`: Environment file used to create Python environment for the notebook
* `full-environment.yaml`: Fully-resolved environment used to execute the notebook

## Metadata

* elements: {C, Cl, F, H, N, O, S}
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

<!-- References -->
[tm]: https://openforcefield.atlassian.net/wiki/spaces/FF/pages/2603909164/Torsion+multiplicity
[chembl]: https://www.ebi.ac.uk/chembl/
[nci]: https://cactus.nci.nih.gov/download/nci/
