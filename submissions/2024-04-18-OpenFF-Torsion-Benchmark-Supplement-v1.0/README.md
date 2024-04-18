# OpenFF Torsion Benchmark Supplement v1.0

## Description

An optimization data set created to expand the benchmarking data for both
existing Sage 2.2.0 proper torsion parameters and new parameters added through
the [torsion multiplicity][tm] project. The molecules in this data set were The
molecules in this data set were selected from the [ChEMBL 33][chembl] database
with additional processing steps described below.

## General Information

* Date: 2024-04-17
* Class: OpenFF Optimization Dataset
* Purpose: Improve proper torsion training data in Sage
* Name: OpenFF Torsion Drive Supplement v1.0
* Number of unique molecules: 51
* Number of filtered molecules: 0
* Number of conformers: 51
* Number of conformers (min, mean, max): 1, 1, 1
* Mean molecular weight: 259.38
* Max molecular weight: 508.31
* Charges: [-1.0, 0.0, 1.0, 2.0]
* Dataset submitter: Brent Westbrook
* Dataset generator: Brent Westbrook

## QCSubmit generation pipeline

* `generate-dataset.ipynb`: This notebook shows how the dataset was prepared
  from the input files: `all.smiles` and `tm-2.2.offxml` (Sage 2.2.0 release
  candidate with [torsion multiplicity][tm] changes applied).
* The list of proper torsion parameter ID and SMILES pairs in `all.smiles` were
  collected by searching the ChEMBL database for all of the molecules matching
  the parameters of interest, fragmenting these molecules with the
  [RECAP][recap] algorithm as implemented in RDKit, clustering them based on
  their 1024-bit [Morgan][morgan] fingerprints using the [DBSCAN][dbscan]
  algorithm, and taking the smallest molecules from these cluster centroids. The
  code used for all of these steps can be found [here][cura].

## QCSubmit Manifest

* `dataset.json.bz2`: Compressed dataset ready for submission
* `dataset.pdf`: Visualization of dataset molecules
* `dataset.smi`: Smiles strings for dataset molecules
* `generate-dataset.ipynb`: Notebook describing dataset generation and submission
* `input-environment.yaml`: Environment file used to create Python environment for the notebook
* `full-environment.yaml`: Fully-resolved environment used to execute the notebook

## Metadata

* elements: H, C, N, O, F, P, S, Cl, Br
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
[cura]: https://github.com/ntBre/cura
