# OpenFF Torsion Coverage Supplement v1.0

## Description

A torsion drive data set created to improve the coverage of both existing Sage
2.1.0 proper torsion parameters and new parameters added through the [torsion
multiplicity][tm] project. The molecules in this data set were selected from the
[ChEMBL 33][chembl] database with additional processing steps described below.

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
* The list of proper torsion parameter ID and SMILES pairs in `all.smiles` were
  collected by searching the ChEMBL database for all of the molecules matching
  the parameters of interest, fragmenting these molecules with the
  [RECAP][recap] algorithm as implemented in RDKit, clustering them based on
  their 1024-bit [Morgan][morgan] fingerprints using the [DBSCAN][dbscan]
  algorithm, and taking the two smallest molecules from these cluster centroids.
  The code used for all of these steps can be found [here][chembl-search].

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
[recap]: https://www.rdkit.org/docs/source/rdkit.Chem.Recap.html
[morgan]: https://www.rdkit.org/docs/cppapi/namespaceRDKit_1_1MorganFingerprints.html
[dbscan]: https://en.wikipedia.org/wiki/DBSCAN
[chembl-search]: https://github.com/ntBre/chembl-search
