# OpenFF Additional Generated ChEMBL Optimizations 4.0

## Description

This optimization dataset adds more coverage for rare parameters, as listed in `inputs/`.
The parameters, as in Sage 2.2.1 are:
- a16, a36, a7
- b15, b23, b24, b29, b40, b42, b44, b47, b49, b54, b55, b59, b62, b63, b64, b65, b67, b69, b74, b76, b77, b78, b80, b81, b82
- i4, i5
- t101, t102, t103, t112, t113, t114, t12, t126, t128, t129, t136, t137, t138a, t141, t141a, t141c
- t154, t158, t164, t165, t167, t30, t31a, t33, t42a, t49, t54, t55, t59, t60, t61, t62, t7, t73, t8, t81, t88, t89.

Molecules were generated according to the process in https://github.com/lilyminium/fragment-chembl-data/commit/1a85d6b296350867a134529aa76930965a710a68
repo. In short, for torsions:

- ChEMBL molecules were split into "elementary" fragments without rotatable bonds
- Elementary fragments were combined with a single bond
- For each torsion, a pool of up to 5000 molecules were initially selected after sorting for low molecular weight
- Up to 250 molecules were selected from this pool by maximising chemical diversity, using the Tanimoto distance of the Morgan fingerprints
- Up to 10 molecules were selected from the second pool by maximising the diversity of coupled torsions through the central bond

For bonds, angles, and impropers:
- For each torsion, a pool of up to 10000 molecules were initially selected after sorting for low molecular weight
- Up to 500 molecules were selected from this pool by maximising chemical diversity, using the Tanimoto distance of the Morgan fingerprints
- Up to 50 molecules were selected from the second pool by maximising the diversity of other parameters applied to the atoms of each parameter


This dataset uses the OpenFF default level of theory (B3LYP-D3BJ/DZVP).
It covers the H, N, P, F, Br, C, I, Cl, S, O elements and -3.0, -2.0, -1.0, 0.0, 1.0, 2.0 charges.
Molecular MW ranges from 32 – 313 Da with mean MW of 126 Da.


## General Information

* Date: 2025-04-14
* Class: OpenFF Optimization Dataset
* Purpose: Improve coverage of low-coverage valence parameters with automatically selected molecules from ChEMBL fragments
* Dataset name: OpenFF Additional Generated ChEMBL Optimizations 4.0
* Number of unique molecules: 1844
* Number of filtered molecules: 3
* Number of conformers: 2429
* Number of conformers per molecule (min, mean, max): 1, 1.32, 6
* Mean molecular weight: 125.99
* Min molecular weight: 32.05
* Max molecular weight: 312.99
* Charges: [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0]
* Dataset generator: Lily Wang
* Dataset submitter: Lily Wang

## QCSubmit Generation Pipeline

* `generate-dataset.ipynb`: This notebooks shows how the dataset was prepared from the
  input files in `inputs`.

## QCSubmit Manifest

### Input files
* `generate-dataset.ipynb`: Notebook describing dataset generation and submission
* `input-env.yaml`: Environment file used to create Python environment for the script
* `full-env.yaml`: Fully-resolved environment used to execute the script

### Output Files
* `dataset.json.bz2`: Compressed dataset ready for submission
* `dataset.pdf`: Visualization of dataset molecules
* `output.smi`: SMILES strings for dataset molecules
* `inputs.tar.gz`: Compressed input files as detailed below
* `inputs/selected-*.smi`: Selected SMILES for particular parameter
* `inputs/*_molecules.pdf`: Generated molecules for particular parameter


## Metadata
* Elements: {F, Br, Cl, P, S, O, N, H, C, I}
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

