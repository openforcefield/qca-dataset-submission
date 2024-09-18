# OpenFF NAGL2 ESP Timing Benchmark ddx v1.0

## Description
A dataset of 1005 molecules from the [OpenFF NAGL2 ESP Timing Benchmark](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2024-09-06-OpenFF-NAGL2-ESP-Timing-Benchmark-v1.0) dataset, calculated at the PBE0/def2-TZVPPD/ddx water level of theory.

The purpose of this dataset is to estimate timing and storage requirements for this methodology, which is a candidate for NAGL2 dataset generation.

## General information
* Date: 2024-09-16
* Class: OpenFF Basic Dataset
* Purpose: ESP/electric field generation
* Name: OpenFF NAGL2 ESP Timing Benchmark ddx v1.0
* Number of unique molecules: 1005
* Number of conformers: 1009
* Number of conformers (min, mean, max): 1, 1.004, 2
* Molecular weight (min, mean, max): 58.10, 155.15, 329.82
* Charges: -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0
* Dataset submitter: Alexandra McIsaac
* Dataset generator: Alexandra McIsaac

## QCSubmit generation pipeline
* `br_subsample_filtered.json`, `esp_subsample_filtered.json`, `i_subsample_filtered.json` sub-sampled datasets pulled from the multi-Br, ESP50k, and I fragment datasets, produced as described in the [OpenFF NAGL2 ESP Timing Benchmark](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2024-09-06-OpenFF-NAGL2-ESP-Timing-Benchmark-v1.0) repo and used as inputs for `generate-dataset.ipynb` 
* `generate-dataset.ipynb` was used to transform the small sub-sampled dataset into a new dataset for submission.

## QCSubmit Manifest
* `dataset.json.bz2`: compressed dataset ready for submission
* `dataset.pdf`: Visualization of dataset molecules
* `dataset.smi`: Smiles strings for dataset molecules
* `generate-dataset.ipynb`: Notebook describing dataset generation and submission
* `input-environment.yaml`: Environment file used to create Python environment for the notebook
* `input-environment-full.yaml`: Fully-resolved environment used to execute the notebook.
* `br_subsample_filtered.json`: subset of multi-Br dataset selected for this benchmark
* `i_subsample_filtered.json`: subset of I fragment dataset selected for this benchmark
* `esp_subsample_filtered.json`: subset of ESP50k dataset selected for this benchmark

## Metadata
* Elements: {Br, F, I, N, Cl, C, P, H, S, O}
* Spec: pbe0/def2-TZVPPD/ddx-water
	* basis: def2-TZVPPD
	* implicit_solvent: {'ddx_model': 'pcm', 'ddx_radii_scaling': 1.1, 'ddx_radii_set': 'uff', 'ddx_solvent_epsilon': 78.4, 'ddx_solvent': 'water'}
	* keywords: {'dft_spherical_points': 590, 'dft_radial_points': 99}
	* maxiter: 200
	* method: pbe0
	* program: psi4
	* SCF properties:
		* dipole
		* quadrupole
		* lowdin_charges
		* mulliken_charges
		* mbis_charges
		* mayer_indices
		* wiberg_lowdin_indices
		* dipole_polarizabilities
