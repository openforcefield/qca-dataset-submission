# OpenFF NAGL2 ESP Timing Benchmark v1.1

## Description
A dataset of 1005 molecules, sub-sampled from the [OpenFF multi-Br ESP Fragment Conformers v1.1](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2023-11-30-OpenFF-multi-Br-ESP-Fragment-Conformers-v1.1-single-point), the [OpenFF Iodine Fragment Opt v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2024-09-10-OpenFF-Iodine-Fragment-Opt-v1.0), and the [OpenFF ESP Fragment Conformers v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2022-01-16-OpenFF-ESP-Fragment-Conformers-v1.0) datasets. Calculated at the PBE0/def2-TZVPPD level of theory and wb97x-D/def2-TZVPPD level of theory, in vacuum and PCM water using the DDX package in Psi4. Wavefunctions were saved for the PBE0 vacuum dataset, but not the others. This is because a) we just want to benchmark how much storage space we'd need, and the different sets of wavefunctions should require the same amount of space and b) we expect there may be some errors due to possible numerical issues with I + diffuse functions + PCM that we may need to debug, and don't want to store potentially problematic wavefunctions.

Sub-sampling was done in `subsample_esp_ds.ipynb` by selecting 945 conformers at random from the ESP Fragment dataset, 53 from the I fragment dataset, and 11 from the multi-Br dataset.

The purpose of this dataset is to estimate timing and storage requirements for this methodology, which is a candidate for NAGL2 dataset generation.

## General information
* Date: 2024-09-18
* Class: OpenFF Basic Dataset
* Purpose: ESP/electric field generation
* Name: OpenFF NAGL2 ESP Timing Benchmark v1.1
* Number of unique molecules: 1005
* Number of conformers: 1009
* Number of conformers (min, mean, max): 1, 1.004, 2
* Molecular weight (min, mean, max): 58.10, 155.15, 329.82
* Charges: -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0
* Dataset submitter: Alexandra McIsaac
* Dataset generator: Alexandra McIsaac

## QCSubmit generation pipeline
* `subsample_esp_ds.ipynb` was used to reduce the parent datasets into a small dataset to use for benchmarking purposes.
* `br_subsample_filtered.json`, `esp_subsample_filtered.json`, `i_subsample_filtered.json` sub-sampled datasets pulled from the multi-Br, ESP50k, and I fragment datasets, produced by `subsample_esp_ds.ipynb`. Used as inputs for `generate-dataset.ipynb`
* `generate-dataset.ipynb` was used to transform the small sub-sampled dataset into a new dataset for submission.

## QCSubmit Manifest
* `dataset.json.bz2`: compressed dataset ready for submission
* `dataset.pdf`: Visualization of dataset molecules
* `dataset.smi`: Smiles strings for dataset molecules
* `subsample_esp_ds.ipynb`: Notebook showing how dataset was pruned from larger parent datasets.
* `generate-dataset.ipynb`: Notebook describing dataset generation and submission
* `input-environment.yaml`: Environment file used to create Python environment for the notebook
* `input-environment-full.yaml`: Fully-resolved environment used to execute the notebook.
* `br_subsample_filtered.json`: subset of multi-Br dataset selected for this benchmark
* `i_subsample_filtered.json`: subset of I fragment dataset selected for this benchmark
* `esp_subsample_filtered.json`: subset of ESP50k dataset selected for this benchmark
* `generate-compute.ipynb`: Notebook used to generate the compute expansion `compute.json`
* `compute.json`: compute expansion specs
* `generate-compute-wb97xd`
* `compute2.json`: second compute expansion specs

## Metadata
* elements: {'H', 'P', 'C', 'N', 'O', 'S', 'Cl', 'I', 'Br', 'F'}
* unique molecules: 1005
* Spec: pbe0/def2-TZVPPD
	* basis: def2-TZVPPD
	* implicit_solvent: None
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
	* Spec: wb97x-d/def2-TZVPPD/ddx-water
			* basis: def2-TZVPPD
			* implicit_solvent: {'ddx_model': 'pcm', 'ddx_radii_scaling': 1.1, 'ddx_radii_set': 'uff', 'ddx_solvent_epsilon': 78.4, 'ddx_solvent': 'water'}
			* keywords: {'dft_spherical_points': 590, 'dft_radial_points': 99, 'dft_bs_radius_alpha': 5.0}
			* maxiter: 200
			* method: wb97x-d
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
	* Spec: wb97x-d/def2-TZVPPD
			* basis: def2-TZVPPD
			* implicit_solvent: None
			* keywords: {'dft_spherical_points': 590, 'dft_radial_points': 99, 'dft_bs_radius_alpha': 5.0}
			* maxiter: 200
			* method: wb97x-d
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
