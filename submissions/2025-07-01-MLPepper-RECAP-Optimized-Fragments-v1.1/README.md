# MLPepper RECAP Optimized Fragments v1.1

## Description
A single point dataset created by combining the [50k ESP from Simon](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2022-01-16-OpenFF-ESP-Fragment-Conformers-v1.0) and 
[Br substituted set from Lily](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2023-11-30-OpenFF-multi-Br-ESP-Fragment-Conformers-v1.1-single-point). 
This was then extended to cover Boron and Silicon by applying the RECAP decomposition scheme in RDKit to the [B & Si pubchem set from SPICE2](https://github.com/openmm/spice-dataset/blob/main/pubchem/pubchem-boron-silicon.hdf5). This dataset now contains additional iodine containing molecules.
Each fragment had 5 conformations generated which were optimised locally using an AIMNET2 model trained to `wb97m-d3`. 

The aim of the dataset is to provide polarised and gas phase electrostatic properties which can be used to generate ML models 
for partial charge prediction. Unlike past datasets the wavefunction will not be saved to recompute the ESP instead we recommend building the ESP 
from the MBIS atomic multipoles which save substantial amount of space. 

An off equilibrium data set will also be generated to enable conformation dependent prediction of charges.

## General Information
 
* Date: 2025-07-01
* Class: OpenFF SinglePoint Dataset
* Purpose: Single point property calculations for charge models, expanded to include iodine
* Name: MLPepper RECAP Optimized Fragments v1.1
* Number of unique molecules: 50618
* Number of filtered molecules: 0
* Number of conformers: 68966
* Number of conformers per molecule (min, mean, max): 1, 1.36, 5
* Mean molecular weight: 150.10
* Max molecular weight: 466.59
* Charges: [-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0]
* Dataset submitter: Jennifer A Clark
* Dataset curator: Josh Horton/ Charlie Adams
* Dataset generator: Josh Horton/ Charlie Adams

## QCSubmit Generation Pipeline
* `main.ipynb`: A jupyter notebook which shows how the optimization dataset was prepared from the input files.

## QCSubmit Manifest

### Input Files
* `main.ipynb`: A jupyter notebook which shows how the optimization dataset was prepared from the input files.
* `ds_info.json`: Metadata information for optimization dataset record imported by `main.ipynb`
### Output Files
* `dataset.json.bz2`: Compressed dataset ready for submission
* `dataset.pdf.zip`: Visualization of dataset molecules compressed due to size
* `dataset.smi`: Smiles strings for dataset molecules

## Metadata
* Elements: {F, P, Br, O, I, N, Si, S, B, Cl, H, C}
* Spec: wb97x-d/def2-tzvpp
	* basis: def2-tzvpp
	* implicit_solvent: None
	* keywords: {'dft_spherical_points': 590, 'dft_radial_points': 99}
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
* Spec: wb97x-d/def2-tzvpp/ddx-water
	* basis: def2-tzvpp
	* implicit_solvent: {'ddx_model': 'pcm', 'ddx_radii_scaling': 1.1, 'ddx_radii_set': 'uff', 'ddx_solvent_epsilon': None, 'ddx_solvent': 'water'}
	* keywords: {'dft_spherical_points': 590, 'dft_radial_points': 99}
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
