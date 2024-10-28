# MLPepper-RECAP-Optimized-Fragments-Add-Iodines-v1.0

## Description

A single point dataset created by combining the [50k ESP from Simon](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2022-01-16-OpenFF-ESP-Fragment-Conformers-v1.0) and 
[Br substituted set from Lily](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2023-11-30-OpenFF-multi-Br-ESP-Fragment-Conformers-v1.1-single-point), filtering by Cl and Br and replacing them successively with iodines, as well as some additional iodines from [Lexie and Lily's set](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2024-09-10-OpenFF-Iodine-Fragment-Opt-v1.0).
Each fragment had 5 conformations generated which were optimised locally using an AIMNET2 model trained to `wb97m-d3`. 
This adds an extension of iodines molecules to the [original mlpepper dataset](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2024-07-26-MLPepper-RECAP-Optimized-Fragments-v1.0).

The aim of the dataset is to provide polarised and gas phase electrostatic properties which can be used to generate ML models 
for partial charge prediction. Unlike past datasets the wavefunction will not be saved to recompute the ESP instead we recommend building the ESP 
from the MBIS atomic multipoles which save substantial amount of space. 

An off equilibrium data set will also be generated to enable conformation dependent prediction of charges. 

## General Information


* Date: 2024-10-11
* Class: OpenFF SinglePoint Dataset
* Purpose: Electrostatic properties for ML prediction models 
* Name: MLPepper RECAP Optimized Fragments v1.0 Add Iodines
* Number of unique molecules: 5733
* Number of filtered molecules: 0
* Number of conformers: 6131
* Number of conformers per molecule (min, mean, max): 1, 1.07, 3
* Mean molecular weight: 278.86
* Max molecular weight: 701.59
* Charges: [-4.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0]
* Dataset submitter: Josh Horton/ Charlie Adams
* Dataset generator: Josh Horton/ Charlie Adams

## QCSubmit generation pipeline

- `create_optimisation.py`: Was used to create the iodine Cl and Br replacements from the original MLPepper, and 
then optimisation dataset.
- `create_singlepoints.py`: Was used to create the singlepoints dataset for the optimised iodine sets.
- `create_dataset.py`: Finally this script combines the resulting datasets into a single point dataset ready for submission.

## QCSubmit Manifest

### Input Files

- `create_optimisation.py`: Script used to make the optimisation dataset for local optimisation. 
- `create_singlepoints.py`: Script to create the singlepoints dataset from the optimised geometries.
- `create_dataset.py`: Script to create the singlepoint dataset from the optimization set, removing any connectivity issues. 

### Output Files
- `dataset.json.bz2`: The basic dataset ready for submission.
- `dataset.pdf`: A pdf file containing molecule 2D structures.
- `dataset.smi`: SMILES for every molecule in the submission.
- `dataset_mlpepper.smi`: SMILES of the original dataset to generate the Iodines

### Metadata

* Number of conformers: 6131
* Number of conformers per molecule (min, mean, max): 1, 1.07, 3
* Mean molecular weight: 278.86
* Max molecular weight: 701.59
* Charges: [-4.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0]
* Elements: {B, Si, I, C, N, Br, O, S, Cl, H, F, P}
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
        * implicit_solvent: {'ddx_model': 'pcm', 'ddx_radii_scaling': 1.1, 'ddx_radii_set': 'uff', 'ddx_solvent_epsilon': 78.4, 'ddx_solvent': 'water'}
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
