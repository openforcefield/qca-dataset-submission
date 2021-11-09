## Description

This is a single point energy calculation of solvated amino acids. More description on how the data is generated can be found at https://github.com/openmm/qmdataset/tree/main/solvated-amino-acids.

## General Information

 - Date: 2021.11.08
 - Class: Basic dataset 
 - Purpose: Energy calculation
 - Collection: BasicDataset
 - Name: SPICE Solvated Amino acids
 - Number of unique molecules: 26
 - Number of filtered molecules      0
 - Number of conformers              1300
 - Number of conformers min mean max 50  50.00 50
 - Set of charges: [-1.0, 0.0, 1.0]
 - Dataset Submitter: Josh Horton/Pavan Behara/David Dotson
 - Dataset Generator: Peter Eastman
 - Dataset Source: https://github.com/openmm/qmdataset/tree/main/dipeptides

## Changelog

Here any information regarding dataset changes are recorded.

## QCSubmit generation pipeline

 - `Dataset_Generation.ipynb`: This notebook shows how the dataset was prepared from the input files. 
 
## QCSubmit Manifest

- `Dataset_Generation.ipynb`: Dataset creation notebook with instructions for submission.
- `dataset.json.xz`: The compressed constrained optimization dataset ready for submission.
- `solvated-amino-acids.smi`: The smiles file of the peptide molecules.
- `solvated-amino-acids.pdf`: A pdf file containing molecule 2D structures.
- `solvated-amino-acids.hdf5`: HDF5 file which contains structures and is the main input for the dataset
 
## Metadata

- elements {'N', 'S', 'O', 'C', 'H'}
- unique molecules 26
- scf properties:
    - dipole
    - quadrupole
    - wiberg_lowdin_indices
    - mayer_indices
    - mbis_charges
- qc spec
    - name: wB97M-D3BJ/def2-TZVPPD
    - method: wB97M-D3BJ
    - basis: def2-TZVPPD
    
