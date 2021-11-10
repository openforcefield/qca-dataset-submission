## Description

This is a single point energy calculation of DES370K molecules. Detailed description on how the data is generated can be found at https://github.com/openmm/qmdataset/tree/main/des370k.

## General Information

 - Date: 2021.11.09
 - Class: Basic dataset 
 - Purpose: Energy calculation
 - Collection: BasicDataset
 - Name: SPICE DES370K Single Points Dataset
 - Number of unique molecules:        3407
 - Number of filtered molecules:      0
 - Number of conformers:              345682
 - Number of conformers min mean max: 8 101.46 4869
 - Set of charges: [-1.0, 0.0, 1.0, 2.0]
 - Dataset Submitter: Josh Horton, Pavan Behara, David Dotson
 - Dataset Generator: Peter Eastman
 - Dataset Source: https://github.com/openmm/qmdataset/tree/main/des370k

## Changelog

Here any information regarding dataset changes are recorded.

## QCSubmit generation pipeline

 - `Dataset_Generation.ipynb`: This notebook shows how the dataset was prepared from the input files. 
 
## QCSubmit Manifest

- `Dataset_Generation.ipynb`: Dataset creation notebook with instructions for submission.
- `dataset.json.xz`: The compressed constrained optimization dataset ready for submission.
- `des370k.smi`: The smiles file of the peptide molecules.
- `des370k.pdf`: A pdf file containing molecule 2D structures.
- `des370k.hdf5.tar.gz`: HDF5 file which contains structures and is the main input for the dataset
 
## Metadata

- elements {'N', 'O', 'Mg', 'H', 'F', 'K', 'Br', 'Na', 'P', 'Cl', 'I', 'Ca', 'S', 'Li', 'C'} 
- unique molecules 3407
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
    
