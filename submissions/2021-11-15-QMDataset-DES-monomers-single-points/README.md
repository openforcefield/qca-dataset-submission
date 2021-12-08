## Description

This is a single point energy calculation of DES monomers. Detailed description on how the data is generated can be found at https://github.com/openmm/qmdataset/tree/main/des370k.

## General Information

 - Date: 2021.11.15
 - Class: Basic dataset 
 - Purpose: Energy calculation
 - Collection: BasicDataset
 - Name: SPICE DES Monomers Single Points Dataset v1.0
 - Number of unique molecules:        374
 - Number of filtered molecules:      0
 - Number of conformers:              18700
 - Number of conformers min mean max: 50  50.00 50
 - Set of charges: [-1.0, 0.0, 1.0]
 - Dataset Submitter: Josh Horton, Pavan Behara, David Dotson
 - Dataset Generator: Peter Eastman
 - Dataset Source: https://github.com/openmm/qmdataset/tree/main/des370k

## Changelog

Here any information regarding dataset changes are recorded.

## QCSubmit generation pipeline

 - `Dataset_Generation.ipynb`: This notebook shows how the dataset was prepared from the input files. 
 
## QCSubmit Manifest

- `Dataset_Generation.ipynb`: Dataset creation notebook with instructions for submission.
- `dataset.json.xz`: The compressed gradient dataset ready for submission.
- `des370k.smi`: The smiles file of the peptide molecules.
- `des370k.pdf`: A pdf file containing molecule 2D structures.
- `des370k.hdf5.tar.gz`: HDF5 file which contains structures and is the main input for the dataset
 
## Metadata

- elements {'I', 'C', 'Br', 'P', 'Cl', 'H', 'S', 'O', 'F', 'N'}
- unique molecules 374
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
    
