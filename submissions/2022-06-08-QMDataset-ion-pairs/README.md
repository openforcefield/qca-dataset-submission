## Description

This is a single point energy calculation of ion pairs of monoatomic species. Detailed description on how the data is generated can be found at https://github.com/openmm/qmdataset/tree/main/ions.


## General Information

 - Date: 2022.06.08
 - Class: Basic dataset 
 - Purpose: Energy calculation
 - Collection: BasicDataset
 - Name: SPICE Ion Pairs Single Points Dataset v1.0
 - Number of unique molecules:        28
 - Number of filtered molecules:      0
 - Number of conformers:              1428
 - Number of conformers min mean max: 51  51.00 51
 - Set of charges: [-2.0, 0.0, 2.0]
 - Dataset Submitter: Josh Horton, Pavan Behara, David Dotson
 - Dataset Generator: Peter Eastman
 - Dataset Source: https://github.com/openmm/qmdataset/tree/main/ions

## Changelog

Here any information regarding dataset changes are recorded.

## QCSubmit generation pipeline

 - `Dataset_Generation.ipynb`: This notebook shows how the dataset was prepared from the input files. 
 
## QCSubmit Manifest

- `Dataset_Generation.ipynb`: Dataset creation notebook with instructions for submission.
- `dataset.json.bz2`: The compressed gradient dataset ready for submission.
- `dataset.smi`: The smiles file of the peptide molecules.
- `dataset.pdf`: A pdf file containing molecule 2D structures.
- `dataset.hdf5`: HDF5 file which contains structures and is the main input for the dataset
 
## Metadata

- elements {'F', 'Cl', 'Li', 'Na', 'Br', 'K', 'I'}
- unique molecules 28
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
