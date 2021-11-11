## Description

This is a single point energy calculation of pubchem set3 (5001-7500) molecules. Detailed description on how the data is generated can be found at https://github.com/openmm/qmdataset/tree/main/pubchem.

## General Information

 - Date: 2021.11.09
 - Class: Basic dataset 
 - Purpose: Energy calculation
 - Collection: BasicDataset
 - Name: SPICE PubChem Set 3 Single Points Dataset v1.0
 - Number of unique molecules:        2446
 - Number of filtered molecules:      0
 - Number of conformers:              122226
 - Number of conformers min mean max: 6  49.97 50
 - Set of charges: [-4.0, -2.0, -1.0, 0.0, 1.0, 2.0]
 - Dataset Submitter: Josh Horton, Pavan Behara, David Dotson
 - Dataset Generator: Peter Eastman
 - Dataset Source: https://github.com/openmm/qmdataset/tree/main/pubchem

## Changelog

Here any information regarding dataset changes are recorded.

## QCSubmit generation pipeline

 - `Dataset_Generation.ipynb`: This notebook shows how the dataset was prepared from the input files. 
 
## QCSubmit Manifest

- `Dataset_Generation.ipynb`: Dataset creation notebook with instructions for submission.
- `dataset.json.bz2`: The compressed constrained optimization dataset ready for submission.
- `pubchem-set3.smi`: The smiles file of the peptide molecules.
- `pubchem-set3.pdf`: A pdf file containing molecule 2D structures.
- `pubchem-set3.hdf5.tar.gz`: HDF5 file which contains structures and is the main input for the dataset
 
## Metadata

- elements {'N', 'C', 'S', 'Cl', 'Br', 'F', 'P', 'I', 'H', 'O'}
- unique molecules 2446
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
