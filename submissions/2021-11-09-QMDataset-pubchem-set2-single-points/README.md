## Description

This is a single point energy calculation of pubchem set2 (2501-5000) molecules. Detailed description on how the data is generated can be found at https://github.com/openmm/qmdataset/tree/main/pubchem.

## General Information

 - Date: 2021.11.09
 - Class: Basic dataset 
 - Purpose: Energy calculation
 - Collection: BasicDataset
 - Name: SPICE PubChem Set 2 Single Points Dataset
 - Number of unique molecules:        2431
 - Number of filtered molecules:      0
 - Number of conformers:              121540
 - Number of conformers min mean max: 40  50.00 50
 - Set of charges: [-6.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0]
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
- `pubchem-set2.smi`: The smiles file of the peptide molecules.
- `pubchem-set2.pdf`: A pdf file containing molecule 2D structures.
- `pubchem-set2.hdf5.tar.gz`: HDF5 file which contains structures and is the main input for the dataset
 
## Metadata

- elements {'H', 'P', 'C', 'Cl', 'Br', 'N', 'F', 'S', 'O', 'I'}
- unique molecules 2431
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
