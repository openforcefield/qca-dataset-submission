## Description

This is a single point energy calculation of pubchem set 5 (10001-12500) molecules. Detailed description on how the data is generated can be found at https://github.com/openmm/qmdataset/tree/main/pubchem.

 v1.2 removes wavefunction storage, as this uses a substantial amount of space on public QCArchive and was not desired for this dataset.

## General Information

 - Date: 2021.11.09
 - Class: Basic dataset 
 - Purpose: Energy calculation
 - Collection: BasicDataset
 - Name: SPICE PubChem Set 5 Single Points Dataset v1.2
 - Number of unique molecules:        2463
 - Number of filtered molecules:      0
 - Number of conformers:              123150
 - Number of conformers min mean max: 50  50.00 50
 - Set of charges: [-5.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0]
 - Dataset Submitter: Josh Horton/Pavan Behara/David Dotson
 - Dataset Generator: Peter Eastman
 - Dataset Source: https://github.com/openmm/qmdataset/tree/main/pubchem

## Changelog

Here any information regarding dataset changes are recorded.

## QCSubmit generation pipeline

 - `Dataset_Generation.ipynb`: This notebook shows how the dataset was prepared from the input files. 
 
## QCSubmit Manifest

- `Dataset_Generation.ipynb`: Dataset creation notebook with instructions for submission.
- `dataset.json.bz2`: The compressed gradient dataset ready for submission.
- `pubchem-set5.smi`: The smiles file of the peptide molecules.
- `pubchem-set5.pdf`: A pdf file containing molecule 2D structures.
- `pubchem-set5.hdf5.tar.gz`: HDF5 file which contains structures and is the main input for the dataset
 
## Metadata

- elements {'F', 'H', 'S', 'Br', 'Cl', 'N', 'P', 'C', 'I', 'O'}
- unique molecules 2463
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
