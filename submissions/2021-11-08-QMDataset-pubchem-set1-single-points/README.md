## Description

This is a single point energy calculation of PubChem set 1 (1-2500) molecules. Detailed description on how the data is generated can be found at https://github.com/openmm/qmdataset/tree/main/pubchem.

v1.2 removes wavefunction storage, as this uses a substantial amount of space on public QCArchive and was not desired for this dataset.

## General Information

 - Date: 2021.11.08
 - Class: Basic dataset 
 - Purpose: Energy calculation
 - Collection: BasicDataset
 - Name: SPICE PubChem Set 1 Single Points Dataset v1.1
 - Number of unique molecules:        2372
 - Number of filtered molecules:      0
 - Number of conformers:              118606
 - Number of conformers min mean max: 6  50.00 100
 - Set of charges: [-8.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0]
 - Dataset Submitter: Josh Horton, Pavan Behara, David Dotson
 - Dataset Generator: Peter Eastman
 - Dataset Source: https://github.com/openmm/qmdataset/tree/main/pubchem

## Changelog

Here any information regarding dataset changes are recorded.

## QCSubmit generation pipeline

 - `Dataset_Generation.ipynb`: This notebook shows how the dataset was prepared from the input files. 
 
## QCSubmit Manifest

- `Dataset_Generation.ipynb`: Dataset creation notebook with instructions for submission.
- `dataset.json.bz2`: The compressed gradient dataset ready for submission.
- `pubchem-set1.smi`: The smiles file of the peptide molecules.
- `pubchem-set1.pdf`: A pdf file containing molecule 2D structures.
- `pubchem-set1.hdf5.tar.gz`: HDF5 file which contains structures and is the main input for the dataset
 
## Metadata

- elements {'O', 'Cl', 'N', 'C', 'P', 'Br', 'S', 'F', 'I', 'H'}
- unique molecules 2372
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
