# OpenFF Unconstrained Protein Fragments

### Description

This is a protein fragment dataset consisting of optimizations on various protein fragments prepared by David Cerutti.
We have 12 central residues capped with a combination of different terminal residues.

Details from David Cerutti, referring to motivations for these selections:

> The way it's set up, we are scanning phi and psi of the central residue with a random selection of ALA, GLY, SER, or VAL to the N- or C-terminus,
> and the customary ACE and NME blocking groups outside of that.

### General Information

 - Date: 2020.10.27
 - Class: OpenFF Optimization
 - Purpose: Optimizations for common amino acids
 - Collection: OptimizationDataset
 - Name: OpenFF Unconstrained Protein Fragments v2.1
 - Number of unique molecules: 185
 - Number of optimizations: 6709
 - Submitter: David Dotson
 
Note, each folder contains molecules saved via mol2 in each confirmation; however, the bond order is incorrect.
We let openeye interpret it by re-saving to PDB first before creating the dataset.

### QCSubmit generation pipeline

 - `Dataset Preparation.ipynb`: This notebook shows how the Optimization dataset was prepared from the input files. 
 
### QCSubmit Manifest

- `Dataset Preparation.ipynb`: Dataset creation notebook with details on decisions made for submission.
- `dataset.json.bz2`: The compressed Optimization dataset ready for submission.
- `protein_dataset.pdf`: A pdf file containing molecule 2D structures.
- `../2020-08-12-OpenFF-Protein-Fragments-version2/Input_files.tar.gz`: Folders containing the input molecule conformations and the (unused) corresponding dihedral restraints.
- `molecules.smi`: SMILES for every molecule in the submission
 
### Metadata

- elements {'C', 'H', 'N', 'O', 'S'}
- unique molecules: 185
- optimizations: 6709
- scf properties:
    - `dipole`
    - `quadrupole`
    - `wiberg_lowdin_indices`
    - `mayer_indices`
- qc spec
    - name: default
    - method: B3LYP-D3BJ
    - basis: DZVP
