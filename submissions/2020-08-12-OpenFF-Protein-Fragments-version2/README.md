## Description

This is the full protein fragment dataset (version2) consisting of constrained optimizations on various protein fragments prepared by David Cerutti.
We have 12 central residues which are capped with a combination of different terminal residues.

Details from David:

> The way it's set up, we are scanning phi and psi of the central residue with a random selection of ALA, GLY, SER, or VAL to the N- or C-terminus, and the customary ACE and NME blocking groups outside of that.

## General Information

 - Date: 12/8/2020
 - Class: OpenFF constrained optimization 
 - Purpose: Optimizations 
 - Collection: OptimizationDataset
 - Name: OpenFF Protein Fragments v2.0
 - Number of unique molecules: 16
 - Submitter: Josh Horton
 
## Changelog

Here any information regarding dataset changes are recorded.

### v2.1

The dataset constraints were found to be incorrect due to a mistake in the indexing.
The indices were assumed to be 0 based but were indexed from 1; this is corrected in version 2.1.
 - Date: 2020.10.06
 - Class: OpenFF constrained optimization 
 - Purpose: Optimizations 
 - Collection: OptimizationDataset
 - Name: OpenFF Protein Fragments v2.1
 - Number of unique molecules: 16
 - Submitter: Josh Horton
 
Note, each folder contains molecules saved via mol2 in each confirmation however the bond order is incorrect, we let openeye interpret it by re-saving to PDB first before creating the dataset.

## QCSubmit generation pipeline

 - `Qcsubmit protein prep.ipynb`: This notebook shows how the optimization dataset was prepared from the input files and constraint files. 
 
## QCSubmit Manifest

- `Qcsubmit protein prep.ipynb`: Dataset creation notebook with instructions for submission.
- `dataset.json`: The optimization dataset ready for submission.
- `fragments.smi`: The smiles file of the biaryl molecules.
- `protein_dataset.pdf`: A pdf file containing molecule 2D structures.
- `Input_files`: Folders containing the input molecule conformations and the corresponding dihedral restraints.
 
## Metadata

- elements {'C', 'H', 'N', 'O', 'S'}
- unique molecules 185
- optimizations 6716
- scf properties:
    - dipole
    - quadrupole
    - wiberg_lowdin_indices
    - mayer_indices
- qc spec
    - name: default
    - method: B3LYP-D3BJ
    - basis: DZVP
    
