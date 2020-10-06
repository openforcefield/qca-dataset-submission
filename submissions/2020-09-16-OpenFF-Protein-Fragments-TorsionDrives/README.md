# OpenFF Protein Fragments TorsionDrives

### Description

This is a protein fragment dataset consisting of torsion drives on various protein fragments prepared by David Cerutti.
We have 12 central residues capped with a combination of different terminal residues.
We drive the following angles for each fragment:
- omega
- phi
- psi
- chi1 (if applicable)
- chi2 (if applicable)

Details from David Cerutti:

> The way it's set up, we are scanning phi and psi of the central residue with a random selection of ALA, GLY, SER, or VAL to the N- or C-terminus,
> and the customary ACE and NME blocking groups outside of that.


### General Information

 - Date: 2020.09.16
 - Class: OpenFF TorsionDrive
 - Purpose: TorsionDrive scans for common amino acids
 - Collection: TorsionDriveDataset 
 - Name: OpenFF Protein Fragments TorsionDrives v1.0
 - Number of unique molecules: 185
 - Submitter: David Dotson
 
Note, each folder contains molecules saved via mol2 in each confirmation; however, the bond order is incorrect.
We let openeye interpret it by re-saving to PDB first before creating the dataset.

### QCSubmit generation pipeline

 - `Dataset Preparation.ipynb`: This notebook shows how the TorsionDrive dataset was prepared from the input files. 
 
### QCSubmit Manifest

- `Dataset Preparation.ipynb`: Dataset creation notebook with details on decisions made for submission.
- `dataset.json.bz2`: The compressed TorsionDrive dataset ready for submission.
- `torsions.pdf`: A pdf file containing molecule 2D structures.
- `Input_files.tar.gz`: Folders containing the input molecule conformations and the (unused) corresponding dihedral restraints.
- `molecules.smi`: SMILES for every molecule in the submission; duplicate molecules for each driven torsion dropped
 
### Metadata

- elements {'C', 'H', 'N', 'O', 'S'}
- unique molecules 185
- scf properties:
    - `dipole`
    - `quadrupole`
    - `wiberg_lowdin_indices`
    - `mayer_indices`
- qc spec
    - name: default
    - method: B3LYP-D3BJ
    - basis: DZVP
