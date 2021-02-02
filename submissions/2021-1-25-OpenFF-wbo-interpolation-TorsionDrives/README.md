# OpenFF WBO TorsionDrives

### Description

This is a torsion drive dataset that probes a range of WBO for different chemitries. The dataset is being used for developing WBO interpolated torsion parameters in OpenFF


### General Information

 - Date: 2021.02.21
 - Class: OpenFF TorsionDrive
 - Purpose: TorsionDrive scans molecules with varying WBO and subsituted functional groups
 - Collection: TorsionDriveDataset
 - Name: OpenFF WBO TorsionDrives v1.0
 - Number of unique molecules: 94
 - Number of torsiondrives: 96
 - Submitter: Jessica Maat


### QCSubmit generation pipeline

 - `Generate_dataset`: This notebook shows how the TorsionDrive dataset was prepared from the input files.

### QCSubmit Manifest

- `Generate_dataset.ipynb`: Dataset creation notebook with details on decisions made for submission.
- `dataset.json.bz2`: The compressed TorsionDrive dataset ready for submission.
- `torsions_wbodataset.pdf`: A pdf file containing molecule 2D structures.
- `molecules_wbo.smi`: SMILES for every molecule in the submission; duplicate molecules for each driven torsion dropped

### Metadata

- unique molecules: 94
- torsiondrives: 96
- scf properties:
    - `dipole`
    - `quadrupole`
    - `wiberg_lowdin_indices`
    - `mayer_indices`
- qc spec
    - name: default
    - method: B3LYP-D3BJ
    - basis: DZVP
