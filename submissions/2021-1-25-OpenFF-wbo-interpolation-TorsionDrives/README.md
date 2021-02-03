# OpenFF WBO TorsionDrives

### Description

This is a torsion drive dataset that probes a range of wiberg bond order for different chemitries. The dataset is being used for developing wiberg bond order interpolated torsion parameters in OpenFF.

The general principal behind this dataset is we start with a base chemical group enyl, styrene, primary amide, secondary amide, tertiary amide, carbamate, urea and carboynl group. For each of these chemical groups, we substitute (1) hydroxy, (2) thiol, (3) carboxylic, (4) primary amine, (5) pronated amine, (6) urea , (7) secondary amine, (8) hydroxyl amine , (9) nitrile, (10) alkene, (11)  sulfone, (12) ethoxy, (13) hydroxide group. The aim to is to substitute chemical groups with varying electron withdrawing and donating properties, which will vary the wiberg bond order of the central torsion bond. This dataset enables exploration of the affects of wiberg bond order on the torsion barrier height for various chemistries.


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
