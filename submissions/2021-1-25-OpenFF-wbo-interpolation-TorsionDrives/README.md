# OpenFF WBO Conjugated Torsion Series v1.0

### Description

This is a torsion drive dataset that probes a range of Wiberg bond orders for different chemistries to better understand the relationship between torsion barrier height and Wiberg bond order. The dataset is being used for developing Wiberg bond order interpolated torsion parameters in OpenFF.

The general principal behind the dataset is to start with a “base” chemical group and substitute chemical groups onto these “base” groups. The “base” groups included in this dataset are enyl, styrene, primary amide, secondary amide, tertiary amide, carbamate, urea and carboynl group. For each of these “base” chemical groups, we substitute (1) hydroxy, (2) thiol, (3) carboxylic, (4) primary amine, (5) pronated amine, (6) urea , (7) secondary amine, (8) hydroxyl amine , (9) nitrile, (10) alkene, (11)  sulfone, (12) ethoxy, (13) hydroxide groups. The aim is to substitute chemical groups with varying electron withdrawing and donating properties, which will vary the Wiberg bond order of the central torsion bond. This dataset enables exploration of the effects of Wiberg bond order on the torsion barrier height for various chemistries.

### General Information

 - Date: 2021.02.21
 - Class: OpenFF TorsionDrive
 - Purpose: TorsionDrive scans molecules with varying WBO and subsituted functional groups
 - Collection: TorsionDriveDataset
 - Name: OpenFF WBO Conjugated Torsion Series v1.0
 - Number of unique molecules: 104
 - Number of torsiondrives: 104
 - Submitter: Jessica Maat


### QCSubmit generation pipeline

 - `Generate_dataset`: This notebook shows how the TorsionDrive dataset was prepared from the input files.

### QCSubmit Manifest

- `Generate_dataset.ipynb`: Dataset creation notebook with details on decisions made for submission.
- `dataset.json.bz2`: The compressed TorsionDrive dataset ready for submission.
- `torsions_wbodataset.pdf`: A pdf file containing molecule 2D structures.
- `molecules_wbo.smi`: SMILES for every molecule in the submission; duplicate molecules for each driven torsion dropped

### Metadata

- unique molecules: 104
- torsiondrives: 104
- scf properties:
    - `dipole`
    - `quadrupole`
    - `wiberg_lowdin_indices`
    - `mayer_indices`
- qc spec
    - name: default
    - method: B3LYP-D3BJ
    - basis: DZVP
