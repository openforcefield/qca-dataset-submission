# OpenFF Gen3 Torsion Set v1.0

### Description

This dataset is a simple-molecule-only dataset, a candidate of Sage torsion parameter training set. 

### General Information

 - Date: 2021.04.09
 - Class: OpenFF TorsionDrive
 - Purpose: For valence parameter optimization 
 - Collection: TorsionDriveDataset
 - Name: OpenFF Gen3 Torsion Set v1.0
 - Number of unique molecules        2433
 - Number of filtered molecules      167
 - Number of torsion drives          4684
 - Number of conformers min mean max 1   1.05 4
 - Dataset Submitter: Hyesu Jang
 - Dataset Generator: Hyesu Jang
 - Set of charges: [-1.0, 0.0, 1.0, 2.0]
 - Mean molecular weight: 154.34
 - Max molecular weight: 514.56
 - Enumerate stereoisomers: True
 - Enumerate tautomers: False
 - Enumerate protomers: False 

### QCSubmit generation pipeline

 - `generate-dataset`: The notebook used to generate the dataset from the input smiles set.

### QCSubmit Manifest

- `generate-dataset`: Dataset creation notebook with details on decisions made for submission.
- `dataset.json.bz2`: The compressed TorsionDrive dataset ready for submission.
- `dataset.pdf`: A pdf file containing molecule 2D structures.
- `dataset.smi`: SMILES for every molecule in the submission; duplicate molecules for each driven torsion dropped

### Metadata
