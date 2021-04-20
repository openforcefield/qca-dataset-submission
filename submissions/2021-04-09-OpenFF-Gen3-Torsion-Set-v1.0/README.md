# OpenFF Gen3 Torsion Set v1.0

### Description

This dataset is a simple-molecule-only dataset, which is a part of Sage torsion parameter training set candidate. 

### General Information

 - Date: 2021.04.09
 - Class: OpenFF TorsionDrive
 - Purpose: For valence parameter optimization 
 - Collection: TorsionDriveDataset
 - Name: OpenFF Gen3 Torsion Set v1.0
 - Number of unique molecules        888
 - Number of filtered molecules      0
 - Number of torsion drives          889
 - Number of conformers min mean max 1   2.61 12
 - Dataset Submitter: Hyesu Jang
 - Dataset Generator: Hyesu Jang
 - Set of charges: [0.0, 1.0]
 - Mean molecular weight: 131.31
 - Max molecular weight: 433.68
 - Enumerate stereoisomers: False
 - Enumerate tautomers: False
 - Enumerate protomers: False 

### QCSubmit generation pipeline

 - `generate-dataset`: The notebook used to generate the dataset from the input smiles set.

### QCSubmit Manifest
- `1.gen-substituent-lists`: Step1. Generation of a list of small and rigid substituents.
- `2.gen-mol-list`: Step2. Generation of small molecules.
- `3.filter-mol-list`: Step3. Elimination of internal hydrogen bond forming molecules.
- `4.select-torsions`: Step4. List molecules matching to each torsion parameter and select torsions.
- `generate-dataset`: Dataset creation notebook with details on decisions made for submission.
- `dataset.json.bz2`: The compressed TorsionDrive dataset ready for submission.
- `dataset.pdf`: A pdf file containing molecule 2D structures.
- `dataset.smi`: SMILES for every molecule in the submission; duplicate molecules for each driven torsion dropped

### Metadata

```
{'basis': 'DZVP',
 'implicit_solvent': None,
 'method': 'B3LYP-D3BJ',
 'program': 'psi4',
 'spec_description': 'Standard OpenFF optimization quantum chemistry '
                     'specification.',
 'spec_name': 'default',
 'store_wavefunction': 'none'}
```

### SCF Properties

```
[<SCFProperties.Dipole: 'dipole'>,
 <SCFProperties.Quadrupole: 'quadrupole'>,
 <SCFProperties.WibergLowdinIndices: 'wiberg_lowdin_indices'>,
 <SCFProperties.MayerIndices: 'mayer_indices'>]
```
