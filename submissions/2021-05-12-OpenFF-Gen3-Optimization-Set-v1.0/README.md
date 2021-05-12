# OpenFF Gen3 Optimization Set v1.0

### Description

This dataset is a simple-molecule-only optimization dataset, which is a part of Sage torsion parameter training set candidate. 

### General Information

 - Date: 2021.05.12
 - Class: OpenFF OptimizationDataset
 - Purpose: For valence parameter optimization 
 - Collection: OptimizationDataset
 - Name: OpenFF Gen3 Optimization Set v1.0
 - Number of unique molecules        3312
 - Number of filtered molecules      237
 - Number of records                 29656
 - Number of conformers min mean max 1   8.95 10
 - Dataset Submitter: Hyesu Jang
 - Dataset Generator: Hyesu Jang
 - Set of charges: [-2.0, -1.0, 0.0, 1.0, 2.0]
 - Mean molecular weight: 132.29
 - Max molecular weight: 433.68
 - Enumerate stereoisomers: True
 - Enumerate tautomers: True
 - Enumerate protomers: True 

### QCSubmit generation pipeline

 - `generate-dataset`: The notebook used to generate the dataset from the input smiles set.

### QCSubmit Manifest
- `generate-dataset`: Dataset creation notebook with details on decisions made for submission.
- `dataset.json.bz2`: The compressed Optimization dataset ready for submission.
- `dataset.pdf`: A pdf file containing molecule 2D structures.
- `dataset.smi`: SMILES for every molecule in the submission; duplicate molecules for each driven torsion dropped

### Metadata
```
{'collection_type': 'OptimizationDataset',
 'creation_date': datetime.date(2021, 5, 12),
 'dataset_name': 'OpenFF Gen3 Optimization Set v1.0',
 'elements': {'F', 'P', 'N', 'S', 'Br', 'H', 'Cl', 'O', 'C'},
 'long_description': 'This dataset is a simple-molecule-only optimization '
                     'dataset. The input molecules are those being scanned in '
                     'OpenFF Gen3 Torsion Set',
 'long_description_url': HttpUrl('https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2021-05-12-OpenFF-Gen3-Optimization-Set-v1.0', scheme='https', host='github.com', tld='com', host_type='domain', path='/openforcefield/qca-dataset-submission/tree/master/submissions/2021-05-12-OpenFF-Gen3-Optimization-Set-v1.0'),
 'short_description': 'OpenFF Gen3 Torsion Set v1.0',
 'submitter': 'hyesujang'}
 ```
 
### QCSpecifications
```
Spec: default
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
