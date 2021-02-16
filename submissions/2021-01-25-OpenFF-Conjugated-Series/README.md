# OpenFF WBO Conjugated Series v1.0

### Description

This is a torsion drive dataset that probes a range of Wiberg bond orders for different chemistries to better understand the relationship between torsion barrier height and Wiberg bond order. The dataset is being used for developing Wiberg bond order interpolated torsion parameters in OpenFF.

The general principal behind the dataset is to start with a “base” chemical group and substitute chemical groups onto these “base” groups. The “base” groups included in this dataset are enyl, styrene, primary amide, secondary amide, tertiary amide, carbamate, urea and carbonyl group. For each of these “base” chemical groups, we substitute (1) hydroxy, (2) thiol, (3) carboxylic, (4) primary amine, (5) protonated amine, (6) urea , (7) secondary amine, (8) hydroxyl amine , (9) nitrile, (10) alkene, (11)  sulfone, (12) ethoxy, (13) hydroxide groups. The aim is to substitute chemical groups with varying electron withdrawing and donating properties, which will vary the Wiberg bond order of the central torsion bond. This dataset enables exploration of the effects of Wiberg bond order on the torsion barrier height for various chemistries.

This dataset enumerates the protomers, tautomers, and stereoisomers of the molecules.

### General Information

 - Date: 2021.02.21
 - Class: OpenFF TorsionDrive
 - Purpose: A series of functional groups to study bond conjugation effects for FF parameter interpolation
 - Collection: TorsionDriveDataset
 - Name: OpenFF WBO Conjugated Series v1.0
 - Number of unique molecules: 487
 - Number of filtered molecules: 1
 - Number of torsiondrives: 787
 - Number of conformers min mean max: 1   9.05 10
 - Dataset Submitter: Jessica Maat
 - Dataset Generator: Trevor Gokey
 - Set of charges: {-1, 0, 1}
 - Mean molecular weight: 115.52
 - Max molecular weight: 182.24
 - Enumerate stereoisomers: True
 - Enumerate tautomers: True
 - Enumerate protomers: True


### QCSubmit generation pipeline

 - `Generate_dataset`: This notebook shows how the TorsionDrive dataset was prepared. All data is contained in the notebook.

### QCSubmit Manifest

- `Generate_dataset.ipynb`: Dataset creation notebook with details on decisions made for submission.
- `dataset.json.bz2`: The compressed TorsionDrive dataset ready for submission.
- `dataset.pdf`: A pdf file containing molecule 2D structures.
- `dataset.smi`: SMILES for every molecule in the submission; duplicate molecules for each driven torsion dropped

### Metadata

```
{'collection_type': 'TorsiondriveDataset',
 'creation_date': datetime.date(2021, 2, 4),
 'dataset_name': 'OpenFF WBO Conjugated Series v1.0',
 'elements': {'N', 'O', 'C', 'H', 'S'},
 'long_description': 'This is a torsion drive dataset that probes a range of '
                     'Wiberg bond orders for different chemistries to better '
                     'understand the relationship between torsion barrier '
                     'height and Wiberg bond order. The dataset is being used '
                     'for developing Wiberg bond order interpolated torsion '
                     'parameters in OpenFF.\n'
                     'The general principal behind the dataset is to start '
                     'with a “base” chemical group and substitute chemical '
                     'groups onto these “base” groups. The “base” groups '
                     'included in this dataset are enyl, styrene, primary '
                     'amide, secondary amide, tertiary amide, carbamate, urea '
                     'and carbonyl group. For each of these “base” chemical '
                     'groups, we substitute (1) hydroxy, (2) thiol, (3) '
                     'carboxylic, (4) primary amine, (5) pronated amine, (6) '
                     'urea , (7) secondary amine, (8) hydroxyl amine , (9) '
                     'nitrile, (10) alkene, (11)  sulfone, (12) ethoxy, (13) '
                     'hydroxide groups. The aim is to substitute chemical '
                     'groups with varying electron withdrawing and donating '
                     'properties, which will vary the Wiberg bond order of the '
                     'central torsion bond. This dataset enables exploration '
                     'of the effects of Wiberg bond order on the torsion '
                     'barrier height for various chemistries.\n'
                     'This dataset enumerates the protomers, tautomers, and '
                     'stereoisomers of the molecules.',
 'long_description_url': HttpUrl('https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/OpenFF-WBO-Conjugated-Series', scheme='https', host='github.com', tld='com', host_type='domain', path='/openforcefield/qca-dataset-submission/tree/master/submissions/OpenFF-WBO-Conjugated-Series'),
 'short_description': 'A series of functional groups to study bond conjugation '
                      'effects for FF parameter interpolation',
 'submitter': 'jmaat'}
```


## QC specifications:

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

## SCF Properties

```
[<SCFProperties.Dipole: 'dipole'>,
 <SCFProperties.Quadrupole: 'quadrupole'>,
 <SCFProperties.WibergLowdinIndices: 'wiberg_lowdin_indices'>,
 <SCFProperties.MayerIndices: 'mayer_indices'>]
 ```
