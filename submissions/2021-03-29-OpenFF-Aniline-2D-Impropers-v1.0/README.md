# OpenFF Amide Torsion Set v1.0

### Description

TODO TODO

### General Information

 - Date: 2021.03.29
 - Class: OpenFF TorsionDrive
 - Purpose: TODO TODO
 - Collection: TorsionDriveDataset
 - Name: OpenFF Amide Torsion Set v1.0
 - Number of unique molecules        XXX
 - Number of filtered molecules      XXX
 - Number of torsion drives          XXX
 - Number of conformers min mean max XXX
 - Dataset Submitter: Simon Boothroyd
 - Dataset Generator: Jessica Maat
 - Set of charges: XXX
 - Mean molecular weight: XXX
 - Max molecular weight: XXX
 - Enumerate stereoisomers: False
 - Enumerate tautomers: False
 - Enumerate protomers: False


### QCSubmit generation pipeline

 - `generate-dataset`: The notebook used to generate the dataset from the input smiles set.

### QCSubmit Manifest

- `generate-smiles.py`: Molecule generation script.
- `generate-dataset`: Dataset creation notebook with details on decisions made for submission.
- `dataset.json.bz2`: The compressed TorsionDrive dataset ready for submission.
- `dataset.pdf`: A pdf file containing molecule 2D structures.
- `dataset.smi`: SMILES for every molecule in the submission; duplicate molecules for each driven torsion dropped

### Metadata

```
{'collection_type': 'TorsiondriveDataset',
 'creation_date': datetime.date(2021, 3, 24),
 'dataset_name': 'OpenFF Amide Torsion Set v1.0',
 'elements': {'H', 'O', 'S', 'C', 'N'},
 'long_description': "This dataset contains a set of 'simple' amides, "
                     'thioamides, and amidines which have been functionalised '
                     'with groups of varying electron-donating and '
                     'electron-withdrawing characteristics. It was curated as '
                     'part of an effort to explore why the cis-trans '
                     'preference of amides is poorly reproduced in OpenFF '
                     '1.3.0, and as an extra source of data to try and remedy '
                     'this.',
 'long_description_url': HttpUrl('https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2021-03-23-OpenFF-Amide-Torsion-Set-v1.0', scheme='https', host='github.com', tld='com', host_type='domain', path='/openforcefield/qca-dataset-submission/tree/master/submissions/2021-03-23-OpenFF-Amide-Torsion-Set-v1.0'),
 'short_description': 'Amides, thioamides and amidines diversely '
                      'functionalized.',
 'submitter': 'simonboothroyd'}
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
