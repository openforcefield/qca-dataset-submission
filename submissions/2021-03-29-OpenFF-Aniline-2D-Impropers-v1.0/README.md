# OpenFF Amide Torsion Set v1.0:

### Description

This dataset contains a set of amides which have substituted groups of varying electron donating and withdrawing properties. This dataset was curated in an effort to improve and understand improper torsions in force fields. We will scan the improper and proper angle simultaneously to better understand the coupling and energetics of these torsions.

### General Information

 - Date: 2021.03.29
 - Class: OpenFF TorsionDrive
 - Purpose: Set to amides which contain different electron donating and withdrawing groups
 - Collection: TorsionDriveDataset
 - Name: OpenFF Amide Torsion Set v1.0
 - Number of unique molecules        23
 - Number of filtered molecules      0
 - Number of torsion drives          23
 - Number of conformers min mean max 1 5.96 10
 - Dataset Submitter: Simon Boothroyd
 - Dataset Generator: Jessica Maat
 - Set of charges: [-1.0, 0.0, 1.0]
 - Mean molecular weight: 143.80
 - Max molecular weight: 256.30
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
 'creation_date': datetime.date(2021, 4, 2),
 'dataset_name': 'OpenFF Amide Torsion Set v1.0',
 'elements': {'N', 'C', 'O', 'H', 'S'},
 'long_description': 'This dataset contains a set of amides which have '
                     'substituted groups of varying electron donating and '
                     'withdrawing properties. This dataset was curated in an '
                     'effort to improve and understand improper torsions in '
                     'force fields. We will scan the improper and proper angle '
                     'simultaneously to better understand the coupling and '
                     'energetics of these torsions.',
 'long_description_url': HttpUrl('https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2021-03-29-OpenFF-Aniline-2D-Impropers-v1.0', scheme='https', host='github.com', tld='com', host_type='domain', path='/openforcefield/qca-dataset-submission/tree/master/submissions/2021-03-29-OpenFF-Aniline-2D-Impropers-v1.0'),
 'short_description': 'Substituted amides with various electronc withdrawing '
                      'and donating groups',
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
