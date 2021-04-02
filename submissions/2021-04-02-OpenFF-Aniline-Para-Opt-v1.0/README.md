# OpenFF Aniline Para Opt v1.0

### Description

This dataset contains a set of aniline derivatives which are para-substituted with groups of varying electron donating 
and withdrawing properties. This dataset was curated in an effort to improve and understand which functional groups 
will allow us to best interpolate between a planer and a pyramidal nitrogen.

### General Information

 - Date: 2021.04.02
 - Class: OpenFF TorsionDrive
 - Purpose: Optimizations of diverse, para-substituted aniline derivatives.
 - Collection: OptimizationDataset
 - Name: OpenFF Aniline Para Opt v1.0
 - Number of unique molecules        50
 - Number of filtered molecules      0
 - Number of conformers              223
 - Number of conformers min mean max 1   4.46 10
 - Dataset Submitter: Simon Boothroyd
 - Dataset Generator: Jessica Maat
 - Set of charges: [-1.0, 0.0, 1.0]
 - Mean molecular weight: 150.42
 - Max molecular weight: 343.84
 - Enumerate stereoisomers: False
 - Enumerate tautomers: False
 - Enumerate protomers: False


### QCSubmit generation pipeline

- `generate-smiles`: A python script which combines an aniline scaffold with diverse functional groups.
- `generate-dataset`: The notebook used to generate the dataset from the input smiles set.

### QCSubmit Manifest

- `generate-smiles.py`: Molecule generation script.
- `generate-dataset`: Dataset creation notebook with details on decisions made for submission.
- `dataset.json.bz2`: The compressed dataset ready for submission.
- `dataset.pdf`: A pdf file containing molecule 2D structures.
- `dataset.smi`: SMILES for every molecule in the submission.

### Metadata

```
{'collection_type': 'OptimizationDataset',
 'creation_date': datetime.date(2021, 4, 2),
 'dataset_name': 'OpenFF Aniline Para Opt v1.0',
 'elements': {'Br', 'C', 'O', 'N', 'S', 'H', 'Cl', 'F'},
 'long_description': 'This dataset contains a set of aniline derivatives which '
                     'are para-substituted with groups of varying electron '
                     'donating and withdrawing properties. This dataset was '
                     'curated in an effort to improve and understand which '
                     'functional groups will allow us to best interpolate '
                     'between a planer and a pyramidal nitrogen.',
 'long_description_url': HttpUrl('https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2021-04-02-OpenFF-Aniline-Para-Opt-v1.0', scheme='https', host='github.com', tld='com', host_type='domain', path='/openforcefield/qca-dataset-submission/tree/master/submissions/2021-04-02-OpenFF-Aniline-Para-Opt-v1.0'),
 'short_description': 'Optimizations of diverse, para-substituted aniline '
                      'derivatives.',
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
