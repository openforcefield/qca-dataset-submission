# OpenFF BCC Refit Study COH v2.0

### Description

This dataset contains references to the quantum chemistry (QC) data used to train the OpenFF SMIRNOFF 
[Sage 2.0.0](https://github.com/openforcefield/openff-sage) forcefield. 

### General Information

<!-- copy outputs from generation notebook into these fields -->

- Date: <date>
- Class: OpenFF Optimization Dataset
- Purpose: <brief description>
- Collection: OptimizationDataset
- Name: <dataset name>
- Number of unique molecules        <...>
- Number of filtered molecules     <...> 
- Number of conformers             <...> 
- Number of conformers min mean max <...>   <...>   <...>
- Dataset Submitter: <name of submitter>
- Dataset Generator: <name of generator>
- Set of charges: <[...,....]>
- Mean molecular weight: <...>
- Max molecular weight: <...>
- Enumerate stereoisomers: <bool>
- Enumerate tautomers: <bool>
- Enumerate protomers: <bool>

### QCSubmit generation pipeline

- `generate-dataset.ipynb`: A notebook which shows how the dataset was prepared from the input files.

### QCSubmit Manifest

- `generate-dataset.ipynb`
- `dataset.json.bz2`: The basic dataset ready for submission.
- `dataset.pdf`: A pdf file containing molecule 2D structures.
- `dataset.smi`: SMILES for every molecule in the submission.
 
### Metadata

<!-- copy metadata outputs from generation notebook here -->

```
{'collection_type': 'OptimizationDataset',
 'creation_date': datetime.date(2021, 7, 2),
 'dataset_name': 'OpenFF BCC Refit Study COH v2.0',
 'elements': {'H', 'O', 'C'},
 'long_description': 'A data set curated for the initial stage of the on-going '
                     'OpenFF study which aims to co-optimize the AM1BCC bond '
                     'charge correction (BCC) parameters against an '
                     'experimental training set of density and enthalpy of '
                     'mixing data points and a QM training set of electric '
                     'field data.\n'
                     '\n'
                     'The initial data set is limited to only molecules '
                     'composed of C, O, H. This limited scope significantly '
                     'reduces the number of BCC parameters which must be '
                     'retrained, thus allowing for easier convergence of the '
                     'initial optimizations.\n'
                     '\n'
                     'The included molecules were combinatorially generated to '
                     'cover a range of alcohol, ether, and carbonyl containing '
                     'molecules.',
 'long_description_url': HttpUrl('https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2021-06-22-OpenFF-BCC-Refit-Study-COH-v2.0', scheme='https', host='github.com', tld='com', host_type='domain', path='/openforcefield/qca-dataset-submission/tree/master/submissions/2021-06-22-OpenFF-BCC-Refit-Study-COH-v2.0'),
 'short_description': 'Optimizations of diverse, para-substituted aniline '
                      'derivatives.',
 'submitter': 'simonboothroyd'}
```

## QC specifications:

<!-- copy spec outputs from generation notebook here -->

```

```

## SCF Properties

<!-- copy scf property outputs from generation notebook here -->

```
[<SCFProperties.Dipole: 'dipole'>,
 <SCFProperties.Quadrupole: 'quadrupole'>,
 <SCFProperties.WibergLowdinIndices: 'wiberg_lowdin_indices'>,
 <SCFProperties.MayerIndices: 'mayer_indices'>]
```
