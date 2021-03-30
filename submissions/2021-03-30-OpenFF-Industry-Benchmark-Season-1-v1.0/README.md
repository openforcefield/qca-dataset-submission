### Description

This dataset is the public counterpart of the OpenFF Industry Benchmark Season 1. Each industry partner has selected a 
range of diverse molecules which represent there current chemical interests. The dataset will be used in conjunction with 
private counterparts also designed by each partner to give an unbiased assessment of the progress and current performance 
of the OpenFF line of force fields in comparison with other contemporary force fields. 


### General Information

 - Date: 2021.03.30
 - Class: OpenFF Optimization
 - Purpose: OpenFF Industry Benchmark Season 1 Public
 - Collection: OptimizationDataset
 - Name: OpenFF Industry Benchmark Season 1 v1.0
 - Number of unique molecules        9068
 - Number of filtered molecules      17
 - Number of optimizations           69672
 - Number of conformers min mean max 1   7.65 10
 - Dataset Submitter: David Dotson / Joshua Horton
 - Dataset Generator: David Dotson / Joshua Horton
 - Set of charges: {-2, -1, 0, 1, 2}
 - Mean molecular weight: 341.97
 - Max molecular weight: 1104.40
 - Enumerate stereoisomers: False
 - Enumerate tautomers: False
 - Enumerate protomers: False


### QCSubmit generation pipeline
 - `generate-dataset`: The notebook used to generate the dataset from the input serialised openff-qcsubmit datasets produced using [openff-benchmark](https://github.com/openforcefield/openff-benchmark).

### QCSubmit Manifest

- `generate-dataset`: Dataset creation notebook with details on molecule filtering. 
- `dataset.json.bz2`: The compressed Optimization dataset ready for submission.
- `dataset.pdf`: A pdf file containing molecule 2D structures.
- `dataset.smi`: SMILES for every molecule in the submission; duplicate molecules for each optimization dropped

### Metadata

```
{'collection_type': 'OptimizationDataset',
 'creation_date': datetime.date(2021, 3, 30),
 'dataset_name': 'OpenFF Industry Benchmark Season 1 v1.0',
 'elements': {'N', 'F', 'Cl', 'C', 'H', 'O', 'Br', 'P', 'S'},
 'long_description': 'The combination of all publicly chossen compound sets by '
                     'industry partners from the OpenFF season 1 industry '
                     'benchmark.',
 'long_description_url': HttpUrl('https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2021-03-30-OpenFF-Industry-Benchmark-Season-1-v1.0', scheme='https', host='github.com', tld='com', host_type='domain', path='/openforcefield/qca-dataset-submission/tree/master/submissions/2021-03-30-OpenFF-Industry-Benchmark-Season-1-v1.0'),
 'short_description': 'The public molecules from the OpenFF Industry '
                      'Benchmark.',
 'submitter': 'jthorton'}
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
