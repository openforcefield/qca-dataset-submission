### Description

This dataset is the MM component of the "OpenFF Industry Benchmark Season 1 v1.1".
It uses the final geometries of the QM results (computed with `b3lyp-d3bj/dzvp`) from that set as the starting conformation for all MM optimizations.

Each industry partner has selected a range of diverse molecules which represent their current chemical interests.
The dataset will be used in conjunction with private counterparts also designed by each partner to give an unbiased assessment of the progress and current performance of the OpenFF line of force fields in comparison with other contemporary force fields.

### General Information

 - Date: 2021.07.22
 - Class: OpenFF Optimization
 - Purpose: OpenFF Industry Benchmark Season 1 Public
 - Collection: OptimizationDataset
 - Name: OpenFF Industry Benchmark Season 1 - MM v1.1
 - Number of unique molecules        9467
 - Number of optimizations           71655
 - Number of conformers min mean max 1   1.00   1
 - Dataset Submitter: David Dotson / Joshua Horton
 - Dataset Generator: David Dotson / Joshua Horton
 - Set of charges: {-2, -1, 0, 1, 2}
 - Mean molecular weight: 358.69
 - Max molecular weight: 1104.40


### QCSubmit generation pipeline

 - `Industry Dataset - MM.ipynb`: preparation notebook using [openff-benchmark](https://github.com/openforcefield/openff-benchmark).

### QCSubmit Manifest

 - `Industry Dataset - MM.ipynb`: preparation notebook
 - `dataset.json.bz2`: compressed Optimization dataset ready for submission
 - `dataset.smi`: SMILES for every molecule in the submission
 
### Metadata

```
{'submitter': 'dotsdl',
 'creation_date': datetime.date(2021, 7, 22),
 'collection_type': 'OptimizationDataset',
 'dataset_name': 'OpenFF Industry Benchmark Season 1 - MM v1.1',
 'short_description': 'The public molecules from the OpenFF Industry Benchmark, MM optimizations following QM.',
 'long_description_url': HttpUrl('https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2021-07-28-OpenFF-Industry-Benchmark-Season-1-MM-v1.1', scheme='https', host='github.com', tld='com', host_type='domain', path='/openforcefield/qca-dataset-submission/tree/master/submissions/2021-07-28-OpenFF-Industry-Benchmark-Season-1-MM-v1.1'),
 'long_description': 'MM combination of all publicly chosen compound sets by industry partners from the OpenFF season 1 industry benchmark.',
 'elements': {'Br', 'C', 'Cl', 'F', 'H', 'N', 'O', 'P', 'S'}}

```

### QC specifications

 
```
{'smirnoff99Frosst-1.1.0': {'method': 'smirnoff99frosst-1.1.0',
  'basis': 'smirnoff',
  'program': 'openmm',
  'spec_name': 'smirnoff99Frosst-1.1.0',
  'spec_description': 'default smirnoff99Frosst-1.1.0 optimization spec',
  'store_wavefunction': 'none',
  'implicit_solvent': None},
 'openff-1.0.0': {'method': 'openff-1.0.0',
  'basis': 'smirnoff',
  'program': 'openmm',
  'spec_name': 'openff-1.0.0',
  'spec_description': 'default openff-1.0.0 optimization spec',
  'store_wavefunction': 'none',
  'implicit_solvent': None},
 'openff-1.1.1': {'method': 'openff-1.1.1',
  'basis': 'smirnoff',
  'program': 'openmm',
  'spec_name': 'openff-1.1.1',
  'spec_description': 'default openff-1.1.1 optimization spec',
  'store_wavefunction': 'none',
  'implicit_solvent': None},
 'openff-1.2.1': {'method': 'openff-1.2.1',
  'basis': 'smirnoff',
  'program': 'openmm',
  'spec_name': 'openff-1.2.1',
  'spec_description': 'default openff-1.2.1 optimization spec',
  'store_wavefunction': 'none',
  'implicit_solvent': None},
 'openff-1.3.0': {'method': 'openff-1.3.0',
  'basis': 'smirnoff',
  'program': 'openmm',
  'spec_name': 'openff-1.3.0',
  'spec_description': 'default openff-1.3.0 optimization spec',
  'store_wavefunction': 'none',
  'implicit_solvent': None},
 'openff-2.0.0-rc.2': {
  'method': 'openff-2.0.0-rc.2',
  'basis': 'smirnoff',
  'program': 'openmm',
  'spec_name': 'openff-2.0.0-rc.2',
  'spec_description': 'default openff-2.0.0-rc.2 optimization spec',
  'store_wavefunction': 'none',
  'implicit_solvent': None},
 'gaff-2.11': {'method': 'gaff-2.11',
  'basis': 'antechamber',
  'program': 'openmm',
  'spec_name': 'gaff-2.11',
  'spec_description': 'default gaff-2.11 optimization spec',
  'store_wavefunction': 'none',
  'implicit_solvent': None}}
```
