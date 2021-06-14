### Description

This dataset is the public counterpart of the OpenFF Industry Benchmark Season 1.
Each industry partner has selected a range of diverse molecules which represent their current chemical interests.
The dataset will be used in conjunction with private counterparts also designed by each partner to give an unbiased assessment of the progress and current performance of the OpenFF line of force fields in comparison with other contemporary force fields.

This v1.1 dataset features corrected Merck (MRK) molecules with explicit hydrogens.
The original v1.0 dataset did not have explicit hydrogens on these molecules, resulting in poor starting conformers that have largely failed to geometry optimize under QM.
The v1.1 dataset was prepared from the v1.0 dataset, excising the MRK molecules and replacing them with the explicit hydrogen variants prepared using the [Season 1 protocol](https://openforcefield.atlassian.net/wiki/spaces/PS/pages/971898891/Optimization+Benchmarking+Protocol+-+Season+1) via `openff-benchmark`.

### General Information

 - Date: 2021.06.04
 - Class: OpenFF Optimization
 - Purpose: OpenFF Industry Benchmark Season 1 Public
 - Collection: OptimizationDataset
 - Name: OpenFF Industry Benchmark Season 1 v1.1
 - Number of unique molecules        9847
 - Number of filtered molecules      7
 - Number of optimizations           77055
 - Number of conformers min mean max 1   7.77 10
 - Dataset Submitter: David Dotson / Joshua Horton
 - Dataset Generator: David Dotson / Joshua Horton
 - Set of charges: {-2, -1, 0, 1, 2}
 - Mean molecular weight: 348.07
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
 'creation_date': datetime.date(2021, 6, 4),
 'dataset_name': 'OpenFF Industry Benchmark Season 1 v1.1',
 'elements': {'Br', 'F', 'P', 'H', 'N', 'S', 'Cl', 'O', 'C'},
 'long_description': 'The combination of all publicly chosen compound sets by '
                     'industry partners from the OpenFF season 1 industry '
                     'benchmark.',
 'long_description_url': HttpUrl('https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2021-06-04-OpenFF-Industry-Benchmark-Season-1-v1.1', scheme='https', host='github.com', tld='com', host_type='domain', path='/openforcefield/qca-dataset-submission/tree/master/submissions/2021-06-04-OpenFF-Industry-Benchmark-Season-1-v1.1'),
 'short_description': 'The public molecules from the OpenFF Industry '
                      'Benchmark.',
 'submitter': 'dotsdl'}
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
