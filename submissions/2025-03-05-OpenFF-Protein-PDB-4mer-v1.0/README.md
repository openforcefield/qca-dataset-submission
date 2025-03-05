# OpenFF Protein PDB 4-mers V1.0

### Description

This dataset is composed of 1,000 4-mer structures which were extracted from PDB entries. The purpose of this dataset is to fill in gaps within the existing protein training data with secondary structures which occur in real proteins. For all scripts used to generated the input files please visit https://github.com/ajfriedman22/4mer_Generation.

### General Information

<!-- copy outputs from generation notebook into these fields -->

- Date: 2025-03-05
- Class: OpenFF Optimization Dataset
- Purpose:Optimizations of 4-mers extracted from the PDB
- Collection: OptimizationDataset
- Name: 2025-03-05-OpenFF-Protein-PDB-4mer-v1.0
- Number of unique molecules        200
- Number of filtered molecules     0 
- Number of conformers             1000 
- Number of conformers min mean max 5 5 5
- Dataset Submitter: Anika Friedman
- Dataset Generator: Anika Friedman
- Set of charges: [-2.0, -1.0, 0.0, 1.0]
- Mean molecular weight: 451.35
- Max molecular weight: 570.64

### QCSubmit generation pipeline

* `generate-dataset.ipynb`: This notebooks shows how the dataset was prepared from the
  input files in `inputs`.

### QCSubmit Manifest

### Input files
* `inputs/*.sdf`: Input SMILES strings for dataset molecules
* `generate-dataset.ipynb`: Notebook describing dataset generation and submission
* `input-env.yaml`: Environment file used to create Python environment for the script

### Output Files
* `dataset.json.bz2`: Compressed dataset ready for submission
* `dataset.pdf`: Visualization of dataset molecules
* `output.smi`: SMILES strings for dataset molecules
 
### Metadata

<!-- copy metadata outputs from generation notebook here -->

```
{'collection_type': 'OptimizationDataset',
 'creation_date': datetime.date(2025, 3, 5),
 'dataset_name': 'OpenFF Protein PDB 4-mers v1.0',
 'elements': {'O', 'C', 'N', 'H'},
 'long_description': 'This dataset aims to eliminate some of the present '
                     'issues in the OpenFF protein force field development '
                     'process. It has proven difficult  to develop a '
                     'self-consistant force field which performs well for '
                     'small molecules and also maintains stable secondary '
                     'structer in proteins. By adding extracted 4-mersfrom '
                     'deposited PDB structures we aim to introduce more '
                     'physically realistic secondarystructure to the QC '
                     'training data to better locate minima within dihedral '
                     'space forprotein residue.\n'
                     '\n'
                     'This dataset includes 200 distinct sequences with 5 '
                     'conformers each extracted from the Top8000 PDB database.',
 'long_description_url': HttpUrl('https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2025-03-05-OpenFF-Protein-PDB-4mer-v1.0', ),
 'short_description': 'Optimizations of 4-mers extracted from the PDB.',
 'submitter': 'anikafriedman'}
```

## QC specifications:

```
{'basis': 'DZVP',
 'implicit_solvent': None,
 'keywords': {},
 'maxiter': 200,
 'method': 'B3LYP-D3BJ',
 'program': 'psi4',
 'spec_description': 'Standard OpenFF optimization quantum chemistry '
                     'specification.',
 'spec_name': 'default',
 'store_wavefunction': 'none'}
```

## SCF Properties
```
['dipole',
                    'quadrupole',
                    'wiberg_lowdin_indices',
                    'mayer_indices']
```