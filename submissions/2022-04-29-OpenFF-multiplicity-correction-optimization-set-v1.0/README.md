# OpenFF multiplicity correction optimization set v1.0

### Description

An optimization dataset created to correct multiplicity issues in the force field. These are in addition to a corresponding torsiondrive set, which is for training of new torsion parameters that do not have pre-existing data in current training datasets. 

### General Information

- Date: 2022.05.24
- Class: OpenFF Optimization Dataset
- Purpose: An optimization dataset created to correct multiplicity issues in the force field.
- Collection: OptimizationDataset
- Name: OpenFF multiplicity correction optimization set v1.0
- Number of unique molecules        99
- Number of filtered molecules      0
- Number of conformers              400
- Number of conformers min mean max 1   4.04 50
- Mean molecular weight: 172.23
- Max molecular weight: 317.32
- Set of Charges: [-1.0, 0.0, 1.0]
- Elements: 'S', 'P', 'O', 'C', 'H', 'N'
- Dataset Submitter: Jessica Maat
- Dataset Generator: Jessica Maat

### QCSubmit generation pipeline

- `Dataset_Generation.ipynb`: A notebook which shows how the dataset was prepared from the input files.

### QCSubmit Manifest

- `Dataset_Generation.ipynb`: A notebook which shows how the dataset was prepared from the input files.
- `dataset.json.bz2`: The torsiondrive dataset ready for submission.
- `dataset.pdf`: A pdf file containing molecule 2D structures with targeted torsions highlighted.
- `dataset.smi`: SMILES for every molecule in the submission.
 

## QC specifications:

```
Spec: default
{'basis': 'DZVP',
 'implicit_solvent': None,
 'keywords': None,
 'maxiter': 200,
 'method': 'B3LYP-D3BJ',
 'program': 'psi4',
 'scf_properties': ['dipole',
                    'quadrupole',
                    'wiberg_lowdin_indices',
                    'mayer_indices'],
 'spec_description': 'Standard OpenFF optimization quantum chemistry '
                     'specification.',
 'spec_name': 'default',
 'store_wavefunction': 'none'}
```
