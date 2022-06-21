# OpenFF multiplicity correction torsion drive data v1.1

### Description

A torsiondrive dataset created to correct multiplicity issues in the force field. These torsions are for training of the new torsion parameters which do not have pre-existing data in the training datasets. 

### General Information

- Date: 2022.05.24
- Class: OpenFF TorsionDrive Dataset
- Purpose: 
- Collection: TorsiondriveDataset
- Name: OpenFF multiplicity correction torsion drive data v1.1
- Number of unique molecules        131
- Number of filtered molecules      0
- Number of torsion drives          131
- Number of conformers min mean max 1   2.87 10
- Mean molecular weight: 173.20
- Max molecular weight: 400.43
- Set of Charges: [-1.0, 0.0, 1.0]
- Elements: 'S', 'P', 'O', 'C', 'H', 'N'
- Dataset Submitter: Jessica Maat, Pavan Behara
- Dataset Generator: Jessica Maat, Pavan Behara

### Changelog
 Additional molecules added to avoid fitting to rings and also to resolve some convergence errors on v1.0.

### QCSubmit generation pipeline

- `Dataset_Generation.ipynb`: A notebook which shows how the dataset was prepared from the input files.

### QCSubmit Manifest

- `Dataset_Generation.ipynb`: A notebook which shows how the dataset was prepared from the input files.
- `dataset-v1.1.json.bz2`: The torsiondrive dataset ready for submission.
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
