# OpenFF benchmark ligand fragments v2.0

### Description

A torsiondrive dataset created from the [OpenFF FEP benchmark dataset](https://github.com/openmm/openmmforcefields/tree/master/openmmforcefields/data/perses_jacs_systems) ligands. The ligands are fragmented using openeye/ambertools and openff-fragmenter before having key torsions driven.
This second version expands the original dataset with new fragments produced by openff-fragmenter after its refactor as the performance is now slightly different in some cases depending on which backend toolkit is used
and fixes some known issues around incorrect stereochemistry. 

### General Information

- Date: 2021.08.10
- Class: OpenFF TorsionDrive Dataset
- Purpose: JACS ligand fragments torsiondrives for torsion optimization 
- Collection: TorsiondriveDataset
- Name: OpenFF-benchmark-ligand-fragments-v2.0
- Number of unique molecules        490
- Number of filtered molecules      1
- Number of torsiondrives              671
- Number of conformers min mean max 1   2.89  4
- Dataset Submitter: Joshua Horton
- Dataset Generator: Joshua Horton
- Set of charges: [-2.0, -1.0, 0.0, 1.0]
- Mean molecular weight: 259.64
- Max molecular weight: 536.44
- Enumerate stereoisomers: False
- Enumerate tautomers: False
- Enumerate protomers: False

### QCSubmit generation pipeline

- `Dataset_Generation.ipynb`: A notebook which shows how the dataset was prepared from the input files.

### QCSubmit Manifest

- `Dataset_Generation.ipynb`
- `dataset.json.bz2`: The torsiondrive dataset ready for submission.
- `dataset.pdf`: A pdf file containing molecule 2D structures with targeted torsions highlighted.
- `dataset.smi`: SMILES for every molecule in the submission.
 
### Metadata

```
{'collection_type': 'TorsionDriveDataset',
 'creation_date': datetime.date(2021, 8, 23),
 'dataset_name': 'OpenFF-benchmark-ligand-fragments-v2.0',
 'elements': {'S', 'N', 'Br', 'C', 'H', 'O', 'Cl', 'F', 'I'},
 'long_description': 'Ligand fragments generated via openff-fragmenter using '
                     'openeye/ambertools for the JACS benchmark systems. These '
                     'fragments are then used to fit bespoke torsion '
                     'parameters for the bespokefit paper.',
 'long_description_url': HttpUrl('https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2021-08-10-OpenFF-JACS-Fragments-v2.0', scheme='https', host='github.com', tld='com', host_type='domain', path='/openforcefield/qca-dataset-submission/tree/master/submissions/2021-08-10-OpenFF-JACS-Fragments-v2.0'),
 'short_description': 'Ligand fragments from the JACS benchmark systems.',
 'submitter': 'JTHorton'}
```

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

