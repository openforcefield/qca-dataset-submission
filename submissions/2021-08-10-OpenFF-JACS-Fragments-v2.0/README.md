# OpenFF benchmark ligand fragments v2.0

### Description

A torsiondrive dataset created from the [OpenFF FEP benchmark dataset](https://github.com/openmm/openmmforcefields/tree/master/openmmforcefields/data/perses_jacs_systems) ligands.
The ligands are fragmented using openeye/ambertools and openff-fragmenter before having key torsions driven. 
Multiple specifications are included alongside the openff default to be used for benchmarking.

This second version expands the original dataset with new fragments produced by openff-fragmenter after its refactor,
as the performance is now slightly different in some cases depending on which backend toolkit is used,
and fixes some known issues around incorrect stereochemistry.
Here we combine fragments from both antechamber and openeye backends into one dataset.

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

- `Dataset_Generation.ipynb`: A notebook which shows how the dataset was prepared from the input files.
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
Spec: gfn0xtb
{'basis': None,
 'implicit_solvent': None,
 'keywords': None,
 'maxiter': 200,
 'method': 'gfn0xtb',
 'program': 'xtb',
 'scf_properties': ['dipole',
                    'quadrupole',
                    'wiberg_lowdin_indices',
                    'mayer_indices'],
 'spec_description': 'A default spec for gn0xtb',
 'spec_name': 'gfn0xtb',
 'store_wavefunction': 'none'}
Spec: gfn1xtb
{'basis': None,
 'implicit_solvent': None,
 'keywords': None,
 'maxiter': 200,
 'method': 'gfn1xtb',
 'program': 'xtb',
 'scf_properties': ['dipole',
                    'quadrupole',
                    'wiberg_lowdin_indices',
                    'mayer_indices'],
 'spec_description': 'A default spec for gfn1xtb',
 'spec_name': 'gfn1xtb',
 'store_wavefunction': 'none'}
Spec: gfn2xtb
{'basis': None,
 'implicit_solvent': None,
 'keywords': None,
 'maxiter': 200,
 'method': 'gfn2xtb',
 'program': 'xtb',
 'scf_properties': ['dipole',
                    'quadrupole',
                    'wiberg_lowdin_indices',
                    'mayer_indices'],
 'spec_description': 'A default spec for gfn2xtb',
 'spec_name': 'gfn2xtb',
 'store_wavefunction': 'none'}
Spec: gfnff
{'basis': None,
 'implicit_solvent': None,
 'keywords': None,
 'maxiter': 200,
 'method': 'gfnff',
 'program': 'xtb',
 'scf_properties': ['dipole',
                    'quadrupole',
                    'wiberg_lowdin_indices',
                    'mayer_indices'],
 'spec_description': 'A default spec for gfnff',
 'spec_name': 'gfnff',
 'store_wavefunction': 'none'}
Spec: ani2x
{'basis': None,
 'implicit_solvent': None,
 'keywords': None,
 'maxiter': 200,
 'method': 'ani2x',
 'program': 'torchani',
 'scf_properties': ['dipole',
                    'quadrupole',
                    'wiberg_lowdin_indices',
                    'mayer_indices'],
 'spec_description': 'A default spec for ani2x',
 'spec_name': 'ani2x',
 'store_wavefunction': 'none'}
Spec: openff-1.0.0
{'basis': 'smirnoff',
 'implicit_solvent': None,
 'keywords': None,
 'maxiter': 200,
 'method': 'openff-1.0.0',
 'program': 'openmm',
 'scf_properties': ['dipole',
                    'quadrupole',
                    'wiberg_lowdin_indices',
                    'mayer_indices'],
 'spec_description': 'A default spec for openff-1.0.0',
 'spec_name': 'openff-1.0.0',
 'store_wavefunction': 'none'}
Spec: openff-1.1.1
{'basis': 'smirnoff',
 'implicit_solvent': None,
 'keywords': None,
 'maxiter': 200,
 'method': 'openff-1.1.1',
 'program': 'openmm',
 'scf_properties': ['dipole',
                    'quadrupole',
                    'wiberg_lowdin_indices',
                    'mayer_indices'],
 'spec_description': 'A default spec for openff-1.1.1',
 'spec_name': 'openff-1.1.1',
 'store_wavefunction': 'none'}
Spec: openff-1.2.1
{'basis': 'smirnoff',
 'implicit_solvent': None,
 'keywords': None,
 'maxiter': 200,
 'method': 'openff-1.2.1',
 'program': 'openmm',
 'scf_properties': ['dipole',
                    'quadrupole',
                    'wiberg_lowdin_indices',
                    'mayer_indices'],
 'spec_description': 'A default spec for openff-1.2.1',
 'spec_name': 'openff-1.2.1',
 'store_wavefunction': 'none'}
Spec: openff-1.3.0
{'basis': 'smirnoff',
 'implicit_solvent': None,
 'keywords': None,
 'maxiter': 200,
 'method': 'openff-1.3.0',
 'program': 'openmm',
 'scf_properties': ['dipole',
                    'quadrupole',
                    'wiberg_lowdin_indices',
                    'mayer_indices'],
 'spec_description': 'A default spec for openff-1.3.0',
 'spec_name': 'openff-1.3.0',
 'store_wavefunction': 'none'}
Spec: openff-2.0.0
{'basis': 'smirnoff',
 'implicit_solvent': None,
 'keywords': None,
 'maxiter': 200,
 'method': 'openff-2.0.0',
 'program': 'openmm',
 'scf_properties': ['dipole',
                    'quadrupole',
                    'wiberg_lowdin_indices',
                    'mayer_indices'],
 'spec_description': 'A default spec for openff-2.0.0',
 'spec_name': 'openff-2.0.0',
 'store_wavefunction': 'none'}
Spec: gaff-2.11
{'basis': 'antechamber',
 'implicit_solvent': None,
 'keywords': None,
 'maxiter': 200,
 'method': 'gaff-2.11',
 'program': 'openmm',
 'scf_properties': ['dipole',
                    'quadrupole',
                    'wiberg_lowdin_indices',
                    'mayer_indices'],
 'spec_description': 'A default spec for gaff-2.11',
 'spec_name': 'gaff-2.11',
 'store_wavefunction': 'none'}
```

