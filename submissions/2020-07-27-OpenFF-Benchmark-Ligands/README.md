### Description

This is a torsiondrive dataset created from the [OpenFF FEP benchmark dataset](https://github.com/openmm/openmmforcefields/tree/master/openmmforcefields/data/perses_jacs_systems). The ligands are fragmented before having key torsions driven.
### General Information
 - Date: 27/7/2020
 - Class: OpenFF torsiondrive 
 - Purpose: Torsiondrives 
 - Collection: TorsionDriveDataset
 - Name: OpenFF-FEP-Benchmark-Ligands
 - Number of Entries: 87
 - Submitter: Josh Horton
 
 ### QCSubmit generation pipeline
 - `Ligand_fragments.ipynb`: This notebook shows how the torsiondrive dataset was prepared from the input sdf files. 
 
 ### QCSubmit Manifest
- `QCSubmit workflow.ipynb`: Dataset creation notebook with instructions for submission.
- `biaryl_dataset.json`: The torsiondrive dataset ready for submission.
- `biaryl_settings.yaml`: The QCSubmit settings used to generate the torsiondrive dataset.
- `biaryls.smi`: The smiles file of the biaryl molecules.
- `biaryls_dataset.pdf`: A pdf file containing molecule 2D structures with targeted torsions for driving highlighted.
 
 ## Metadata
 Each dataset is prepared twice, with and without non rotor ring substituents
 ### BACE
Without
- elements {'C', 'H', 'N', 'O', 'S'}
- unique molecules 23
- torsiondrives 38

 With
- elements {'C', 'Cl', 'F', 'H', 'N', 'O', 'S'}
- unique molecules 29
- torsiondrives 44

### CDK2
Without
- elements {'Br', 'C', 'Cl', 'H', 'N', 'O', 'S'}
- unqiue molecules 23
- torsiondrives 28

With
- elements {'Br', 'C', 'Cl', 'H', 'N', 'O', 'S'}
- unqiue molecules 25
- torsiondrives 32

### JNK1
Without
- elements {'Br', 'C', 'Cl', 'H', 'N', 'O', 'S'}
- unqiue molecules 28
- torsiondrives 39

With
- elements {'Br', 'C', 'Cl', 'H', 'N', 'O', 'S'}
- unqiue molecules 45
- torsiondrives 61

### MCL1
Without
- elements {'C', 'Cl', 'F', 'H', 'N', 'O', 'S'}
- unqiue molecules 58
- torsiondrives 59

With
- elements {'C', 'Cl', 'F', 'H', 'N', 'O', 'S'}
- unqiue molecules 105
- torsiondrives 107

### P38a
Without
- elements {'C', 'Cl', 'F', 'H', 'N', 'O', 'S'}
- unqiue molecules 49
- torsiondrives 72

With
- elements {'C', 'Cl', 'F', 'H', 'N', 'O', 'S'}
- unqiue molecules 53
- torsiondrives 80


### PTP1B
Without
- elements {'Br', 'C', 'Cl', 'F', 'H', 'N', 'O', 'S'}
- unqiue molecules 58
- torsiondrives 78

With
- elements {'Br', 'C', 'Cl', 'F', 'H', 'N', 'O', 'S'}
- unqiue molecules 61
- torsiondrives 81


### Thrombin
Without
- elements {'Br', 'C', 'Cl', 'F', 'H', 'I', 'N', 'O'}
- unqiue molecules 22
- torsiondrives 25

With
- elements {'Br', 'C', 'Cl', 'F', 'H', 'I', 'N', 'O'}
- unqiue molecules 35
- torsiondrives 38


### TYK2
Without
- elements {'C', 'Cl', 'F', 'H', 'N', 'O'}
- unqiue molecules 21
- torsiondrives 44

With
- elements {'C', 'Cl', 'F', 'H', 'N', 'O'}
- unqiue molecules 22
- torsiondrives 45 

### TOTAL
Without
- elements {'Br', 'C', 'Cl', 'F', 'H', 'I', 'N', 'O', 'S'}
- unqiue molecules 274
- torsiondrives 375

With
- elements {'Br', 'C', 'Cl', 'F', 'H', 'I', 'N', 'O', 'S'}
- unqiue molecules 368
- torsiondrives 481

- scf properties:
    - dipole
    - quadrupole
    - wiberg_lowdin_indices
    - mayer_indices
- qc spec
    - name: default
    - method: B3LYP-D3BJ
    - basis: DZVP
    

 