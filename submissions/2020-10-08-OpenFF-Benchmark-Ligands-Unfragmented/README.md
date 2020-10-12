## Description

This is a torsiondrive dataset created from the [OpenFF FEP benchmark dataset](https://github.com/openmm/openmmforcefields/tree/master/openmmforcefields/data/perses_jacs_systems).
The ligands are *not* fragmented before having key torsions driven.

## General Information

- Date: 2020.10.08
- Class: OpenFF TorsionDrive 
- Purpose: TorsionDrive scans for the JACS benchmark ligand set, unfragmented primarily for MM, ML potentials
- Collection: TorsionDriveDataset
- Name: OpenFF Benchmark Ligands - Unfragmented v1.0
- Number of unique molecules: 174
- Number of torsiondrives: 1255
- Submitter: David Dotson

 
## QCSubmit generation pipeline

- `Dataset Preparation.ipynb`: This notebook shows how the torsiondrive dataset was prepared from the input SDF files. 
- `dataset.json.bz2`: The dataset submission itself, for consumption by QCSubmit for submission to the public QCArchive.
- `molecules.pdf`: Visualization of all torsions exercised in this dataset.
- `torsiondrivefactory_settings.yaml`: Complete settings used for the TorsionDrive protocol.


## QCSubmit Manifest

- `Ligand_fragments.ipynb`: Dataset creation notebook with instructions for submission.
- `dataset.json`: The torsiondrive dataset ready for submission.
- `fragment_settings.yaml`: The QCSubmit settings used to generate the torsiondrive dataset.
- `fragments.smi`: The smiles file of the fragmented molecules.
- `fragments.pdf`: A pdf file containing molecule 2D structures with targeted torsions for driving highlighted.


## Metadata

- elements {'C', 'Cl', 'F', 'H', 'N', 'O', 'S'}
- unqiue molecules 174
- torsiondrives 1255

- scf properties:
    - dipole
    - quadrupole
    - wiberg_lowdin_indices
    - mayer_indices
    
- qc spec
    - openff-1.0.0:
        - name: openff-1.0.0
        - method: openff-1.0.0
        - basis: smirnoff
        - program: openmm
    - openff-1.1.0:
        - name: openff-1.1.0
        - method: openff-1.1.0
        - basis: smirnoff
        - program: openmm
    - openff-1.2.0:
        - name: openff-1.2.0
        - method: openff-1.2.0
        - basis: smirnoff
        - program: openmm
    - openff-1.2.1:
        - name: openff-1.2.1
        - method: openff-1.2.1
        - basis: smirnoff
        - program: openmm
    - ani2x:
        - name: ani2x
        - method: openff-1.2.1
        - basis: null 
        - program: torchani
