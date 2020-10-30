### Description

This is a torsiondrive dataset created from the [OpenFF FEP benchmark dataset](https://github.com/openmm/openmmforcefields/tree/master/openmmforcefields/data/perses_jacs_systems). The ligands are fragmented before having key torsions driven.
### General Information
 - Date: 27/7/2020
 - Class: OpenFF torsiondrive 
 - Purpose: Torsiondrives 
 - Collection: TorsionDriveDataset
 - Name: OpenFF-FEP-Benchmark-Ligands
 - Number of Entries: 368
 - Submitter: Josh Horton
 
 ### QCSubmit generation pipeline
 - `Ligand_fragments.ipynb`: This notebook shows how the torsiondrive dataset was prepared from the input sdf files. 
 
 ### QCSubmit Manifest
- `Ligand_fragments.ipynb`: Dataset creation notebook with instructions for submission.
- `dataset.json`: The torsiondrive dataset ready for submission.
- `fragment_settings.yaml`: The QCSubmit settings used to generate the torsiondrive dataset.
- `fragments.smi`: The smiles file of the fragmented molecules.
- `fragments.pdf`: A pdf file containing molecule 2D structures with targeted torsions for driving highlighted.
 
 ## Metadata
- elements {'Br', 'C', 'Cl', 'F', 'H', 'I', 'N', 'O', 'S'}
- unqiue molecules 368
- torsiondrives 481

- scf properties:
    - dipole
    - quadrupole
    - wiberg_lowdin_indices
    - mayer_indices
    
- qc spec
    - default:
        - name: default
        - method: B3LYP-D3BJ
        - basis: DZVP
        - program: psi4
    - openff-1.0.0:
        - name: openff-1.0.0
        - method: openff-1.0.0
        - basis: smirnoff
        - program: openmm
    - gaff-2.11:
        - name: gaff-2.11
        - method: gaff-2.11
        - basis: antechamber
        - program: openmm
    - ani2x:
        - name: ani2x
        - method: ani2x
        - basis: None
        - program: torchani
    - ani2x-v2:
        - name: ani2x
        - method: ani2x
        - basis: None
        - program: torchani
    - gfn2xtb:
        - name: gfn2xtb
        - method: gfn2xtb
        - basis: None
        - program: xtb
    - gfn1xtb:
        - name: gfn1xtb 
        - method: gfn1xtb
        - basis: None
        - program: xtb
    - gfnff:
        - name: gfnff
        - method: gfnff
        - basis: None
        - program: xtb
     
    
    

 
