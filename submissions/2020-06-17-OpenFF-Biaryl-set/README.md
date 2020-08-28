### Description

This is a TorsionDrive dataset consisting of biaryl torsions provided by Christopher Rowley. Originally used to benchmark parsley, but could also be useful for fitting. 

### General Information
 - Date: 17/6/2020
 - Class: OpenFF torsiondrive 
 - Purpose: Torsiondrives 
 - Collection: TorsionDriveDataset
 - Name: Openff-Rowley-Biaryl-set
 - Number of Entries: 87
 - Submitter: Josh Horton
 
 ### QCSubmit generation pipeline
 - `QCSubmit workflow.ipynb`: This notebook shows how the torsiondrive dataset was prepared from the input files. 
 
 ### QCSubmit Manifest
- `QCSubmit workflow.ipynb`: Dataset creation notebook with instructions for submission.
- `biaryl_dataset.json`: The torsiondrive dataset ready for submission.
- `biaryl_settings.yaml`: The QCSubmit settings used to generate the torsiondrive dataset.
- `biaryls.smi`: The smiles file of the biaryl molecules.
- `biaryls_dataset.pdf`: A pdf file containing molecule 2D structures with targeted torsions for driving highlighted.
 
 ### Metadata
- elements {'N', 'S', 'C', 'H', 'O'}
- unique molecules 87
- torsiondrives 87
- scf properties:
    - dipole
    - quadrupole
    - wiberg_lowdin_indices
    - mayer_indices
- qc spec
    - name: default
    - method: B3LYP-D3BJ
    - basis: DZVP
    

 