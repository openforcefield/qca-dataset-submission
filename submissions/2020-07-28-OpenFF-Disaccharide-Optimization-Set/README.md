### Description

This is a Optiimization dataset containing the first 1000 disaccharide constrained optimizations.



### General Information 

- Date: 2020-09-01
- Class: OpenFF optimization 
- Purpose: a constrained optimization dataset for disaccharides
- Collection: OptimizationDataset
- Name: OpenFF-Disaccharides-v1.0
- Number of Entries: 976
- Submitter: Joshua Horton
 

### Manifest

- `Dataset Prep.ipynb`: A notebook showing how qcsubmit was used to build the dataset from the mdgx json files.
- `dataset.json`: The optimization dataset ready for submission.
- `optimization_settings.yaml`: The QCSubmit settings used to generate the optimization dataset.
- `Disaccharide.smi`:  The smiles file of the dataset.
- `Disaccharide.pdf`: A pdf file containing molecule 2D structures.


### Metadata

- elements {'C', 'H', 'O'}
- unique molecules 976
- optimizations 976
- scf properties:
    - dipole
    - quadrupole
    - wiberg_lowdin_indices
    - mayer_indices
- qc spec
    - name: default
    - method: B3LYP-D3BJ
    - basis: DZVP
