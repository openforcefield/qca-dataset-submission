## Description

This is a supplement to single point energy calculation of DES370K (initial submission here, https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2021-11-08-QMDataset-DES370K-single-points). This dataset includes rerun of calculations failed due to MBIS charge calculation convergence errors. QCSpecification excludes MBIS charges in scf_properties on this set. Detailed description on how the initial dataset is generated can be found at https://github.com/openmm/qmdataset/tree/main/des370k.

## General Information

 - Date: 2022.02.18
 - Class: Basic dataset 
 - Purpose: Energy calculation
 - Collection: BasicDataset
 - Name: SPICE DES370K Single Points Dataset Supplement v1.0
 - Number of unique molecules:        93
 - Number of filtered molecules:      0
 - Number of conformers:              3631
 - Number of conformers min mean max: 1 39.04 2677
 - Set of charges: [-1.0, 0.0, 1.0]
 - Dataset Submitter: Josh Horton, Pavan Behara, David Dotson
 - Dataset Generator: Peter Eastman
 - Dataset Source: https://github.com/openmm/qmdataset/tree/main/des370k

## Changelog

Here any information regarding dataset changes are recorded.

## QCSubmit generation pipeline

 - `Dataset_Generation.ipynb`: This notebook shows how the dataset was prepared from the input files. 
 
## QCSubmit Manifest

- `Dataset_Generation.ipynb`: Dataset creation notebook with instructions for submission.
- `dataset.json.xz`: The compressed gradient dataset ready for submission.
- `des370k_supplement.smi`: The smiles file of the peptide molecules.
- `des370k_supplement.pdf`: A pdf file containing molecule 2D structures.
- `failed_record_ids_SPICE_DES370K_Single_Points_Dataset_v1.0.txt`: Failed records from dataset SPICE DES370K Single Points Dataset v1.0
 
## Metadata

- elements {'F', 'H', 'Cl', 'S', 'I', 'Br', 'N', 'Li', 'O', 'C', 'Na'}
- unique molecules 93
- scf properties:
    - dipole
    - quadrupole
    - wiberg_lowdin_indices
    - mayer_indices
- qc spec
    - name: wB97M-D3BJ/def2-TZVPPD
    - method: wB97M-D3BJ
    - basis: def2-TZVPPD
    
