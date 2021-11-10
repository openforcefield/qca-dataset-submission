## Description

Single point energies of final geometries of TorsionNet500 dataset from https://github.com/PfizerRD/TorsionNet with openff default qc specification for torsion benchmarking purposes. 

## General Information

 - Date: 2021.11.09
 - Class: Basic dataset 
 - Purpose: Energy calculation
 - Collection: BasicDataset
 - Name: TorsionNet500 single points
 - Number of unique molecules:        500
 - Number of filtered molecules:      0
 - Number of conformers:              12000
 - Number of conformers min mean max: 24  24.00 24
 - Set of charges: [0.0]
 - Dataset Submitter: Pavan Behara/Josh Horton
 - Dataset Source: https://github.com/PfizerRD/TorsionNet/tree/main/data

## Changelog

Here any information regarding dataset changes are recorded.

## QCSubmit generation pipeline

 - `Dataset_Generation.ipynb`: This notebook shows how the dataset was prepared from the input files. 
 
## QCSubmit Manifest

- `Dataset_Generation.ipynb`: Dataset creation notebook with instructions for submission.
- `dataset.json.bz2`: The compressed constrained optimization dataset ready for submission.
- `TorsionNet500.smi`: The smiles file of the peptide molecules.
- `TorsionNet500.pdf`: A pdf file containing molecule 2D structures.
- `TorsionNet500_qm_opt_geometries.sdf.gz`: File which contains structures and is the main input for the dataset
 
## Metadata

- elements {'H', 'O', 'F', 'S', 'N', 'Cl', 'C'} 
- unique molecules 500
- scf properties:
    - dipole
    - quadrupole
    - wiberg_lowdin_indices
    - mayer_indices
    - mbis_charges
- qc spec
    - name: B3LYP-D3BJ/DZVP
    - method: B3LYP-D3BJ
    - basis: DZVP
