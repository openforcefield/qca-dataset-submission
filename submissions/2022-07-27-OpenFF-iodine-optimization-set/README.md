# QMDataset Dipeptides

### Description

Optimization set created from Gen1 and Gen2 molecules containing iodine. Earlier datasets with iodine containing molecules and openff-default QC specification of B3LYP-D3BJ/DZVP have a discrepancy in the auxiliary basis set of Iodine, discovered by Bill Swope and Trevor Gokey. This was corrected since psi4-1.4rc3. 

### General Information

- Date: 2022-07-27
- Class: OpenFF OptimizationDataset
- Purpose: QM dataset of force field fitting
- Collection: OptimizationDataset
- Name: OpenFF Iodine Chemistry Optimization Dataset v1.0
- Number of unique molecules        68
- Number of filtered molecules      0
- Number of conformers              250
- Number of conformers min mean max 1   3.68 22
- Mean molecular weight: 318.83
- Max molecular weight: 440.19
- Charges: [0.0]
- Dataset Submitter: Pavan Behara
- Dataset Generator: Pavan Behara

### QCSubmit generation pipeline

- `Dataset-Generation.ipynb`: A notebook which shows how the dataset was prepared from the input files. 
- `dataset.json.bz2`: The optimization dataset ready for submission.
- `dipeptides.pdf`: A pdf file containing molecule 2D structures.
- `dipeptides.smi`: SMILES for every molecule in the submission.
 
## Metadata

- elements {'C', 'F', 'O', 'H', 'Br', 'Cl', 'N', 'I', 'S'}
- unique molecules 68
- scf properties:
    - dipole
    - quadrupole
    - wiberg_lowdin_indices
    - mayer_indices
    - mbis_charges
- qc spec
    - name: default
    - method: b3lyp-d3bj
    - basis: dzvp

