# QMDataset Dipeptides

### Description

Optimization set created from Gen1, Gen2, and other opt datasets that contain molecules with Iodine. Earlier datasets with Iodine containing molecules and openff-default QC specification of B3LYP-D3BJ/DZVP have a discrepancy in the auxiliary basis set of Iodine, discovered by Bill Swope and Trevor Gokey. This was corrected since psi4-1.4rc3. 

Molecules with Iodine from the following sets are added and the coordinates are translated by 2 bohrs so that the records won't get populated by old data on QCA: 
["OpenFF Discrepancy Benchmark 1",
"OpenFF Gen 2 Opt Set 2 Coverage",
"OpenFF Gen 2 Opt Set 4 eMolecules Discrepancy",
"SMIRNOFF Coverage Set 1",
"OpenFF Ehrman Informative Optimization v0.2",
"FDA optimization dataset 1",
"Kinase Inhibitors: WBO Distributions"]


### General Information

- Date: 2022-07-27
- Class: OpenFF OptimizationDataset
- Purpose: QM dataset of force field fitting
- Collection: OptimizationDataset
- Name: OpenFF Iodine Chemistry Optimization Dataset v1.0
- Number of unique molecules        99
- Number of filtered molecules      0
- Number of conformers              327
- Number of conformers min mean max 1   3.30 22
- Mean molecular weight: 318.53
- Max molecular weight: 533.92
- Charges: [-1.0, 0.0, 1.0]
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

