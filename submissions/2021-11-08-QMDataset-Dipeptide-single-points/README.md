# QMDataset Dipeptides

### Description

Dipeptide dataset for ML applications.

### General Information

- Date: 2021.11.08
- Class: OpenFF Dataset
- Purpose: QM dataset of dipeptides for ML
- Collection: Dataset
- Name: SPICE Dipeptides single points
- Number of unique molecules        677
- Number of filtered molecules      0
- Number of conformers              33850
- Number of conformers min mean max 50  50  50
- Dataset Submitter: Josh Horton/Pavan Behara/David Dotson
- Dataset Generator: Peter Eastman
- Dataset Source: https://github.com/openmm/qmdataset/tree/main/dipeptides
- Set of charges: [-2.0, -1.0, 0.0, 1.0, 2.0]
- Mean molecular weight: 313.72
- Max molecular weight: 445.51
- Enumerate stereoisomers: False
- Enumerate tautomers: False
- Enumerate protomers: False

### QCSubmit generation pipeline

- `Dataset-Generation.ipynb`: A notebook which shows how the dataset was prepared from the input files. 

### QCSubmit Manifest

- `Dataset-Generation.ipynb`
- `dataset.json.bz2`: The basic dataset ready for submission.
- `dipeptides.pdf`: A pdf file containing molecule 2D structures.
- `dipeptides.smi`: SMILES for every molecule in the submission.
 

