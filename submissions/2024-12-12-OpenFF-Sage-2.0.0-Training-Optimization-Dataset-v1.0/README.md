# OpenFF Sage 2.0.0 Training Optimization v1.0

### Description

A quantum chemical (QC) dataset curated to train [OpenFF 2.0.0 Sage](https://github.com/openforcefield/openff-sage) forcefield, with reparametrized Lennard-Jones (LJ) and valence parameters, the latter relevent to this dataset. This QC dataset with the OpenFF default level of theory, B3LYP-D3BJ/DZVP, is used to benchmark Sage geometries and energetics. These optimized conformer geometries where used in conjunction with the QC dataset used to train one dimensional torsional profiles. This Generation 2 dataset increases chemical diversity when compared to Generation 1, which are of value to our industry partners. Large molecules (>20 heavy atoms) were also included, including more flexible molecules and a greater degree of conformational variation which provide intramolecular interactions.

### General Information

- Date: 2024 12 12
- Class: OpenFF Optimization Dataset
- Purpose: B3LYP-D3BJ/DZVP conformers applicable to drug-like molecules for OpenFF 2.0.0 Sage
- Collection: OptimizationDataset
- Name: OpenFF Sage 2.0.0 Training Optimization v1.0
- Number of unique molecules       1025
- Number of filtered molecules     0 
- Number of conformers             3663
- Number of conformers min mean max 1.00, 3.53, 10.00
- Mean molecular weight: 261.38
- Max molecular weight: 544.64
- Set of charges: -2.0 -1.0 0.0 1.0
- Dataset Submitter: Jennifer A. Clark
- Dataset Curator: Simon Boothroyd
- Dataset Generator: Hyesu Jang

### QCSubmit generation pipeline

- `generate-combined-dataset.ipynb`: A notebook which shows how the dataset was prepared from the input files.

### QCSubmit Manifest

- `generate-combined-dataset.ipynb`
- `dataset.json.bz2`: The basic dataset ready for submission.
- `dataset.pdf`: A pdf file containing molecule 2D structures.
- `dataset.smi`: SMILES for every molecule in the submission.
 
### Metadata

* Elements: {F, I, N, C, P, Cl, S, Br, O, H}
* QC Specifications: default
  * basis: DZVP
  * implicit_solvent: None
  * keywords: {}
  * maxiter: 200
  * method: B3LYP-D3BJ
  * program: psi4
  * SCF Properties:
    * dipole
    * quadrupole
    * wiberg_lowdin_indices
    * mayer_indices