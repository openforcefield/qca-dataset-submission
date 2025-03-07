# OpenFF Protein PDB 4-mers V1.0

### Description

This dataset is composed of 1,000 4-mer peptide structures which were extracted from PDB entries in the Top8000 database. The purpose of this dataset is to fill in gaps within the existing protein training data with secondary structures which occur in real proteins. For all scripts used to generated the input files please visit https://github.com/ajfriedman22/4mer_Generation. The dataset was computed using the B3LYP-D3BJ method and the DZVP basis set.

### General Information
- Date: 2025-03-05
- Class: OpenFF Optimization Dataset
- Purpose: Optimizations of peptide 4-mers extracted from the PDB
- Collection: OptimizationDataset
- Name: OpenFF Protein PDB 4-mers v1.0
- Number of unique molecules        200
- Number of filtered molecules     0 
- Number of conformers             1000 
- Number of conformers min mean max 5 5 5
- Dataset Submitter: Anika Friedman
- Dataset Generator: Anika Friedman
- Set of charges: [-2.0, -1.0, 0.0, 1.0]
- Mean molecular weight: 451.35
- Max molecular weight: 570.64

### QCSubmit generation pipeline

* `generate-dataset.ipynb`: This notebooks shows how the dataset was prepared from the
  input files in `inputs`.

### QCSubmit Manifest

#### Input files
* `inputs.tar.bz2`: Input SDF files for dataset molecules
* `generate-dataset.ipynb`: Notebook describing dataset generation and submission
* `input-env.yaml`: Environment file used to create Python environment for the script

#### Output Files
* `dataset.json.bz2`: Compressed dataset ready for submission
* `dataset.pdf`: Visualization of dataset molecules
* `dataset.smi`: SMILES strings for dataset molecules
 
#### Metadata
* Elements: {'O', 'C', 'N', 'H'}
* Spec: default
  * basis: DZVP
  * implicit_solvent: None
  * keywords: {}
  * maxiter: 200
  * method: B3LYP-D3BJ
  * program: psi4
  * SCF properties:
    * dipole
    * quadrupole
    * wiberg_lowdin_indices
    * mayer_indices
