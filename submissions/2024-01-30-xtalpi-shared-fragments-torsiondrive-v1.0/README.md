# XtalPi Shared Fragments TorsiondriveDataset v1.0

## Description

A dataset containing representative fragments shared by XtalPi
used in fitting the XFF force field
(DOI: [10.1021/acs.jctc.3c00920](https://doi.org/10.1021/acs.jctc.3c00920)).
Conformers are the post-optimization geometries shared by XtalPi.
Each conformer will be converged according to the 'GAU_LOOSE' criteria.

MOL2 files were first converted into a PyArrow dataset
using the Python script included.
This dataset is read in the Jupyter notebook to create the submission.

### Change log
* v1.0 Initial submission

## General information

* Date: 2024-01-30
* Class: OpenFF Optimization Dataset
* Purpose: Geometry optimization
* Name: XtalPi Shared Fragments TorsiondriveDataset v1.0
* Number of unique molecules: 198
* Number of conformers (min, mean, max): 1 12 30
* Molecular weight (min, mean, max): 78.50 201.15 314.38
* Charges: -2, -1, 0
* Dataset submitter: Lily Wang
* Dataset generator: Lily Wang


## QCSubmit generation pipeline

* `generate-dataset.ipynb`: This notebook shows how the dataset was prepared from the input files.


## QCSubmit Manifest

* `dataset.json.bz2`: Compressed dataset ready for submission
* `dataset.pdf`: Visualization of dataaset molecules
* `dataset.smi`: Smiles strings for dataset molecules
* `generate-dataset.ipynb`: Notebook describing dataset generation and submission
* `input-environment.yaml`: Environment file used to create Python environment for the notebook
* `full-environment.yaml`: Fully-resolved environment used to execute the notebook.
* `convert-to-td-datset.py`: Python script used to parse MOL2 files.


## Metadata

* elements: {'C', 'H', 'Cl', 'Br', 'S', 'O', 'F', 'N', 'P'}
* unique molecules: 198
* Spec: default
    * SCF properties:
        * dipole
        * quadrupole
        * wiberg_lowdin_indices
        * mayer_indices
    * QC Spec:
        * name: default
        * method: B3LYP-D3BJ
        * basis: DZVP
