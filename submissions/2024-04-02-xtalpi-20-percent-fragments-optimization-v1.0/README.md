# XtalPi 20-percent Fragments OptimizationDataset v1.0

## Description

A dataset containing additional representative fragments shared by XtalPi
used in fitting the XFF force field
(DOI: [10.1021/acs.jctc.3c00920](https://doi.org/10.1021/acs.jctc.3c00920)). This is a continuation of the
[XtalPi Shared Fragments OptimizationDataset v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2024-01-30-xtalpi-shared-fragments-optimization-v1.0) dataset.
Where this previous dataset contained ~3% the original XFF fitting dataset,
this dataset contains an additional 20%. 
Conformers are the post-optimization geometries shared by XtalPi.

MOL2 files were first converted into a PyArrow dataset
using the Python script included.
This dataset is read in the Jupyter notebook to create the submission.


### Change log
* v1.0 Initial submission

## General information

* Date: 2024-04-02
* Class: OpenFF Optimization Dataset
* Purpose: Geometry optimization
* Name: XtalPi 20-percent Fragments OptimizationDataset v1.0
* Number of unique molecules: 10069
* Number of conformers: 128180
* Number of conformers (min, mean, max): 1 13 30
* Molecular weight (min, mean, max): 42.04, 189.06, 387.45
* Charges: -2, -1, 0, 1, 2
* Dataset submitter: Lily Wang
* Dataset generator: Lily Wang


## QCSubmit generation pipeline

* `generate-dataset.ipynb`: This notebook shows how the dataset was prepared from the input files.

## QCSubmit Manifest

* `dataset.json.bz2`: Compressed dataset ready for submission
* `dataset.pdf.tar.gz`: Visualization of dataset molecules, compressed
* `dataset.smi`: Smiles strings for dataset molecules
* `generate-dataset.ipynb`: Notebook describing dataset generation and submission
* `input-environment.yaml`: Environment file used to create Python environment for the notebook
* `full-environment.yaml`: Fully-resolved environment used to execute the notebook.
* `convert-to-opt-datset.py`: Python script used to parse MOL2 files.


## Metadata

* elements: {'Cl', 'P', 'Br', 'I', 'H', 'C', 'B', 'Si', 'O', 'N', 'F', 'S'}
* unique molecules: 10069
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