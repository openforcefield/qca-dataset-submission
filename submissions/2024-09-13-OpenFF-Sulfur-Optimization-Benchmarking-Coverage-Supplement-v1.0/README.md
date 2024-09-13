# OpenFF Sulfur Optimization Benchmarking Coverage Supplement v1.0

## Description

An optimization data set created to improve the benchmarking coverage of
sulfonic and phosphonic acids, sulfone, sulfonate, sulfinyl, sulfoximine,
sulfonamides, thioether, and 1,3-thiazole groups. The molecules in this data
set were manually selected from a subset of the smallest matching structures in
the ChEMBL 34 database.

## General Information

* Date: 2024-09-13
* Class: OpenFF Optimization Dataset
* Purpose: Improve coverage in Sage
* Name: OpenFF Sulfur Optimization Benchmarking Coverage Supplement v1.0
* Number of unique molecules: 
* Number of filtered molecules: 
* Number of conformers: 
* Number of conformers per molecule (min, mean, max): , , 
* Mean molecular weight: 
* Max molecular weight: 
* Charges: 
* Dataset submitter: Brent Westbrook
* Dataset generator: Brent Westbrook

## QCSubmit Generation Pipeline

* `generate-dataset.py`: This script shows how the dataset was prepared from the input file `bench.smi`.
* The list of labels and SMILES pairs in `bench.smi` were collected by searching
the ChEMBL database for all of the molecules matching the SMIRKS patterns
corresponding to the labels in `sulfur.dat`. The code used for all of these
steps can be found
[here](https://github.com/ntBre/curato/tree/64261e2261e5b3109223c7fbe8ef5d866937fd13).

## QCSubmit Manifest

* `generate-dataset.py`: Script describing dataset generation and submission
* `input-environment.yaml`: Environment file used to create the Python environment for the script
* `full-environment.yaml`: Fully-resolved environment used to execute the script
* `opt.toml`: Experimental [qcaide](https://github.com/ntBre/qcaide) input file for defining
variables used throughout the QCA submission process
* `dataset.json.bz2`: Compressed dataset ready for submission
* `dataset.pdf`: Visualization of dataset molecules
* `output.smi`: SMILES strings for dataset molecules

## Metadata

* elements: 
* Spec: 
