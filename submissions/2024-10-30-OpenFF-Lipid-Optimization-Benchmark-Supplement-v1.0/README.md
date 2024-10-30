# OpenFF Lipid Optimization Benchmark Supplement v1.0

## Description

An optimization data set created to improve the training coverage of lipid-like
molecules in Sage. The molecules in this data set were selected from the [LIPID
MAPS](https://www.lipidmaps.org/) database via
[cura](https://github.com/ntBre/curato/tree/6039bea5c64f8cd6b374fd12b5fa3b355898d98b)
after being clustered on the 2048-bit, radius-3 Morgan fingerprint from RDKit by
the `LeaderPicker.LazyBitVectorPick` algorithm, also from RDKit, with a distance
threshold of 0.708. The candidate molecules were further restricted to those
with between 3 and 100 heavy atoms and containing only the elements Cl, P, Br,
I, H, C, O, N, F, and S. Candidates with InChIKeys matching existing Sage
training or benchmarking data were also filtered out.

## General Information

* Date: 2024-10-30
* Class: OpenFF Optimization Dataset
* Purpose: Improve testing coverage in Sage
* Name: OpenFF Lipid Optimization Benchmark Supplement v1.0
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

* `generate-dataset.py`: This script shows how the dataset was prepared from the input file `input.smi`.
* `main.py`: This script shows how the dataset was prepared from the initial cura database.


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
