# OpenFF Organometallics Exploratory Optimization Dataset

## Description

An optimization dataset created to test the OpenFF and QCArchive infrastructure
for calculations involving organometallic molecules. The molecules in this
dataset were extracted from the `OpenEye SMILES` entries in the [Chemical
Component Dictionary](https://www.wwpdb.org/data/ccd) mmCIF file. These were
filtered to remove molecules with radical electrons and to include only
molecules with the desired metal atoms: Pd, Fe, Zn, Mg, Cu, Li, Rh, Ir, Pt, Ni,
Cr, and Ag. These were further filtered to retain only molecules with at least
10 atoms, an absolute charge of less than 4, and those not present in any of our
existing training data. From this candidate set, the molecules were sorted based
on their number of atoms, and the smallest 100 were selected. Of these, 56 were
further removed by errors in the dataset preparation process, leaving 44
molecules.

## General Information

* Date: 2024-12-03
* Class: OpenFF Optimization Dataset
* Purpose: Provide training data for metal-containing molecules
* Name: OpenFF Organometallics Exploratory Optimization Dataset
* Number of unique molecules: 44
* Number of filtered molecules: 55
* Number of conformers: 239
* Number of conformers per molecule (min, mean, max): 1, 5.43, 10
* Mean molecular weight: 424.38
* Max molecular weight: 741.40
* Charges: [0.0, 1.0, 2.0, 3.0]
* Dataset submitter: Brent Westbrook
* Dataset generator: Brent Westbrook

## QCSubmit Generation Pipeline

* `main.py`: This script shows how the dataset was prepared from `components.cif`, retrieved
from the CCD, and `inchis.dat`, which contains the InCHI keys of our existing
training data.


## QCSubmit Manifest

* `main.py`: Script describing dataset generation and submission
* `input-environment.yaml`: Environment file used to create the Python environment for the script
* `full-environment.yaml`: Fully-resolved environment used to execute the script
* `opt.toml`: Experimental [qcaide](https://github.com/ntBre/qcaide) input file for defining
variables used throughout the QCA submission process
* `dataset.json.bz2`: Compressed dataset ready for submission
* `dataset.pdf`: Visualization of dataset molecules
* `output.smi`: SMILES strings for dataset molecules

## Metadata
* Elements: {F, P, O, C, Zn, N, Ni, Pt, S, Pd, Mg, Br, Rh, Fe, H, Cl, B, Li}
* Spec: BP86/def2-TZVP
	* basis: def2-TZVP
	* implicit_solvent: None
	* keywords: {}
	* maxiter: 200
	* method: BP86
	* program: psi4
	* SCF properties:
		* dipole
		* quadrupole
		* wiberg_lowdin_indices
		* mayer_indices
