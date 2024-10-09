# OpenFF Lipid Optimization Training Supplement v1.0

## Description

An optimization data set created to improve the training coverage of lipid-like
molecules in Sage. The molecules in this data set were selected from the [LIPID
MAPS](https://www.lipidmaps.org/) database via
[cura](https://github.com/ntBre/curato/tree/6039bea5c64f8cd6b374fd12b5fa3b355898d98b)
after being fragmented by the [XFF
fragmentation](https://github.com/XtalPi-XFF/2023_XFF_paper/tree/b84ef6079d24ebd2e86f78a495e3257b375255fa/fragmentation)
algorithm, clustered on the 2048-bit, radius-3 Morgan fingerprint from RDKit by
the `LeaderPicker.LazyBitVectorPick` algorithm, also from RDKit, with a distance
threshold of 0.6559. Placeholder dummy atoms in the fragments were generally
filled with H atoms, and `S(=O)(=O)*` groups were replaced with `S(=O)([O-])`.
The candidate fragments were further restricted to those with between 3 and 70
heavy atoms and containing only the elements Cl, P, Br, I, H, C, O, N, F, and S.
Candidates with InChIKeys matching existing Sage training or benchmarking data
were also filtered out.

## General Information

* Date: 2024-10-08
* Class: OpenFF Optimization Dataset
* Purpose: Improve training coverage in Sage
* Name: OpenFF Lipid Optimization Training Supplement v1.0
* Number of unique molecules: 3971
* Number of filtered molecules: 32
* Number of conformers: 29770
* Number of conformers per molecule (min, mean, max): 1, 7.50, 10
* Mean molecular weight: 291.73
* Max molecular weight: 1016.29
* Charges: -4.0, -1.0, 0.0, 1.0
* Dataset submitter: Brent Westbrook
* Dataset generator: Brent Westbrook

## QCSubmit Generation Pipeline

* `generate-dataset.py`: This script shows how the dataset was prepared from the input file `input.smi`.
* `main.py`: This script shows how the dataset was prepared from the initial cura database.
* `inchis.dat`: The list of InChIKeys present in existing Sage training and
  benchmarking data, used by `main.py`.

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
* Elements: {I, Br, O, H, P, C, N, Cl, F, S}
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
