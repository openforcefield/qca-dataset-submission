# OpenFF Cresset Additional Coverage TorsionDrives v4.0

## Description

Molecules contributed by Cresset to address lack of data coverage for particular torsion drives.
These molecules are contributed to address five failure cases:

1) Sage 2.1 t48a is fit only to one molecule,
where the complementary torsion contributes most of the profile

2) Sage 2.1 t17 may benefit from splitting and a different shape for non-symmetric rings;
the current n=3 shape sums to a constant

3) Sage 2.1 t19 is mostly trained to terminal methyls.
It has a functional form where the n=1 term dominates unexpectedly,
instead of the more expected equal n=3 contributions.
More data with non-terminal methyls is added

4) Sage 2.1 t18 covers amide-adjacent torsions but is not trained to many.

5) Sage 2.1 t105 covers an O linker with an sp2 or sp3 terminus.
While the sp3 profiles match the QM well, the sp2 profiles look too stiff.
This adds more data.

## General Information

* Date: 2025-02-12
* Class: OpenFF TorsionDrive Dataset
* Purpose: Improve torsiondrive coverage
* Name: OpenFF Cresset Additional Coverage TorsionDrives v4.0
* Number of unique molecules: 70
* Number of driven torsions: 82
* Number of filtered molecules: 0
* Number of conformers: 171
* Number of conformers per molecule (min, mean, max): 1, 2.09, 5
* Mean molecular weight: 145.98
* Max molecular weight: 280.75
* Charges: [-1.0, 0.0, 1.0]
* Dataset generator: Matthew Habgood
* Dataset submitter: Lily Wang

## QCSubmit Generation Pipeline

* `generate-dataset.ipynb`: This notebooks shows how the dataset was prepared from the
  input files in `inputs`.

## QCSubmit Manifest

### Input files
* `inputs/*.smi`: Input SMILES strings for dataset molecules
* `inputs/*.pdf`: PDFs of molecules separated by input file type
* `generate-dataset.ipynb`: Notebook describing dataset generation and submission
* `input-env.yaml`: Environment file used to create Python environment for the script
* `full-env.yaml`: Fully-resolved environment used to execute the script

### Output Files
* `dataset.json.bz2`: Compressed dataset ready for submission
* `dataset.pdf`: Visualization of dataset molecules
* `output.smi`: SMILES strings for dataset molecules

## Metadata
* Elements: {S, Br, C, N, Cl, O, F, H}
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
