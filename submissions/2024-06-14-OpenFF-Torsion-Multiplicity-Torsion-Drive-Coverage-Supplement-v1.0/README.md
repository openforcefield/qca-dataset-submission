# OpenFF Torsion Multiplicity Torsion Drive Coverage Supplement v1.0

## Description
A torsion drive data set created to improve the coverage of both existing Sage
2.2.0 proper torsion parameters and new parameters added through the torsion
multiplicity project. The molecules in this data set were partly selected from
the ChEMBL 33 database and partly generated manually to match parameters not
covered by ChEMBL.

## General Information

* Date: 2024-06-14
* Class: OpenFF TorsionDrive Dataset
* Purpose: Improver proper torsion coverage in Sage
* Name: OpenFF Torsion Multiplicity Torsion Drive Coverage Supplement v1.0
* Number of unique molecules: 58
* Number of filtered molecules: 0
* Number of conformers: 59
* Number of conformers per molecule (min, mean, max): 1, 3.08, 10
* Mean molecular weight: 174.43
* Max molecular weight: 401.33
* Charges: [0.0, 1.0, 2.0]
* Dataset submitter: Brent Westbrook
* Dataset generator: Brent Westbrook

## QCSubmit Generation Pipeline

* `generate-dataset.ipynb`: This notebook shows how the dataset was prepared
  from the input files: `dataset.smi`, `ff.offxml`, and `test.toml`.
* The list of proper torsion parameter ID and SMILES pairs in `dataset.smi` were
  collected by searching the ChEMBL database for all of the molecules matching
  the parameters of interest. The code used for these steps can be found
  [here][frag]. Some of these (those for parameters `t146j`, `t144j`, `t117k`,
  `t116i`, and `t142j`) were then truncated manually to remove large functional
  groups far from the target dihedral. Finally, the last 20 molecules in
  `dataset.smi` were designed by hand to match their corresponding parameter
  because these parameters were not matched by any molecules in ChEMBL.

## QCSubmit Manifest

### Input Files
* `generate-dataset.ipynb`: Notebook describing dataset generation and submission
* `dataset.smi`: Sequence of parameter ID, mapped SMILES, dihedral index tuples processed by the notebook
* `ff.offxml`: Draft force field with Sage 2.2.0 proper torsions split to ensure single multiplicities
* `test.toml`: Experimental input file for defining variables used throughout the QCA submission process
* `input-environment.yaml`: Environment file used to create Python environment for the notebook
* `full-environment.yaml`: Fully-resolved environment used to execute the notebook

### Output Files
* `dataset.json.bz2`: Compressed dataset ready for submission
* `dataset.pdf`: Visualization of dataset molecules
* `output.smi`: Smiles strings for dataset molecules

## Metadata
* Elements: {N, Br, H, P, Cl, O, C, S}
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

<!-- References -->
[frag]: https://github.com/ntBre/valence-fitting/tree/c1e98fb20e7a4c9622ff031d8b23fb0b1846be7d/02_curate-data/frag
