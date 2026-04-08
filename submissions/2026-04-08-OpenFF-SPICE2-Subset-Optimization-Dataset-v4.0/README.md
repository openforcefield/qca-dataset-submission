# OpenFF SPICE2 Subset Optimization Dataset v4.0


### Description

This optimization dataset targets molecules with suspicious conformers identified in the SPICE2 set and is intended
to improve geometry-level coverage for challenging charged chemistries. Details of how the SMILES strings were 
filtered can be found at <https://github.com/jaclark5/Test_SPICE/tree/main/Filter_SPICE_Suspicious>. Molecules were
ingested from `spice2_smiles.csv`, converted to OpenFF molecules with undefined stereochemistry allowed when needed,
and expanded with up to 10 conformers per molecule using the QCSubmit standard conformer generation workflow.

Computations use the OpenFF default optimization level of theory, B3LYP-D3BJ/DZVP (Psi4), with the standard OpenFF
optimization QC specification. The resulting set spans elements {P, S, C, Br, Cl, N, F, H, I, O}, includes formal
charges from -8.0 to +2.0, and covers a broad molecular-weight range (min 34.08, mean 288.23, max 1091.26 Da). This
submission is designed for force field training and validation workflows that require robust optimization data for
highly charged and chemically diverse SPICE2-derived molecules.


### General Information

- Date: 2026-04-08
- Class: OpenFF Optimization Dataset
- Purpose: Improve optimization coverage for suspicious SPICE2 conformers across challenging charged chemistries
- Dataset Type: optimization
- Name: OpenFF SPICE2 Subset Optimization Dataset v4.0
- Number of unique molecules: 626
- Number of filtered molecules: 2
- Number of conformers: 2350
- Number of conformers per molecule (min, mean, max): 1, 3.75, 10
- Mean molecular weight: 288.23
- Min molecular weight: 34.08
- Max molecular weight: 1091.26
- Charges: [-8.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0]
- Elements: {P, S, C, Br, Cl, N, F, H, I, O}
- Dataset Submitter: Jennifer A Clark
- Dataset Curator: Jennifer A Clark

### QC Specifications
- Spec: default
	 - basis: DZVP
	 - driver: gradient
	 - implicit_solvent: None
	 - keywords: {}
	 - maxiter: 200
	 - method: B3LYP-D3BJ
	 - program: psi4
	- SCF properties:
		- dipole
		- quadrupole
		- wiberg_lowdin_indices
		- mayer_indices
		- lowdin_charges
		- mulliken_charges

### QCSubmit Manifest
- `generate-dataset.ipynb`: A python notebook which shows how the dataset was prepared from the input files.
- `environment.yaml`: Conda environment file to perform this workflow
- `environment_full.yaml`: All installed packages with versions for successful completion of this workflow
- `dataset.json.bz2`: A compressed json file of the target dataset
- `dataset.pdf` : Visualization of the molecules in the dataset 
- `dataset.smi`: SMILES for every molecule in the exported dataset
- `spice2_smiles.csv`: Source SMILES file used as notebook input
