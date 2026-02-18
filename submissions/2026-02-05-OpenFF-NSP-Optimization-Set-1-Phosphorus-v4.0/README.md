# OpenFF NSP Optimization Set 1 Phosphorus v4.0


### Description

Optimization dataset to probe coverage of nitrogen, sulfur, and phosphorus (i.e., NSP) functional groups in general, and phosphorus in particular. Molecules were curated from PubChem datasets, retaining those with HAC < 40 and applying additional filters; the preprocessing steps are detailed at https://github.com/pavankum/NSP_sets. This dataset includes a broader set of molecules, many of which may not be drug-like, but are informative for differentiating force field parameter ranges. In addition to the default OpenFF QC specification (B3LYP-D3BJ/DZVP), an SPICE-level QC specification (ωB97M-D3BJ/def2-TZVPPD) is included, since many molecules have charged states. Set1 prioritizes neutral molecules.

### General Information

- Date: 2026-02-05
- Class: OpenFF Optimization Dataset
- Purpose: Assess coverage of various NSP chemistries
- Dataset Type: optimization
- Name: OpenFF NSP Optimization Set 1 Phosphorus v4.0
- Number of unique molecules: 430
- Number of filtered molecules: 0
- Number of conformers: 2910
- Number of conformers per molecule (min, mean, max): 1, 6.77, 10
- Mean molecular weight: 278.25
- Min molecular weight: 90.10
- Max molecular weight: 586.68
- Charges: [0.0]
- Elements: {H, Br, N, S, C, F, I, Cl, P, O}
- Dataset Submitter: Pavan Behara
- Dataset Curator: Pavan Behara

### QC Specifications

- Spec: default
	 - basis: DZVP
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
- Spec: WB97M-D3BJ/def2-TZVPPD
	 - basis: def2-TZVPPD
	 - implicit_solvent: None
	 - keywords: {}
	 - maxiter: 200
	 - method: WB97M-D3BJ
	 - program: psi4
	- SCF properties:
		- dipole
		- quadrupole
		- wiberg_lowdin_indices
		- mayer_indices
		- lowdin_charges
		- mulliken_charges

### QCSubmit generation pipeline
- `generate-dataset.ipynb`: A python notebook which shows how the dataset was prepared from the input files.


### QCSubmit Manifest
- `generate-dataset.ipynb`: A python notebook which shows how the dataset was prepared from the input files.
- `environment.yaml`: Conda environment file to perform this workflow
- `environment_full.yaml`: All installed packages with versions for successful completion of this workflow
- `dataset.json.bz2`: A compressed json file of the target dataset
- `dataset.pdf` : Visualization of the molecules in the dataset 
- `set1-P-smiles.smi`: Input smiles for the dataset generation, along with their PubChem compound IDs
