# OpenFF NSP Optimization Set 2 Phosphorus v4.0


### Description

Optimization dataset to probe coverage of nitrogen, sulfur, and phosphorus (i.e., NSP) functional groups in general, and phosphorus in particular. Molecules were curated from PubChem datasets, retaining those with HAC < 40 and applying additional filters; the preprocessing steps are detailed at https://github.com/pavankum/NSP_sets. This dataset includes a broader set of molecules, many of which may not be drug-like, but are informative for differentiating force field parameter ranges. In addition to the default OpenFF QC specification (B3LYP-D3BJ/DZVP), an SPICE-level QC specification (ωB97M-D3BJ/def2-TZVPPD) is included, since many molecules have charged states. Set2 prioritizes as many varied charged states as possible with valid valences for the atoms.

### General Information

- Date: 2026-02-17
- Class: OpenFF Optimization Dataset
- Purpose: Assess coverage of various NSP chemistries
- Dataset Type: optimization
- Name: OpenFF NSP Optimization Set 2 Phosphorus v4.0
- Number of unique molecules: 846
- Number of filtered molecules: 1
- Number of conformers: 6115
- Number of conformers per molecule (min, mean, max): 1, 7.23, 10
- Mean molecular weight: 304.30
- Min molecular weight: 90.10
- Max molecular weight: 598.71
- Charges: [-6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 6.0]
- Elements: {C, I, H, N, Cl, P, S, F, O, Br}
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


### QCSubmit Manifest
- `generate-dataset.ipynb`: A python notebook which shows how the dataset was prepared from the input files.
- `environment.yaml`: Conda environment file to perform this workflow
- `environment_full.yaml`: All installed packages with versions for successful completion of this workflow
- `dataset.json.bz2`: A compressed json file of the target dataset
- `dataset.pdf` : Visualization of the molecules in the dataset 
- `set2-P-smiles.smi`: Input smiles for the dataset generation, along with their PubChem compound IDs
