# OpenFF NSP Optimization Set 1 N v1.0


### Description

Optimization dataset to probe coverage of NSP functional groups in general, and in particular nitrogen, in this dataset. Molecules were curated from PubChem datasets, keeping molecules with HAC < 40 and other filters specified in the FILTER.OE file. This includes a broader set of molecules, many of which may not be drug-like, but would be informative for differentiating force field parameter ranges. On top of the default OpenFF spec of B3LYP-D3BJ/DZVP, another QC specification with a triple-zeta basis is added, wb97M-D3BJ/def2-TZVP, since many of the molecules have charged states.

### General Information

- Date: 2026-02-05
- Class: OpenFF Optimization Dataset
- Purpose: Assess coverage of various NSP chemistries
- Dataset Type: optimization
- Name: OpenFF NSP Optimization Set 1 N v1.0
- Number of unique molecules: 609
- Number of filtered molecules: 0
- Number of conformers: 3432
- Number of conformers per molecule (min, mean, max): 1, 5.64, 10
- Mean molecular weight: 242.01
- Min molecular weight: 67.05
- Max molecular weight: 555.24
- Charges: [-1.0, 0.0, 1.0, 2.0]
- Dataset Submitter: Pavan Behara
- Dataset Curator: Pavan Behara

## Metadata
- Elements: {O, N, Cl, S, P, C, F, Br, H, I}
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
- Spec: WB97M-D3BJ/def2-TZVP
	 - basis: def2-TZVP
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

- `generate-dataset.ipynb`
- `environment.yml`: Conda environment file to perform this workflow
- `environment_full.yml`: All installed packages with versions for successful completion of this workflow
- `dataset.json.bz2`: A compressed json file of the target dataset
- `dataset.pdf` : Visualization of the molecules in the dataset 
- `set1-N-smiles.smi`: Input smiles for the dataset generation, along with their PubChem compound IDs
- `pubchem_preprocessing_scripts`: Some utility scripts for processing PubChem data
- `nitrogen_summary_updated.md`: Set of SMARTS patterns used to query PubChem. 
