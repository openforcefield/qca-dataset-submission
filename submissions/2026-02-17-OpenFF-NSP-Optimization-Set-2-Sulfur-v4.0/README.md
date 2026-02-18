# OpenFF NSP Optimization Set 2 Sulfur v4.0


### Description

Optimization dataset to probe coverage of nitrogen, sulfur, and phosphorus (i.e., NSP) functional groups in general, and sulfur in particular. Molecules were curated from PubChem datasets, retaining those with HAC < 40 and applying additional filters; the preprocessing steps are detailed at https://github.com/pavankum/NSP_sets. This dataset includes a broader set of molecules, many of which may not be drug-like, but are informative for differentiating force field parameter ranges. In addition to the default OpenFF QC specification (B3LYP-D3BJ/DZVP), an SPICE-level QC specification (ωB97M-D3BJ/def2-TZVPPD) is included, since many molecules have charged states. Set 2 prioritizes as many varied charged states as possible with valid valences for the atoms. Some of these may include charged states that get protonated immediately in solvent, but they were retained for experiments to check the quality of charge equilibration in MLIPs as some of these may be used as reference for bespoke torsion fits. 


### General Information
- Date: 2026-02-17
- Class: OpenFF Optimization Dataset
- Purpose: Assess coverage of various NSP chemistries
- Dataset Type: optimization
- Name: OpenFF NSP Optimization Set 2 Sulfur v4.0
- Number of unique molecules: 1213
- Number of filtered molecules: 0
- Number of conformers: 8018
- Number of conformers per molecule (min, mean, max): 1, 6.61, 10
- Mean molecular weight: 293.47
- Min molecular weight: 86.16
- Max molecular weight: 597.88
- Charges: [-6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0]
- Elements: {O, P, I, C, F, N, Br, S, H, Cl}
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
- `set2-S-smiles.smi`: Smiles before pKa normalization, which includes charges [-8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0]
- `set2-S-smiles-pka-normalized.smi` : Input smiles used for the dataset generation, along with their PubChem compound IDs. This set was generated after pKa normalization of set2-S-smiles.smi with Openeye filter tool.
