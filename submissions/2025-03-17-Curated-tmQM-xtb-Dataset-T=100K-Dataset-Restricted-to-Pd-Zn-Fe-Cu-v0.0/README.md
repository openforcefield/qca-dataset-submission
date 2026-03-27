# Curated tmQM-xtb Dataset: T=100K Dataset Restricted to Pd, Zn, Fe, Cu v0.0

### Description

This dataset was generated starting from the tmQM dataset (release 13Aug2024, https://github.com/uiocompcat/tmQM) containing 108541 unique molecules; each molecule was evaluated using gfn2-xtb, and then a short MD simulation performed to provide additional configurations of the molecules. Further details can be found in the hosting repository: https://zenodo.org/records/15015000. This dataset contains 23134 unique systems with 208141 total configurations.  This dataset is limited to systems that contain transition metals Pd, Zn, Fe, or Cu, and also only contain elements C, H, P, S, O, N, F, Cl, or Br with charges: {-1,0,+1}. Run with the BP86/def2-TZVP for singlepoint calculations at geometries from molecular dynamics using gfn-xtb, generating the properties: 'energy', 'gradient', 'dipole', 'quadrupole', 'wiberg_lowdin_indices', 'mayer_indices', 'lowdin_charges', 'dipole_polarizabilities', 'mulliken_charges'. When available, SMILES strings were provided from the repository, xyz2mol_tm, which provides SMILES from the CSD database were sanitized and SMILES strings generated with their proximity based SMILES production method using Huckel theory. These SMILES strings can be imported into RDKit for initial visualization, but will not reflect the coordinate geometries presented from tmQm.


### General Information

- Date: 2025-03-17
- Purpose: BP86/def2-TZVP Conformers for single metal complexes with Pd, Fe, Zn, Cu, Mg, Li and change of {-1,0,+1}
- Dataset Type: singlepoint
- Name: Curated tmQM-xtb Dataset: T=100K Dataset Restricted to Pd, Zn, Fe, Cu v0.0
- Number of unique molecules: 23,134
- Number of filtered molecules: 0
- Number of Conformers: 225,068
- Number of conformers (min mean max): 1, 9, 10
- Number of Conformers with Multiplicity Variation: 675,204
- Number of multiplicities per molecule (min mean max): 3, 3, 3
- Molecular Weight (min mean max): 95 590 2541
- Set of charges: -1.0, 0.0, 1.0
- Dataset Submitter: Jennifer A. Clark
- Dataset Curator: Christopher R. Iacovella

### QCSubmit generation pipeline

- `generate-dataset.ipynb`: A python notebook which shows how the dataset was prepared from the input files.

### QCSubmit Manifest

- `generate-dataset.ipynb`
- `environment.yml`: Conda environment file to perform this workflow
- `environment_full.yml`: All installed packages with versions for successful completion of this workflow
- `scaffold.py`: Python script to export a QCFractal dataset to a json file
- `scaffold.json.bz2`: A compressed json file of the target dataset
 
### Metadata

* Elements: {'Br', 'C', 'Cl', 'Cu', 'F', 'Fe', 'H', 'N', 'O', 'P', 'Pd', 'S', 'Zn'}
* QC Specifications: BP86/def2-TZVP
  * program: psi4
  * method: BP86
  * basis: def2-TZVP
  * driver: gradient
  * implicit_solvent: None
  * keywords: {}
  * maxiter: 500
  * SCF Properties:
    * dipole
    * quadrupole
    * wiberg_lowdin_indices
    * mayer_indices
    * lowdin_charges
    * dipole_polarizabilities
    * mulliken_charges