# tmQM xtb Dataset T=100K Pd Zn Fe Cu low mw v0.0

### Description

This dataset was generated starting from the tmQM dataset (release 13Aug2024, https://github.com/uiocompcat/tmQM). 
This dataset contains 6,829 unique systems with 471,097 total configurations/spin states below 600 Da.  The molecules are 
limited to containing transition metals d, Zn, Fe, or Cu, and also only contain elements C, H, P, S, O, N, F, Cl, 
or Br with charges: {-1,0,+1}. The metal is restricted to greater than four coordination sites, except for Cu and Zn 
must be greater than or equal to two. Each molecule was preprocessed using gfn2-xtb, and then a short MD simulation
performed to provide ~30 additional configurations of per molecules at three different spin states, 1, 3, and 5. This
singlepoint dataset was then run with the BP86/def2-TZVP for with those geometries from molecular dynamics using
 gfn-xtb. Each configuration is reported with the following properties: 'energy', 'gradient', 'dipole', 'quadrupole',
'wiberg_lowdin_indices', 'mayer_indices', 'lowdin_charges' 'dipole_polarizabilities', 'mulliken_charges'. SMILES
strings where generated from tmos (https://github.com/openforcefield/tmos) when possible. These SMILES strings can be
imported into RDKit for initial visualization, but will not reflect the coordinate geometries presented from tmQm.

### General Information

- Date: 2025-08-14
- Purpose: BP86/def2-TZVP Conformers for single metal complexes with Pd, Fe, Zn, Cu, and change of {-1,0,+1}, MW <= 600 Da, generally high coordinate, and 30 geometry samples
- Dataset Type: singlepoint
- Name: tmQM xtb Dataset T=100K Pd Zn Fe Cu low mw v0.0
- Number of unique molecules: 6,829
- Number of filtered molecules: 0
- Number of Conformers: 471,097
- Number of conformers (min mean max): 30, 68, 88
- Molecular Weight (min mean max): 95 455 600
- Set of charges: -1.0, 0.0, 1.0
- Dataset Submitter: Jennifer A. Clark
- Dataset Curator: Christopher R. Iacovella

### QCSubmit generation pipeline

- `generate-dataset.ipynb`: A python notebook which shows how the dataset was prepared from the input files.

### QCSubmit Manifest

- `generate-dataset.ipynb`
- `environment.yml`: Conda environment file to perform this workflow
- `environment_full.yml`: All installed packages with versions for successful completion of this workflow
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