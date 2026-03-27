# tmQM xtb Dataset T=100K low-mw high-coordinate xtb mult=3 v0.0

### Description

This dataset was generated starting from an adaptation of the tmQM dataset (DOI: 10.5281/zenodo.14920177; 
https://zenodo.org/records/17042449). This dataset contains 6,905 unique systems with 206,240 total 
configurations / spin state combinations  below 600 Da.  The molecules are limited to containing transition 
metals Pd, Zn, Fe, or Cu, and also only contain elements Br, C, H, P, S, O, N, F, Cl, or Br with charges: 
{-1,0,+1}. The metal is restricted to greater than three coordination sites for Pd, four for Fe, 
and one for Cu and Zn. Each molecule was preprocessed using gfn2-xtb, and then a short MD simulation
performed to provide a maximum of 30 off-optimum configurations in addition to the minimized geometry per molecules at 
a multiplicity of 3. Using the geometries generated with gfn-xtb, this singlepoint dataset was then run with the DFT method
BP86/def2-TZVP and a multiplicity of either 1, 3, or 5. 
Each configuration is reported with the following properties: 'energy', 'gradient', 'dipole', 'quadrupole',
'wiberg_lowdin_indices', 'mayer_indices', 'lowdin_charges' 'dipole_polarizabilities', 'mulliken_charges'. SMILES
strings where generated from tmos (https://github.com/openforcefield/tmos) when possible. These SMILES strings can be
imported into RDKit for initial visualization, but will not reflect the coordinate geometries presented from tmQm.

### General Information

- Date: 2025-09-02
- Purpose: BP86/def2-TZVP Conformers for single metal complexes with Pd, Fe, Zn, Cu, and charge of {-1,0,+1} run in xtb with a multiplicity of 3 and in DFT run with a multiplicity of 1, 3, or 5. MW <= 600 Da, generally high coordinate, and a max of 30 geometry samples
- Dataset Type: singlepoint
- Name: tmQM xtb Dataset T=100K low-mw high-coordinate xtb mult=3 v0.0
- Number of unique molecules: 6,905
- Number of filtered molecules: 0
- Number of Conformers: 206,240
- Number of conformers (min mean max): 3, 30, 31
- Molecular Weight (min mean max): 190 463 600
- Set of charges: -1.0, 0.0, 1.0
- Dataset Submitter: Jennifer A. Clark
- Dataset Curator: Christopher R. Iacovella

### Changelog

- 2025-09-02: Dataset is initially created
- 2025-10-09: Update specs to use the "COM" as the reference for SCF multipole moments.
- 2026-01-13: Update dataset metadata to properly reflect the dataset contents

### QCSubmit generation pipeline

- `generate_dataset.ipynb`: A python notebook which shows how the dataset was prepared from the input files.

### QCSubmit Manifest

- `generate_dataset.ipynb`
- `environment.yml`: Conda environment file to perform this workflow
- `environment_full.yml`: All installed packages with versions for successful completion of this workflow
- `scaffold.json.bz2`: A compressed json file of the target dataset
- `update_dataset_spec_2025_10_09.ipynb`: Update specs to use the "COM" as the reference for SCF multipole moments.
- `scaffold_2.json.bz2`: A compressed json file of the dataset produced by update_dataset_spec_2025_10_09.ipynb
- `update_dataset_2026_01_13.ipynb`: Update metadata of the dataset 
 
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