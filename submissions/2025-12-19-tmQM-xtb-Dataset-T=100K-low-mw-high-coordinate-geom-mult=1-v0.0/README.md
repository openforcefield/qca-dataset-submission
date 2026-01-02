# tmQM xtb Dataset T=100K low-mw high-coordinate geom-mult=1 v0.0

### Description

This dataset was generated starting from an adaptation of the tmQM dataset where each molecule was 
preprocessed using gfn2-xtb, and then a short MD simulation performed to provide 10 off-optimum 
configurations. More details of how these conformers were generated can be found on Zenodo 
(https://doi.org/10.5281/zenodo.17042449). This dataset contains 10,234 unique systems with 102,338 
total configurations / spin states, with molecular weights between 95 Da and 600 Da, with a mean of 
462 Da.  The molecules are limited to containing transition metals Pd, Zn, Fe, or Cu, and also only
contain elements Br, C, H, P, S, O, N, F, or Cl with charges: {-1,0,+1}. The metal is 
restricted to greater than three coordination sites for Pd, four for Fe, and one for Cu and Zn. 
This singlepoint dataset was run with the BP86/def2-TZVP for those geometries from molecular dynamics 
using gfn-xtb. Each configuration is reported with the following properties: 'energy', 'gradient', 'dipole', 
'quadrupole', 'wiberg_lowdin_indices', 'mayer_indices', 'lowdin_charges', 'dipole_polarizabilities', 'mulliken_charges'.

Each entry name corresponds to a molecule label from the original HDF5 files, combined with a conformer index in the
format `{label}_{id}`. In these HDF5 files, each molecule (identified by a unique `label`) contains multiple conformers,
representing different sets of coordinates. Conformers are sequentially numbered starting from 0 up to `N` for each 
molecule.

### General Information

- Date: 2025-12-19
- Purpose: BP86/def2-TZVP Conformers for single metal complexes with Pd, Fe, Zn, Cu, and charge of {-1,0,+1} and multiplicity of 1. MW <= 600 Da, generally high coordinate, and 10 geometry samples
- Dataset Type: singlepoint
- Name: tmQM xtb Dataset T=100K low-mw high-coordinate geom-mult=1 v0.0
- Number of unique molecules: 10,234
- Number of filtered molecules: 0
- Number of Conformers: 102,338
- Number of conformers (min mean max): 9 10 11
- Molecular Weight (min mean max): 95 462 600
- Multiplicities: 1
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
- `tmqm_m1.txt`: List of structure labels in the HDF5 input to add to the dataset
 
### Metadata

* Elements: {'Br', 'C', 'Cl', 'Cu', 'F', 'Fe', 'H', 'N', 'O', 'P', 'Pd', 'S', 'Zn'}
* QC Specifications: BP86/def2-TZVP
  * program: psi4
  * method: BP86
  * basis: def2-TZVP
  * driver: gradient
  * implicit_solvent: None
  * keywords: see `maxiter`, `SCF Properties` below
  * maxiter: 500
  * SCF Properties:
    * dipole
    * quadrupole
    * wiberg_lowdin_indices
    * mayer_indices
    * lowdin_charges
    * dipole_polarizabilities
    * mulliken_charges