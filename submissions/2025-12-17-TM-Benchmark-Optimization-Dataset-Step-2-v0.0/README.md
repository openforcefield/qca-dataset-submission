# TM Benchmark Optimization Dataset Step 2 v0.0

### Description

  This dataset includes single metal complexes with: {'Pd', 'Fe', 'Zn', 'Mg', 'Cu', 'Li'}, and the non-metals:
  {'C', 'H', 'P', 'S', 'O', 'N', 'F', 'Cl', 'Br'}, with a complex charge of {-1,0,+1}. Additionally, there are 
  some organic molecules for benchmarking purposes. All molecules are taken from those completed optimizations 
  in the optimization dataset "TM Benchmark Optimization Dataset Step 1 v0.0". These complexes are optimized 
  at a higher level of theory with SCS-MP2 and SCS-MP3 with and without frozen core. The molecular weight min,
  mean, and max are 81, 420, and 1026, respectively. There are 75 unique molecules, each transition metal 
  complex is submitted with multiplicities that completed on the previous step to assess the spin state.

### Change Log

**2026-13-03**: Added `generate-dataset_2_update_spec.ipynb` to update specifications and basis sets to leverage analytical gradients in psi4.

### General Information

- Date: 2025-12-17
- Class: OpenFF Optimization Dataset
- Purpose: Diverse set of conformers for single metal complexes with Pd, Fe, Zn, Cu, Mg, Li and charge of {-1,0,+1}, with some organic molecules, all undergoing step 2, high level of theory final optimization.
- Dataset Type: optimization
- Name: TM Benchmark Optimization Dataset Step 2 v0.0
- Number of unique molecules:   163
- Number of filtered molecules: 0
- Number of Conformers: 163
- Number of conformers (min mean max): 1 1 1
- Number of multiplicities per molecule (min mean max): 1 2 3
- Molecular Weight (min mean max): 81 421 1026
- Set of charges: -1, 0, 1
- Dataset Submitter: Jennifer A. Clark
- Dataset Curator: Jennifer A. Clark

### QCSubmit generation pipeline

- `generate-dataset.ipynb`: A python notebook which shows how the dataset was prepared from the input files.
- `generate-dataset_2_update_spec.ipynb`: A python notebook that removes the records in the previous notebook, and creates new ones with revised specs

### QCSubmit Manifest

- `generate-dataset.ipynb`
- `generate-dataset_2_update_spec.ipynb`
- `environment.yml`: Conda environment file to perform this workflow
- `environment_full.yml`: All installed packages with versions for successful completion of this workflow
- `scaffold.json.bz2`: A compressed json file of the target dataset
- `scaffold_mp2_mp3.json.bz2`: A compressed json file of the target dataset after running `generate-dataset_2_update_spec.ipynb`

### Metadata
* Elements: {'Br', 'C', 'Cl', 'Cu', 'F', 'Fe', 'H', 'Li', 'Mg', 'N', 'O', 'P', 'Pd', 'S', 'Zn'}
* QC Specification: mp2/def2-qzvppd
  * program: psi4
  * method: mp2
  * basis: def2-qzvppd
  * implicit_solvent: None
  * maxiter: 500
  * reference: uhf
  * print: 3
  * properties_origin: ['COM']
  * freeze_core: False
  * SCF Properties:
    * dipole
    * quadrupole
    * wiberg_lowdin_indices
    * mayer_indices
    * lowdin_charges
    * mulliken_charges
  * Function Kwargs
    * Properties
      * dipole_polarizabilities
* QC Specification: mp2/def2-qzvppd-fc
  * program: psi4
  * method: mp2
  * basis: def2-qzvppd
  * implicit_solvent: None
  * maxiter: 500
  * reference: uhf
  * print: 3
  * properties_origin: ['COM']
  * freeze_core: True
  * SCF Properties:
    * dipole
    * quadrupole
    * wiberg_lowdin_indices
    * mayer_indices
    * lowdin_charges
    * mulliken_charges
  * Function Kwargs
    * Properties
      * dipole_polarizabilities
* QC Specification: mp3/def2-tzvppd
  * program: psi4
  * method: mp3
  * basis: def2-tzvppd
  * implicit_solvent: None
  * maxiter: 500
  * reference: uhf
  * print: 3
  * properties_origin: ['COM']
  * freeze_core: False
  * SCF Properties:
    * dipole
    * quadrupole
    * wiberg_lowdin_indices
    * mayer_indices
    * lowdin_charges
    * mulliken_charges
  * Function Kwargs
    * Properties
      * dipole_polarizabilities
* QC Specification: mp3/def2-tzvppd-fc
  * program: psi4
  * method: mp3
  * basis: def2-tzvppd
  * implicit_solvent: None
  * maxiter: 500
  * reference: uhf
  * print: 3
  * properties_origin: ['COM']
  * freeze_core: True
  * SCF Properties:
    * dipole
    * quadrupole
    * wiberg_lowdin_indices
    * mayer_indices
    * lowdin_charges
    * mulliken_charges
  * Function Kwargs
    * Properties
      * dipole_polarizabilities
* GeometricProcedure
  * tmax: 0.3,
  * check: 0,
  * qccnv: False,
  * reset: True,
  * trust: 0.1,
  * molcnv: False,
  * enforce: 0.0,
  * epsilon: 1e-05,
  * maxiter: 500,
  * coordsys: 'dlc',
  * constraints: {},
  * convergence_set: 'GAU'
