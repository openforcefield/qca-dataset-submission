# TM Benchmark Optimization Dataset Step 1 v0.0

### Description

This dataset includes single metal complexes with: {'Pd', 'Fe', 'Zn', 'Mg', 'Cu', 'Li'}, and the non-metals: {'C', 'H', 'P', 'S', 'O', 'N', 'F', 'Cl', 'Br'}, with a complex charge of {-1,0,+1}. Additionally, there are some organic molecules for benchmarking purposes. These complexes are optimized using ROHF SOS-MP2 / def2-mSVP and  ROHF / STO-3G in suit with the literature. There are two ROHF SOS-MP2 / def2-mSVP specifications with and without frozen core. The molecular weight min, mean, and max are 81, 445, and 1026, respectively. There are 81 unique molecules, each tmc is submitted with 3 different multiplicities to assess the spin state.

Geometries were sources from the PDB CCD and multiple sources in the literature, the DOIs include: 10.1134/S0022476620090103, 10.1021/acs.inorgchem.7b03000, 10.1016/j.molstruc.2022.132506, 10.1107/S2053229619001396, 10.1021/om0492045, 10.1107/S0108270113021148, 10.1016/j.inoche.2013.06.007, and 10.1016/j.ijbiomac.2023.125847.

### General Information

- Date: 2025-04-03
- Class: OpenFF Optimization Dataset
- Purpose: Diverse set of conformers for single metal complexes with Pd, Fe, Zn, Cu, Mg, Li and charge of {-1,0,+1}, with some organic molecules for benchmarking purposes at highest level of theory
- Dataset Type: optimization
- Name: TM Benchmark Optimization Dataset Step 1 v0.0
- Number of unique molecules:   81
- Number of filtered molecules: 0
- Number of Conformers: 81
- Number of conformers (min mean max): 1 1 1
- Number of multiplicities per molecule (min mean max): 1 2 3
- Molecular Weight (min mean max): 81 445 1026
- Set of charges: -1, 0, 1
- Dataset Submitter: Jennifer A. Clark
- Dataset Curator: Jennifer A. Clark

### QCSubmit generation pipeline

- `generate-dataset.ipynb`: A python notebook which shows how the dataset was prepared from the input files.

### QCSubmit Manifest

- `generate-dataset.ipynb`
- `CS-A_primary.json`: A file containing the subset of transition metal complexes from the PDB CCD that we will include with our dataset 
- `charged_organic.json`: A file containing the subset of charged organic molecules from the PDB CCD that we will include with our dataset
- `cifs_misc/*`: A directory of cif files taken from the literature. DOIs provided in generate-dataset.ipynb
- `mol2_misc/*`: A directory of mol2 files taken from the literature. DOIs provided in generate-dataset.ipynb
- `dataset.pdf`: A pdf file containing molecule 2D structures.
- `environment.yml`: Conda environment file to perform this workflow
- `environment_full.yml`: All installed packages with versions for successful completion of this workflow
- `scaffold.py`: Python script to export a QCFractal dataset to a json file
- `scaffold.json.bz2`: A compressed json file of the target dataset
 
### Metadata

* Elements: {'Br', 'C', 'Cl', 'Cu', 'F', 'Fe', 'H', 'Li', 'Mg', 'N', 'O', 'P', 'Pd', 'S', 'Zn'}
* QC Specification: sos-mp2/def2-msvp
  * program: psi4
  * method: sos-mp2
  * basis: def2-svp
  * implicit_solvent: None
  * maxiter: 500
  * reference: rohf
  * opt_coordinates: both
  * scf_type: df
  * mp2_type: df
  * print: 3
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
* QC Specification: sos-mp2/def2-msvp FC
  * program: psi4
  * method: sos-mp2
  * basis: def2-svp
  * implicit_solvent: None
  * maxiter: 500
  * reference: rohf
  * opt_coordinates: both
  * scf_type: df
  * mp2_type: df
  * print: 3
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
* QC Specification: hf/sto-3g
  * program: psi4
  * method: hf
  * basis: sto-3g
  * implicit_solvent: None
  * maxiter: 500
  * reference: rohf
  * opt_coordinates: both
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
