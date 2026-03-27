# tmQM xtb Dataset T=100K low-mw high-coordinate geom-mult=3 v0.0

### Description

This dataset was generated starting from an adaptation of the tmQM dataset (https://zenodo.org/records/17983516). 
This dataset contains 6,885 unique systems with 68,794 total configurations / spin states below 600 Da.  The molecules are 
limited to containing transition metals Pd, Zn, Fe, or Cu, and also only contain elements Br, C, H, P, S, O, N, F, Cl, 
or Br with charges: {-1,0,+1}. The metal is restricted to greater than three coordination sites for Pd, four for Fe, 
and one for Cu and Zn. Each molecule was preprocessed using gfn2-xtb, and then a short MD simulation
performed to provide 20 off-optimum configurations. This singlepoint dataset was then run with the BP86/def2-TZVP 
for with those geometries from molecular dynamics using gfn-xtb. Each configuration is reported with the following
properties: 'energy', 'gradient', 'dipole', 'quadrupole', 'wiberg_lowdin_indices', 'mayer_indices', 'lowdin_charges'
'dipole_polarizabilities', 'mulliken_charges'.

### General Information

- Date: 2025-12-19
- Purpose: BP86/def2-TZVP Conformers for single metal complexes with Pd, Fe, Zn, Cu, and change of {-1,0,+1} and multiplicity of 3. MW <= 600 Da, generally high coordinate, and 20 geometry samples
- Dataset Type: singlepoint
- Name: tmQM xtb Dataset T=100K low-mw high-coordinate geom-mult=3 v0.0
- Number of unique molecules: 6,885
- Number of filtered molecules: 0
- Number of Conformers: 68,794
- Number of conformers (min mean max): 3, 30, 31
- Molecular Weight (min mean max): 95 462 600
- Multiplicities: 3
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
- `tmqm_m3.txt`: List of structure labels in the HDF5 input to add to the dataset
 
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