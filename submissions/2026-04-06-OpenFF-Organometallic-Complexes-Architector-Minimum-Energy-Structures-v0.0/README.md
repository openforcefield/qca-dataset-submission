# OpenFF Organometallic Complexes Architector Minimum Energy Structures v0.0

### Description

This dataset was generated using [architector](https://github.com/lanl/Architector/tree/Secondary_Solvation_Shell), the
details of which can be found in the Zenodo record (https://zenodo.org/records/19372923). This dataset contains 40,486 
unique systems/configurations below 1005 Da.  The molecules are limited to containing transition metals Pd, Zn, Fe, Cu, 
Li, or Mg and also only contain elements Br, C, H, P, S, O, N, F, Cl, or Br with overall charges: {-1,0,+1}. The metal
coordination for each metal center is between 1 and 12. Each molecule was preprocessed using gfn2-xtb as implemented 
in the architector package. This optimization dataset was then run with the BP86/def2-TZVP. Each configuration is
reported with the following properties: 'energy', 'gradient', 'dipole', 'quadrupole', 'wiberg_lowdin_indices', 
'mayer_indices', 'lowdin_charges', 'dipole_polarizabilities', 'mulliken_charges'.

### General Information

- Date: 2026-04-06
- Purpose: BP86/def2-TZVP single metal complex optimizations with Pd, Fe, Zn, Cu, Li, and Mg and charge of {-1,0,+1} and multiplicities ranging from 1 to 6 and coordination ranging from 1 to 12, and MW <= 1005 Da.
- Dataset Type: optimization
- Name: OpenFF Organometallic Complexes Architector Minimum Energy Structures v0.0
- Number of unique molecules: 40,486
- Number of filtered molecules: 0
- Number of Conformers: 40,486
- Number of conformers (min mean max): 1, 1, 1
- Molecular Weight (min mean max): 22 392 1005
- Multiplicities: 1, 2, 3, 4, 5, 6
- Coordination Numbers: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12]
- Oxidation States: [0, 1, 2, 3, 4, 5, 6, 7]
- Set of charges: -1.0, 0.0, 1.0
- Dataset Submitter: Jennifer A. Clark
- Dataset Curator: Jennifer A. Clark

### QCSubmit generation pipeline

- `generate_dataset.ipynb`: A python notebook which shows how the dataset was prepared from the input files.

### QCSubmit Manifest

- `generate_dataset.ipynb`
- `environment.yml`: Conda environment file to perform this workflow
- `environment_full.yml`: All installed packages with versions for successful completion of this workflow
- `scaffold.json.bz2`: A compressed json file of the target dataset
 
### Metadata

* Elements: Br, C, Cu, Fe, H, Li, Mg, N, O, P, Pd, S, Zn
* Spec: BP86/def2-TZVP
    * program: geometric
    * keywords:
       * tmax: 0.3
       * check: 0
       * qccnv: False
       * reset: True
       * trust: 0.1
       * molcnv: False
       * enforce: 0.0
       * epsilon: 1e-05
       * maxiter: 300
       * coordsys: dlc
       * convergence_set: GAU
    * qc_specification:
       * program: psi4
       * driver: SinglepointDriver.deferred
       * method: bp86
       * basis: def2-tzvp
       * keywords: {'maxiter': 500, 'scf_properties': ['dipole', 'quadrupole', 'wiberg_lowdin_indices', 'mayer_indices', 'lowdin_charges', 'mulliken_charges'], 'function_kwargs': {'properties': ['dipole_polarizabilities']}}
* SCF properties:
       * dipole
       * quadrupole
       * wiberg_lowdin_indices
       * mayer_indices
       * lowdin_charges
       * mulliken_charges