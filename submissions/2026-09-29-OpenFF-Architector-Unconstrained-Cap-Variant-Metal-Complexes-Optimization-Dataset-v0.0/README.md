# OpenFF Architector Unconstrained Cap Variant Metal Complexes Optimization Dataset v0.0

### Description

This dataset was generated using [architector](https://github.com/lanl/Architector/tree/Secondary_Solvation_Shell), the
details of the HDF5 file can be found in the Zenodo record (https://zenodo.org/records/23041053). This dataset contains
22,343 unique systems/configurations below 1140 Da using the same keys as in the HDF5 as entry labels. The molecules
are limited to containing transition metals Pd, Zn, Fe, Cu, Li, or Mg with ligands capped with either methyl groups, fluorine,
or hydrogen (representing electronically neutral, withdrawing, and donating) and also only contain elements C, H, P, S, O,
N, F, Cl, or Br with overall charges: {-1,0,+1}. The metal coordination for each metal center is between 1 and 12. Each
molecule was preprocessed using gfn2-xtb as implemented in the architector package. This optimization dataset was then 
run with the BP86/def2-TZVP. Each configuration is reported with the following properties: 'energy', 'gradient', 'dipole', 
'quadrupole', 'wiberg_lowdin_indices', 'mayer_indices', 'lowdin_charges', 'lowdin_spins', 'dipole_polarizabilities', 
'mulliken_charges'. These structures are loosely optimized without constraints.

### General Information

- Date: 2026-09-29
- Purpose: Metal complexes with ligands capped with methyl, fluroine, or hydrogen, optimized with BP86/def2-TZVP.
- Dataset Type: optimization
- Name: OpenFF Architector Unconstrained Cap Variant Metal Complexes Optimization Dataset v0.0
- Number of unique molecules: 22,343
- Number of filtered molecules: 0
- Number of Conformers: 22,343
- Number of conformers (min mean max): 1, 1, 1
- Molecular Weight (min mean max): 20 339 1136
- Multiplicities: 1, 2, 3, 4, 5, 6
- Coordination Numbers: {8: 3817, 6: 2926, 7: 2733, 3: 2721, 4: 2575, 5: 2470, 9: 1847, 2: 1572, 12: 625, 10: 618, 1: 439}
- Metals: {'Pd': 6913, 'Fe': 5161, 'Zn': 3784, 'Cu': 2615, 'Mg': 2609, 'Li': 1261}
- Oxidation States: {1: 6580, 2: 6137, 3: 4589, 4: 3646, 0: 1391}
- Set of charges: -1.0, 0.0, 1.0
- Dataset Submitter: Jennifer A. Clark
- Dataset Curator: Jennifer A. Clark

### QCSubmit generation pipeline

- `generate_dataset.ipynb`: A python notebook which shows how the dataset was prepared from the input files.

### QCSubmit Manifest

- `generate_dataset.ipynb`
- `environment.yml`: Conda environment file to perform this workflow
- `environment_full.yaml`: All installed packages with versions for successful completion of this workflow
- `scaffold.json.bz2`: A compressed json file of the original target dataset
 
### Metadata

* Elements: Br, C, Cu, F, Fe, H, Li, Mg, N, O, P, Pd, S, Zn
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
       * converge: ['energy', '1e-3', 'grms', '0.2', 'gmax', '1.0', 'drms', '15', 'dmax', '30']
       * coordsys: dlc
       * convergence_set: GAU
    * qc_specification:
       * program: psi4
       * driver: SinglepointDriver.deferred
       * method: bp86
       * basis: def2-tzvp
       * keywords: {'maxiter': 500, 'reference': 'uks', 'scf_properties': ['dipole', 'quadrupole', 'wiberg_lowdin_indices', 'mayer_indices', 'lowdin_charges', 'lowdin_spins', 'mulliken_charges'], 'function_kwargs': {'properties': ['dipole_polarizabilities']}, 'properties_origin': ['COM']}
       * protocols: {'wavefunction': <WavefunctionProtocolEnum.none: 'none'>, 'stdout': True, 'error_correction': {'default_policy': True, 'policies': None}, 'native_files': <NativeFilesProtocolEnum.none: 'none'>}
    * SCF properties:
           * dipole
           * quadrupole
           * wiberg_lowdin_indices
           * mayer_indices
           * lowdin_charges
           * lowdin_spins
           * mulliken_charges
* Spec: BP86/def2-TZVP cuEST
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
       * converge: ['energy', '1e-3', 'grms', '0.2', 'gmax', '1.0', 'drms', '15', 'dmax', '30']
       * coordsys: dlc
       * convergence_set: GAU
    * qc_specification:
       * program: psi4
       * driver: SinglepointDriver.deferred
       * method: bp86
       * basis: def2-tzvp
       * keywords: {'maxiter': 500, 'reference': 'uks', 'use_cuest': True, 'scf_properties': ['dipole', 'quadrupole', 'wiberg_lowdin_indices', 'mayer_indices', 'lowdin_charges', 'lowdin_spins', 'mulliken_charges'], 'function_kwargs': {'properties': ['dipole_polarizabilities']}, 'properties_origin': ['COM']}
       * protocols: {'wavefunction': <WavefunctionProtocolEnum.none: 'none'>, 'stdout': True, 'error_correction': {'default_policy': True, 'policies': None}, 'native_files': <NativeFilesProtocolEnum.none: 'none'>}
    * SCF properties:
           * dipole
           * quadrupole
           * wiberg_lowdin_indices
           * mayer_indices
           * lowdin_charges
           * lowdin_spins
           * mulliken_charges