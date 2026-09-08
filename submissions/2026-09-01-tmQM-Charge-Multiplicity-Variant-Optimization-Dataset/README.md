# tmQM Charge Multiplicity Variant Optimization Dataset v0.0

### Description

This dataset was generated starting from the tmQM dataset (release 13Aug2024, https://github.com/uiocompcat/tmQM) containing 108541 unique molecules;
each molecule was evaluated using gfn2-xtb, and then a short MD simulation performed to provide additional configurations of the molecules. Further details
can be found in the hosting repository: https://zenodo.org/records/15059433. This dataset contains 23134 unique transition metal complexes with one Pd, Zn, 
Fe, or Cu, and also only contain elements C, H, P, S, O, N, F, Cl, or Br with charges: {-1,0,+1}. Run with the BP86/def2-TZVP for loose optimizations generating
the properties: 'energy', 'gradient', 'dipole', 'quadrupole', 'wiberg_lowdin_indices', 'mayer_indices', 'lowdin_charges', 'lowdin_spins', 'dipole_polarizabilities',
'mulliken_charges'.

### General Information

- Date: 2026-09-01
- Purpose: BP86/def2-TZVP optimizations for tmQM-derived Pd, Zn, Fe, and Cu complexes across enumerated charge and multiplicity variants with charges of {-1,0,+1}.
- Dataset Type: optimization
- Name: tmQM Charge Multiplicity Variant Optimization Dataset v0.0
- Number of unique molecules: 23,134
- Number of filtered molecules: 0
- Number of Conformers: 23,134
- Number of conformers (min mean max): 1, 1, 1
- Number of Charge/Multiplicity variations per molecule: 103032
- Number of Charge/Multiplicity variations per molecule (min mean max): 1 4 9
- Molecular Weight (min mean max): 95 589 2541
- Metals: {'Cu': 3118, 'Pd': 9362, 'Zn': 6395, 'Fe': 4259}
- Charges: Counter({0.0: 40964, 1.0: 36530, -1.0: 25538})
- Multiplicities: Counter({2.0: 30099, 1.0: 26513, 3.0: 17613, 4.0: 16030, 6.0: 7809, 5.0: 4968})
- Entry name format: tmQM HDF5 label suffixed as `-charge=<total charge>-m<multiplicity>`
- Dataset Submitter: Jennifer A. Clark
- Dataset Curator: Jennifer A. Clark

### QCSubmit generation pipeline

- `generate_dataset.ipynb`: A python notebook which shows how the dataset was prepared from the input files and submitted to QCArchive.
- `enumerate_charge_multiplicity.py`: Helper functions used to enumerate accessible charge and multiplicity variants from the tmQM input structures.
- Input structures were loaded from `tmqm_xtb_dataset_PdZnFeCu_T100_v1.1.hdf5`, derived from the tmQM release dated 13Aug2024 and mirrored in the linked Zenodo record.
- Conformers were not re-enumerated; each submitted entry uses the input coordinates from tmQM for a single charge/multiplicity variant.
- Enumerated stereochemistry: False.
- Enumerated tautomers: False.
- Enumeration settings: total charges limited to {-1, 0, +1}; multiplicities were enumerated from the oxidation-state lookup in `enumerate_charge_multiplicity.py` for the supported metals.

### QCSubmit Manifest

- `README.md`: Submission description, metadata mirror, and QC specification summary
- `generate_dataset.ipynb`: Notebook describing dataset generation and submission
- `enumerate_charge_multiplicity.py`: Charge and multiplicity enumeration helper used by the notebook
- `environment.yml`: Conda environment file to perform this workflow
- `environment_full.yaml`: All installed packages with versions for successful completion of this workflow
- `scaffold.json.bz2`: A compressed json file of the original target dataset
 
### Metadata

* Elements: Pd, Zn, S, P, Fe, O, Cl, N, Br, F, Cu, H, C
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
       * implicit_solvent: none
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