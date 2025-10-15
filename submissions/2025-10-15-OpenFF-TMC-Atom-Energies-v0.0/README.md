# OpenFF TMC Atom Energies v0.0

### Description

This dataset generates the per atom energies for use in calculating formation energies. The model chemistries chosen are: B3LYP-D3BJ/DZVP (OpenFF default), BP86/def2-TZVP (used in OpenFF tmQM), wB97M-V/def2-TZVPD (OMol25), wB97M-D3BJ/def2-TZVPPD (SPICE). Elements include C, H, P, S, O, N, F, Cl, B, Li, Na, K, Mg, Ca, Br, Pd, Fe, Zn, Cu, Rh, Ir, Pt, Ni, Cr, Ag, Ti. Each atom is run at a charge of -1, 0, 1 for nonmetals, 0, 1 for alkali metals, 0, 1, 2 for alkaline earth metals, and -2, -1, 0, 1, 2 for transition metals. Each element / charge combination is run for valid multiplicities up to 6 with a reference of uks. The molecular weight varies with elemental values from 1.008 for H to 195.085 for Pt. Entry names for this dataset are: `"{element symbol}_{charge}_{multiplicity}"`

### General Information

- Date: 2025-10-15
- Purpose: Element energies at various multiplicities, and formal charges of +1, 0, and -1. Elements include C, H, P, S, O, N, F, Cl, B, Li, Na, K, Mg, Ca, Br, Pd, Fe, Zn, Cu, Rh, Ir, Pt, Ni, Cr, Ag, Ti at the following model chemistries: B3LYP-D3BJ/DZVP, BP86/def2-TZVP, wB97M-V/def2-TZVPD, wB97M-D3BJ/def2-TZVPPD.
- Dataset Type: singlepoint
- Name: OpenFF TMC Atom Energies v0.0
- Number of elements: 25
- Number of filtered molecules: 0
- Number of Conformers: 25
- Number of conformers (min mean max): 1, 1, 1
- Molecular Weight (min max): 1.008, 195.085
- Set of charges: -2, -1, 0, 1, 2
- Dataset Submitter: Jennifer A. Clark
- Dataset Curator: Jennifer A. Clark

### QCSubmit generation pipeline

- `generate-dataset.ipynb`: A python notebook which shows how the dataset was prepared from the input files.

### QCSubmit Manifest

- `generate-dataset.ipynb`
- `environment.yml`: Conda environment file to perform this workflow
- `environment_full.yml`: All installed packages with versions for successful completion of this workflow
- `scaffold.json.bz2`: A compressed json file of the target dataset
 
### Metadata

* Elements: {"C", "H", "P", "S", "O", "N", "F", "Cl", "B", "Li", "Na", "K", "Mg", "Ca", "Br", "Pd", "Fe", "Zn", "Cu", "Rh", "Ir", "Pt", "Ni", "Cr", "Ag", "Ti"}

* QC Specifications: B3LYP-D3BJ/DZVP
  * excluded elements: 'Ir', 'Pt'
  * program: psi4
  * method: B3LYP-D3BJ
  * basis: DZVP
  * reference: 'uks'
  * driver: energy
  * implicit_solvent: None
  * keywords: {}
  * maxiter: 500

* QC Specifications: BP86/def2-TZVP
  * program: psi4
  * method: BP86
  * basis: def2-TZVP
  * reference: 'uks'
  * driver: energy
  * implicit_solvent: None
  * keywords: {}
  * maxiter: 500

* QC Specifications: wB97M-V/def2-TZVPD
  * program: psi4
  * method: wB97M-V
  * basis: def2-TZVPD
  * reference: 'uks'
  * driver: energy
  * implicit_solvent: None
  * keywords: {}
  * maxiter: 500

* QC Specifications: wB97M-D3BJ/def2-TZVPPD
  * program: psi4
  * method: wB97M-D3BJ
  * basis: def2-TZVPPD
  * reference: 'uks'
  * driver: gradient
  * implicit_solvent: None
  * keywords: {}
  * maxiter: 500

* SCF Properties:
  * dipole_polarizabilities