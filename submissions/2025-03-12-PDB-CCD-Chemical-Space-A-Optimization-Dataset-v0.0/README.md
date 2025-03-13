# PDB CCD Chemical Space A: Optimization Dataset v0.0

### Description

This dataset includes molecules from the PDB CCD that contain a single metals: {'Pd', 'Fe', 'Zn', 'Mg', 'Cu', 'Li'}, and the non-metals: {'C', 'H', 'P', 'S', 'O', 'N', 'F', 'Cl', 'Br'}, with a complex charge of {-1,0,+1}. These complexes are minimized using gfn-xtb.

### General Information

- Date: 2025-03-12
- Class: OpenFF Optimization Dataset
- Purpose: GFN2-xTB Conformers for single metal complexes with Pd, Fe, Zn, Cu, Mg, Li and change of {-1,0,+1}
- Dataset Type: optimization
- Name: PDB CCD Chemical Space A: Optimization Dataset v0.0
- Number of unique molecules:   137
- Number of filtered molecules: 0
- Number of Conformers: 136
- Number of conformers (min mean max): 1 1 2
- Molecular Weight (min mean max): 71, 565, 1177
- Set of charges: -1.0, 0.0, 1.0
- Dataset Submitter: Jennifer A. Clark
- Dataset Curator: Jennifer A. Clark

### QCSubmit generation pipeline

- `generate-dataset.py`: A python script which shows how the dataset was prepared from the input files.
- `output.txt`: A text file containing the printed output of `generate-combined-dataset.py`.

### QCSubmit Manifest

- `generate-combined-dataset.py`
- `sanitize_tm.py`: A python module to perform same corrective actions on RDKit molecules
- `CS-A_primary.json`: A file containing the subset of the PDB CCD that we will include with our dataset 
- `dataset.pdf`: A pdf file containing molecule 2D structures.
- `dataset.smi`: SMILES for every molecule in the submission.
 
### Metadata

* Elements: {'F', 'Fe', 'Pd', 'C', 'Mg', 'O', 'N', 'P', 'Cl', 'Zn', 'S', 'Cu', 'H'}
* QC Specifications: default
  * program: xtb
  * method: gfn2xtb
  * basis: None
  * implicit_solvent: None
  * keywords: {}
  * maxiter: 200
  * SCF Properties:
    * dipole
    * quadrupole
    * wiberg_lowdin_indices
    * mayer_indices
* GeometricProcedure
  * tmax: 0.3,
  * check: 0,
  * qccnv: False,
  * reset: True,
  * trust: 0.1,
  * molcnv: False,
  * enforce: 0.0,
  * epsilon: 1e-05,
  * maxiter: 300,
  * coordsys: 'dlc',
  * constraints: {},
  * convergence_set: 'GAU'