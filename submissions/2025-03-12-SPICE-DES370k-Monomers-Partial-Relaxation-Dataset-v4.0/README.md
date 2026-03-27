# SPICE DES370k Monomers Partial Relaxation Dataset v4.0

## Description
A dataset containing all entries of the `DES370k Monomers` subset of the SPICE dataset, 
optimized for until the energy converges within 1e-4 Ha, at the OpenFF default level of theory (B3LYP-D3BJ/DZVP). 
Intended to be used to expand SPICE to include geometries closer to the QM local minimum.
Detailed description on how the original dataset is generated can be found at https://github.com/openmm/spice-dataset/tree/main/des370k.

## General information
* Name: SPICE DES370k Monomers Partial Relaxation Dataset v4.0
* Number of unique molecules: 376
* Number of conformers: 18700
* Number of conformers (min, mean, max): 50.00, 50.00, 50.00
* Molecular weight (min, mean, max): 16.04, 95.89, 284.78
* Charges: -1.0 0.0 1.0
* Dataset submitter: Alexandra McIsaac
* Dataset generator: Alexandra McIsaac

## QCSubmit generation pipeline
* `generate-dataset.ipynb`: Notebook used to generate dataset

## QCSubmit Manifest
* `dataset.json.bz2`: Compressed dataset ready for submission
* `dataset.pdf`: Visualization of dataset molecules
* `dataset.smi`: Smiles strings for dataset molecules
* `generate-dataset.ipynb`: Notebook used to generate dataset
* `input_environment.yaml`: Environment file used to create Python environment for the notebook
* `input_environment_full.yaml`: Fully-resolved environment used to execute the notebook.

## Metadata
* Elements: {H, C, Cl, I, F, O, S, Br, N, P}
* Spec: default
  * basis: DZVP
  * implicit_solvent: None
  * keywords: {}
  * maxiter: 200
  * method: B3LYP-D3BJ
  * program: psi4
  * SCF properties:
    * dipole
    * quadrupole
    * wiberg_lowdin_indices
    * mayer_indices
* GeometricProcedure
  * program: geometric 
  * coordsys: dlc 
  * enforce: 0.0 
  * epsilon: 1e-05 
  * reset: True 
  * qccnv: False 
  * molcnv: False 
  * check: 0 
  * trust: 0.1 
  * tmax: 0.3 
  * maxiter: 300 
  * convergence_set: CUSTOM 
  * constraints: {} 
  * converge: ['energy', '1e-4', 'grms', '1', 'gmax', '1', 'drms', '1', 'dmax', '1']
