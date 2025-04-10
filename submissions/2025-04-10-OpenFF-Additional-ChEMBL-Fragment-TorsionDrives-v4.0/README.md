# OpenFF Additional Generated ChEMBL TorsionDrives 4.0

## Description

This dataset generates additional torsiondrives from molecules for torsions in Sage 2.2.1 with low or no existing coverage in QCArchive.
Low coverage torsions are identified from the [profile-qc-data](https://github.com/lilyminium/profile-qc-data) repo.
Rare torsions that run through a ring were not included.
Input files are obtained from the process in the [fragment-chembl-data](https://github.com/lilyminium/fragment-chembl-data/commit/d94c9fc4be6945c4e680b2bb5b31a2693c4b772b) repo and copied to `inputs/`.

In short:

- ChEMBL molecules were split into "elementary" fragments without rotatable bonds
- Elementary fragments were combined with a single bond
- For each torsion, a pool of up to 5000 molecules were initially selected after sorting for low molecular weight
- Up to 250 molecules were selected from this pool by maximising chemical diversity, using the Tanimoto distance of the Morgan fingerprints
- Up to 10 molecules were selected from the second pool by maximising the diversity of coupled torsions through the central bond



## General Information

* Date: 2025-04-10
* Class: OpenFF TorsionDrive Dataset
* Purpose: Improve torsiondrive coverage of low-coverage torsions with automatically selected molecules from ChEMBL fragments
* Dataset name: OpenFF Additional Generated ChEMBL TorsionDrives 4.0
* Number of unique molecules: 270
* Number of driven torsions: 275
* Number of filtered molecules: 17
* Number of conformers: 352
* Number of conformers per molecule (min, mean, max): 1, 1.28, 4
* Mean molecular weight: 124.98
* Max molecular weight: 312.99
* Charges: [-1.0, 0.0, 1.0, 2.0]
* Dataset generator: Lily Wang
* Dataset submitter: Lily Wang

## QCSubmit Generation Pipeline

* `generate-dataset.ipynb`: This notebooks shows how the dataset was prepared from the
  input files in `inputs`.

## QCSubmit Manifest

### Input files
* `inputs/generate-smiles.ipynb`: Notebook showing how molecules were generated and checked for each torsion.
* `generate-dataset.ipynb`: Notebook describing dataset generation and submission
* `input-env.yaml`: Environment file used to create Python environment for the script
* `full-env.yaml`: Fully-resolved environment used to execute the script

### Output Files
* `dataset.json.bz2`: Compressed dataset ready for submission
* `dataset.pdf`: Visualization of dataset molecules
* `output.smi`: SMILES strings for dataset molecules
* `inputs/selected-torsions-*.smi`: Selected SMILES for particular torsion
* `inputs/*_molecules.pdf`: Generated molecules for particular torsion


## Metadata
* Elements: {O, Cl, Br, C, I, P, F, H, N, S}
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

