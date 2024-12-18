# OpenFF Sage 2.0.0 Torsion Drive Training Dataset v1.0

## Description

A quantum chemical (QC) dataset curated to train the OpenFF 2.0.0 Sage torsion potentials. This QC dataset with the OpenFF default level of theory, B3LYP-D3BJ/DZVP, is used to benchmark Sage geometries and energetics. These optimized conformer geometries where used to train one dimensional torsional profiles. This Generation 2 dataset increases chemical diversity when compared to Generation 1, which are of value to our industry partners. Large molecules (>20 heavy atoms) were also included, including more flexible molecules and a greater degree of conformational variation which provide intramolecular interactions. This is the complete optimization dataset used for training OpenFF 2.0.0 Sage, consisting of the following datasets: 

'OpenFF Gen 2 Torsion Set 1 Roche', 
'OpenFF Gen 2 Torsion Set 2 Coverage', 'OpenFF Gen 2 Torsion Set 3 Pfizer Discrepancy', 'OpenFF Gen 2 Torsion Set 4 eMolecules  - Discrepancy', 'OpenFF Gen 2 Torsion Set 5 Bayer' and 'OpenFF Gen 2 Torsion Set 6 supplemental 2'. The `HydrogenBondFilter(method='baker-hubbard')` filter was applied, and the following record IDs were dropped due to issues with ForceBalance: 6098580, 2703504, 2703505, 18045478. Further information can be found in the curation scripts for the linked repositories.

## General Information

* Date: 2024-12-17
* Class: OpenFF TorsionDrive Dataset
* Purpose: Complete set of training data for OpenFF 2.0.0 Sage
* Name: OpenFF Sage 2.0.0 Torsion Drive Training Dataset v1.0
* Number of unique molecules: 562
* Number of filtered molecules: 0
* Number of driven torsions: 713
* Number of conformers: 563
* Number of conformers (min, mean, max): 1.00, 1.00, 2.00
* Molecular weight (min, mean, max): 46.07, 224.91, 503.42
* Charges: -1.0, 0.0, 1.0
* Submitter: Jennifer A Clark
* Dataset generator: Hyesu Jang

## QCSubmit Generation Pipeline

* `generate-combined-dataset.py`: A python script which shows how the dataset was prepared from the input files.


## QCSubmit Manifest

* `generate-combined-dataset.py`: 
* `dataset.json.bz2`: The basic dataset ready for submission.
* `dataset.pdf`: A pdf file containing molecule 2D structures.
* `dataset.smi`: SMILES for every molecule in the submission.


## Metadata

* Elements: {Br, C, Cl, F, H, I, N, O, P, S}
* QC Specifications: default
  * basis: DZVP
  * implicit_solvent: None
  * keywords: {}
  * maxiter: 200
  * method: B3LYP-D3BJ
  * program: psi4
  * SCF Properties:
    * dipole
    * quadrupole
    * wiberg_lowdin_indices
    * mayer_indices
