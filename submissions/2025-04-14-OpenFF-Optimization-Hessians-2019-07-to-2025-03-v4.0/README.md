# OpenFF Optimization Hessians 2019-07 to 2025-03 v4.0

## Description

Hessian single points for the final molecules in the OpenFF datasets listed below at the B3LYP-D3BJ/DZVP level of theory. These are used for calculating MSM starting points in force field fits. The molecules here include the S, H, O, Br, F, N, P, Cl, I, C elements and the charge states {-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0}. They range from 16-1425 Da (mean 224) and 4-99 heavy atoms. Records were filtered for successful completion, no connectivity changes, a non-2D structure (where all Z-coordinates are 0), and whether RDKit can parse the molecule with valid valence. The datasets included are:
 - OpenFF Optimization Set 1
 - SMIRNOFF Coverage Set 1
 - OpenFF VEHICLe Set 1
 - OpenFF Discrepancy Benchmark 1
 - OpenFF Ehrman Informative Optimization v0.2
 - Pfizer discrepancy optimization dataset 1
 - FDA optimization dataset 1
 - Kinase Inhibitors: WBO Distributions
 - OpenFF Gen 2 Opt Set 1 Roche
 - OpenFF Gen 2 Opt Set 2 Coverage
 - OpenFF Gen 2 Opt Set 3 Pfizer Discrepancy
 - OpenFF Gen 2 Opt Set 4 eMolecules Discrepancy
 - OpenFF Gen 2 Opt Set 5 Bayer
 - OpenFF Sandbox CHO PhAlkEthOH v1.0
 - OpenFF Industry Benchmark Season 1 v1.1
 - OpenFF Gen2 Optimization Dataset Protomers v1.0
 - OpenFF Protein Capped 1-mers 3-mers Optimization Dataset v1.0
 - OpenFF Iodine Chemistry Optimization Dataset v1.0
 - XtalPi Shared Fragments OptimizationDataset v1.0
 - XtalPi 20-percent Fragments OptimizationDataset v1.0
 - OpenFF Torsion Benchmark Supplement v1.0
 - OpenFF Torsion Multiplicity Optimization Training Coverage Supplement v1.0
 - OpenFF Torsion Multiplicity Optimization Benchmarking Coverage Supplement v1.0
 - OpenFF Iodine Fragment Opt v1.0
 - OpenFF Sulfur Optimization Training Coverage Supplement v1.0
 - OpenFF Sulfur Optimization Benchmarking Coverage Supplement v1.0
 - OpenFF Lipid Optimization Training Supplement v1.0
 - OpenFF Lipid Optimization Benchmark Supplement v1.0
 - OpenFF Cresset Additional Coverage Optimizations v4.0
 - OpenFF Protein PDB 4-mers v4.0

Any molecules in the below datasets had all molecules containing iodine filtered out, as those records were problematic.
 - OpenFF Discrepancy Benchmark 1
 - OpenFF Gen 2 Opt Set 2 Coverage
 - OpenFF Gen 2 Opt Set 4 eMolecules Discrepancy
 - SMIRNOFF Coverage Set 1
 - OpenFF Ehrman Informative Optimization v0.2
 - FDA optimization dataset 1
 - Kinase Inhibitors: WBO Distributions
 - OpenFF Gen 2 Torsion Set 2 Coverage 2
 - OpenFF Gen 2 Torsion Set 4 eMolecules Discrepancy 2

## General information

* Date: 2025-04-14
* Class: OpenFF Basic Dataset
* Purpose: Hessian dataset generation for MSM parameters, from OpenFF datasets from 2019-07 to 2025-03
* Name: OpenFF Optimization Hessians 2019-07 to 2025-03 v4.0
* Number of unique molecules: 63446
* Number of conformers: 297934
* Number of conformers (min, mean, max): 1, 4.83, 275
* Molecular weight (min, mean, max): 16.04, 224.36, 1425.34
* Charges: [-1.0, 0.0, 1.0]
* Dataset generator: Lily Wang
* Dataset submitter: Lily Wang


## QCSubmit generation pipeline

* `generate-dataset.ipynb`: This notebook shows how the dataset was prepared from the input files.


## QCSubmit Manifest

* `optimizations.json.tar.gz`: Compressed combined input OptimizationResultCollection
* `dataset.json.bz2`: Compressed dataset ready for submission
* `dataset.pdf.tar.gz`: Compressed visualization of dataset molecules
* `dataset.smi`: Smiles strings for dataset molecules
* `generate-dataset.ipynb`: Notebook describing dataset generation and submission
* `input-env.yaml`: Environment file used to create Python environment for the notebook
* `full-env.yaml`: Fully-resolved environment used to execute the notebook.


## Metadata

* elements: {S, H, O, Br, F, N, P, Cl, I, C}
* unique molecules: 63446
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

