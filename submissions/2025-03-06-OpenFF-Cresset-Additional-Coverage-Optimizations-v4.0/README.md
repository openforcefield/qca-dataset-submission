# OpenFF Cresset Additional Coverage Optimizations v4.0

## Description

Molecules contributed by Cresset to address lack of data coverage for particular torsion drives. This is the optimization dataset counterpart to [the torsiondrive dataset](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2025-02-12-OpenFF-Cresset-Additional-Coverage-TorsionDrives-v4.0). The input SMILES files are the same files as in the torsiondrive counterpart.

These molecules are contributed to address five failure cases:

1) Sage 2.1 t48a is fit only to one molecule,
where the complementary torsion contributes most of the profile

2) Sage 2.1 t17 may benefit from splitting and a different shape for non-symmetric rings;
the current n=3 shape sums to a constant

3) Sage 2.1 t19 is mostly trained to terminal methyls.
It has a functional form where the n=1 term dominates unexpectedly,
instead of the more expected equal n=3 contributions.
More data with non-terminal methyls is added

4) Sage 2.1 t18 covers amide-adjacent torsions but is not trained to many.

5) Sage 2.1 t105 covers an O linker with an sp2 or sp3 terminus.
While the sp3 profiles match the QM well, the sp2 profiles look too stiff.
This adds more data.

## General Information

* Date: 2025-03-06
* Class: OpenFF Optimization Dataset
* Purpose: Improve torsiondrive coverage with optimizations
* Name: OpenFF Cresset Additional Coverage Optimizations v4.0
* Number of unique molecules: 70
* Number of conformers: 393
* Number of conformers (min, mean, max): 1.00, 5.61, 10.00
* Molecular weight (min, mean, max): 58.08, 144.98, 280.75
* Charges: [-1.0, 0.0, 1.0]
* Dataset generator: Matthew Habgood
* Dataset submitter: Lily Wang

## QCSubmit Generation Pipeline

* `generate-dataset.ipynb`: This notebooks shows how the dataset was prepared from the
  input files in `inputs`. The input files haven't been added to this repo, but are the same as those in [the torsiondrive dataset](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2025-02-12-OpenFF-Cresset-Additional-Coverage-TorsionDrives-v4.0).

## QCSubmit Manifest

### Input files
* `generate-dataset.ipynb`: Notebook describing dataset generation and submission
* `input-env.yaml`: Environment file used to create Python environment for the script
* `full-env.yaml`: Fully-resolved environment used to execute the script

### Output Files
* `dataset.json.bz2`: Compressed dataset ready for submission
* `dataset.pdf`: Visualization of dataset molecules
* `output.smi`: SMILES strings for dataset molecules

## Metadata
* Elements: {O, C, F, S, H, N, Br, Cl}
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

