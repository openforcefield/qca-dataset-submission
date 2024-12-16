# OpenFF Sage 2.0.0 Training Optimization v1.0

### Description

A quantum chemical (QC) dataset curated to train [OpenFF 2.0.0 Sage](https://github.com/openforcefield/openff-sage) forcefield, with reparametrized Lennard-Jones (LJ) and valence parameters, the latter relevant to this dataset. This QC dataset with the OpenFF default level of theory, B3LYP-D3BJ/DZVP, is used to benchmark Sage geometries and energetics. These optimized conformer geometries where used in conjunction with the QC dataset used to train one dimensional torsional profiles. This Generation 2 dataset increases chemical diversity when compared to Generation 1, which are of value to our industry partners. Large molecules (>20 heavy atoms) were also included, including more flexible molecules and a greater degree of conformational variation which provide intramolecular interactions.

This is the complete optimization dataset used for training OpenFF 2.0.0 Sage, consisting of the following datasets:

 - [OpenFF Gen 2 Opt Set 1 Roche](https://github.com/openforcefield/qca-dataset-submission/tree/0e6e6da930118e2a2d6402b93c3e3e93830600cc/submissions/2020-03-20-OpenFF-Gen-2-Optimization-Set-1-Roche)
 - [OpenFF Gen 2 Opt Set 2 Coverage](https://github.com/openforcefield/qca-dataset-submission/tree/0e6e6da930118e2a2d6402b93c3e3e93830600cc/submissions/2020-03-20-OpenFF-Gen-2-Optimization-Set-2-Coverage)
 - [OpenFF Gen 2 Opt Set 3 Pfizer Discrepancy](https://github.com/openforcefield/qca-dataset-submission/tree/0e6e6da930118e2a2d6402b93c3e3e93830600cc/submissions/2020-03-20-OpenFF-Gen-2-Optimization-Set-3-Pfizer-Discrepancy)
 - [OpenFF Gen 2 Opt Set 4 eMolecules  - Discrepancy](https://github.com/openforcefield/qca-dataset-submission/tree/0e6e6da930118e2a2d6402b93c3e3e93830600cc/submissions/2020-03-20-OpenFF-Gen-2-Optimization-Set-4-eMolecules-Discrepancy)
 - [OpenFF Gen 2 Opt Set 5 Bayer](https://github.com/openforcefield/qca-dataset-submission/tree/0e6e6da930118e2a2d6402b93c3e3e93830600cc/submissions/2020-03-20-OpenFF-Gen-2-Optimization-Set-5-Bayer)

The following filters were applied to those datasets:

 - `RecordStatusFilter(status=RecordStatusEnum.complete)`
 - `ConnectivityFilter(tolerance=1.2)`
 - `UndefinedStereoFilter()`
 - `ElementFilter(allowed_elements=["H", "C", "N", "O", "S", "P", "F", "Cl", "Br", "I"])`
 - `ConformerRMSDFilter(max_conformers=10)`

 Further information can be found in the curation scripts for the linked repositories.

### General Information

- Date: 2024-12-12
- Class: OpenFF Optimization Dataset
- Purpose: Complete set of training data for OpenFF 2.0.0 Sage
- Dataset Type: optimization
- Name: OpenFF Sage 2.0.0 Training Optimization Dataset v1.0
- Number of unique molecules:   1025
- Number of filtered molecules: 0 
- Number of conformers:         3663
- Number of conformers (min mean max): 1.00, 3.53, 10.00
- Mean molecular weight: 261.38
- Max molecular weight: 544.64
- Set of charges: -2.0, -1.0, 0.0, 1.0
- Dataset Submitter: Jennifer A. Clark
- Dataset Curator: Simon Boothroyd
- Dataset Generator: Hyesu Jang

### QCSubmit generation pipeline

- `generate-combined-dataset.py`: A python script which shows how the dataset was prepared from the input files.
- `output.txt`: A text file containing the printed output of `generate-combined-dataset.py`.

### QCSubmit Manifest

- `generate-combined-dataset.py`
- `dataset.json.bz2`: The basic dataset ready for submission.
- `dataset.pdf`: A pdf file containing molecule 2D structures.
- `dataset.smi`: SMILES for every molecule in the submission.
 
### Metadata

* Elements: {F, I, N, C, P, Cl, S, Br, O, H}
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