# OpenFF Sage 2.1.0 Optimization Training Dataset v1.0

## Description

A quantum chemical (QC) dataset of single point energies and properties was curated to train the 
[OpenFF 2.1.0 Sage](https://github.com/openforcefield/sage-2.1.0/) forcefield, with reparametrized Lennard-Jones (LJ) 
and valence parameters. This QC dataset with the OpenFF default level of theory, B3LYP-D3BJ/DZVP, is used to train Sage parameters. 
This is the complete optimization dataset includes calculations used for training OpenFF 2.0.0 Sage in addition to other, consisting 
of the following datasets:

 - ['OpenFF Gen 2 Opt Set 1 Roche'](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2020-03-20-OpenFF-Gen-2-Optimization-Set-1-Roche)
 - ['OpenFF Gen 2 Opt Set 2 Coverage'](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2020-03-20-OpenFF-Gen-2-Optimization-Set-2-Coverage)
 - ['OpenFF Gen 2 Opt Set 3 Pfizer Discrepancy'](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2020-03-20-OpenFF-Gen-2-Optimization-Set-3-Pfizer-Discrepancy)
 - ['OpenFF Gen 2 Opt Set 4 eMolecules Discrepancy'](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2020-03-12-OpenFF-Gen-2-Torsion-Set-4-eMolecules-Discrepancy)
 - ['OpenFF Gen 2 Opt Set 5 Bayer'](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2020-03-20-OpenFF-Gen-2-Optimization-Set-5-Bayer)
 - ['OpenFF Gen2 Optimization Dataset Protomers v1.0'](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2021-12-21-OpenFF-Gen2-Optimization-Set-Protomers)
 - ['OpenFF Iodine Chemistry Optimization Dataset v1.0'](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2022-07-27-OpenFF-iodine-optimization-set)
 - ['OpenFF Optimization Set 1'](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2019-05-16-Roche-Optimization_Set)
 - ['SMIRNOFF Coverage Set 1'](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2020-03-20-OpenFF-Gen-2-Optimization-Set-2-Coverage)
 - ['OpenFF Aniline Para Opt v1.0'](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2021-04-02-OpenFF-Aniline-Para-Opt-v1.0)

These combined datasets were filtered with:

- `ElementFilter(allowed_elements=['H', 'C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br'])`
- `RecordStatusFilter(status=RecordStatusEnum.complete)`
- `ConnectivityFilter(tolerance=1.2)`
- `UnperceivableStereoFilter()`
- `ConformerRMSDFilter(max_conformers=12)`
- `ChargeCheckFilter()`

 Further information can be found in the curation scripts for the linked repositories.

## General Information

* Date: 2025-05-22
* Class: OpenFF Optimization Dataset
* Purpose: Complete set of optimization training data for OpenFF 2.1.0 Sage
* Dataset Type: optimization
* Name: OpenFF Sage 2.1.0 Optimization Training Dataset v1.0
* Number of unique molecules: 1701
* Number of filtered molecules: 0
* Number of conformers: 5580
* Number of conformers (min, mean, max): 1.00, 3.17, 17.00
* Molecular weight (min, mean, max): 16.04, 235.40, 544.64
* Charges: -3.0, -2.0, -1.0, 0.0, 1.0, 2.0
* Submitter: Jennifer A Clark
* Dataset curator: Pavan Behara
* Dataset generator: Chaya Stern, Hyesu Jang, Jessica Maat, and Pavan Behara

## Generation Pipeline

* `main.ipynb`: A jupyter notebook which shows how the dataset was prepared from the input files.

## Manifest

* `main.ipynb`: A jupyter notebook which shows how the dataset was prepared from the input files.
* `ds_info.json`: Metadata information for dataset record imported by `main.ipynb`
* `scaffold.json.bz2`: The basic dataset information read by qcportal.external.scaffold.
* `conda_env.yaml`: A file to reproduce the conda env used to generate this dataset.
* `conda_env_full.yaml`: A file to reproduce the conda env used to generate this dataset.


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
