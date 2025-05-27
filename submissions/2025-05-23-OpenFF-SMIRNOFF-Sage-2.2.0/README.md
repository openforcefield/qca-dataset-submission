# OpenFF SMIRNOFF Sage 2.2.0 v1.0

## Description

A quantum chemical (QC) optimization and torsiondrive datasets generated at the OpenFF default level of theory, B3LYP-D3BJ/DZVP, and
 curated to train parameters in [OpenFF 2.2.0 Sage](https://github.com/openforcefield/sage-2.2.0/) with improved small ring
  internal angles and sulfamide geometries. 

### Optimization Datasets

Targets were curated from the following datasets:

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
 - `HydrogenBondFilter(method='baker-hubbard')`
 - `ConnectivityFilter(tolerance=1.2)`
 - `UnperceivableStereoFilter()`
 - `ChargeCheckFilter()`
  
### Torsiondrive Datasets

Targets were curated from the following datasets:

 - [OpenFF Gen 2 Torsion Set 1 Roche 2](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2020-03-23-OpenFF-Gen-2-Torsion-Set-1-Roche-2)
 - [OpenFF Gen 2 Torsion Set 2 Coverage 2](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2020-03-23-OpenFF-Gen-2-Torsion-Set-2-Coverage-2)
 - [OpenFF Gen 2 Torsion Set 3 Pfizer Discrepancy 2](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2020-03-23-OpenFF-Gen-2-Torsion-Set-3-Pfizer-Discrepancy-2)
 - [OpenFF Gen 2 Torsion Set 4 eMolecules Discrepancy 2](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2020-03-23-OpenFF-Gen-2-Torsion-Set-4-eMolecules-Discrepancy-2)
 - [OpenFF Gen 2 Torsion Set 5 Bayer 2](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2020-03-26-OpenFF-Gen-2-Torsion-Set-5-Bayer-2)
 - [OpenFF Gen 2 Torsion Set 6 Supplemental 2](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2020-03-26-OpenFF-Gen-2-Torsion-Set-6-supplemental-2)
 - [SMIRNOFF Coverage Torsion Set 1](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2019-07-01-smirnoff99Frost-coverage-torsion)
 - [OpenFF Group1 Torsions](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2019-05-01-OpenFF-Group1-Torsions)
 - [OpenFF Group1 Torsions 2](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2020-01-31-OpenFF-Group1-Torsions-2)
 - [OpenFF Group1 Torsions 3](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2020-02-10-OpenFF-Group1-Torsions-3)
 - [Pfizer Discrepancy Torsion Dataset 1](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2019-09-07-Pfizer-discrepancy-torsion-dataset-1)
 - [OpenFF Gen3 Torsion Set v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2021-04-09-OpenFF-Gen3-Torsion-Set-v1.0)
 - [OpenFF Amide Torsion Set v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2021-03-23-OpenFF-Amide-Torsion-Set-v1.0)
 - [OpenFF WBO Conjugated Series v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2021-01-25-OpenFF-Conjugated-Series)
 - [OpenFF DANCE 1 eMolecules t142 v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2020-06-01-DANCE-1-eMolecules-t142-selected)
 
 These combined datasets were filtered with:
 
 - `ElementFilter(allowed_elements=['H', 'C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br'])`
 - `RecordStatusFilter(status=RecordStatusEnum.complete)`
 - `HydrogenBondFilter(method='baker-hubbard')`
 - `ConnectivityFilter(tolerance=1.2)`
 - `UnperceivableStereoFilter()`
 - `ChargeCheckFilter()`

## General Information

* Date: 2025-05-23
* Purpose: Complete set of training data for OpenFF 2.2.0 Sage
* Name: OpenFF SMIRNOFF Sage 2.2.0
* Dataset submitter: Jennifer A Clark
* Dataset curator: Pavan Behara

* Class: OpenFF Optimization Dataset
* Dataset Type: optimization
* Number of unique molecules: 1691
* Number of filtered molecules: 0
* Number of conformers: 5126
* Number of conformers (min, mean, max): 1.00, 3.03, 12.00
* Molecular weight (min, mean, max): 16.04, 236.01, 544.64
* Charges: -3.0, -2.0, -1.0, 0.0, 1.0, 2.0
* Dataset generator: Chaya Stern, Hyesu Jang, Jessica Maat, and Pavan Behara

* Class: OpenFF TorsionDrive Dataset
* Dataset Type: torsiondrive
* Number of unique molecules: 956
* Number of filtered molecules: 0
* Number of driven torsions: 1290
* Number of conformers: 982
* Number of conformers (min, mean, max): 1, 1, 3
* Molecular weight (min, mean, max): 46.07, 185.48, 503.42
* Charges: -1.0, 0.0, 1.0
* Dataset generator: Simon Boothroyd, John Chodera, Trevor Gokey, Hyesu Jang, Yudong Qiu, Bryon Tjanaka

## Generation Pipeline

* `opt_main.ipynb`: A jupyter notebook which shows how the optimization dataset was prepared from the input files.
* `td_main.ipynb`: A jupyter notebook which shows how the torsiondrive dataset was prepared from the input files.

## Manifest

* `opt_main.ipynb`: A jupyter notebook which shows how the optimization dataset was prepared from the input files.
* `td_main.ipynb`: A jupyter notebook which shows how the torsiondrive dataset was prepared from the input files.
* `opt_ds_info.json`: Metadata information for optimization dataset records imported by `opt_main.ipynb`
* `td_ds_info.json`: Metadata information for torsiondrive dataset records imported by `td_main.ipynb`
* `scaffold_opt.json.bz2`: The basic optimization dataset information read by qcportal.external.scaffold.
* `scaffold_td.json.bz2`: The basic torsiondrive dataset information read by qcportal.external.scaffold.
* `environment.yaml`: A file to reproduce the conda env used to generate this dataset.
* `environment_full.yaml`: A file to reproduce the conda env used to generate this dataset.

## Optimization Metadata

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

## Torsiondrive Metadata

* Elements: {Br, C, Cl, F, H, N, O, P, S}
* Program: torsiondrive
* Optimization Specification: geometric
* QC Specifications: default
  * program: psi4
  * method: B3LYP-D3BJ
  * basis: DZVP
  * implicit_solvent: None
  * keywords: {'maxiter': 200, 'scf_properties': ['dipole', 'quadrupole', 'wiberg_lowdin_indices', 'mayer_indices']}
  * protocols: {}
  * maxiter: 200
  * SCF Properties:
    * dipole
    * quadrupole
    * wiberg_lowdin_indices
    * mayer_indices