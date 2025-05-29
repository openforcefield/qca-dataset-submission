# OpenFF SMIRNOFF Sage 2.0.0

## Description

A quantum chemical (QC) dataset of optimization targets was generated at the OpenFF default level of theory, B3LYP-D3BJ/DZVP, and curated to train the 
parameters of the [OpenFF 2.0.0 Sage](https://github.com/openforcefield/openff-sage) forcefield. This Generation 2 dataset increases chemical diversity when 
compared to Generation 1, which are of value to our industry partners. Large molecules (>20 heavy atoms) were also included, including more flexible molecules 
and a greater degree of conformational variation which provide intramolecular interactions.

Further information can be found in the curation scripts for the linked repositories in the following sections.

### Optimization Datasets

This is the complete optimization dataset used for training OpenFF 2.0.0 Sage, consisting of molecules from following datasets:

 - [OpenFF Gen 2 Opt Set 1 Roche](https://github.com/openforcefield/qca-dataset-submission/tree/0e6e6da930118e2a2d6402b93c3e3e93830600cc/submissions/2020-03-20-OpenFF-Gen-2-Optimization-Set-1-Roche)
 - [OpenFF Gen 2 Opt Set 2 Coverage](https://github.com/openforcefield/qca-dataset-submission/tree/0e6e6da930118e2a2d6402b93c3e3e93830600cc/submissions/2020-03-20-OpenFF-Gen-2-Optimization-Set-2-Coverage)
 - [OpenFF Gen 2 Opt Set 3 Pfizer Discrepancy](https://github.com/openforcefield/qca-dataset-submission/tree/0e6e6da930118e2a2d6402b93c3e3e93830600cc/submissions/2020-03-20-OpenFF-Gen-2-Optimization-Set-3-Pfizer-Discrepancy)
 - [OpenFF Gen 2 Opt Set 4 eMolecules Discrepancy](https://github.com/openforcefield/qca-dataset-submission/tree/0e6e6da930118e2a2d6402b93c3e3e93830600cc/submissions/2020-03-20-OpenFF-Gen-2-Optimization-Set-4-eMolecules-Discrepancy)
 - [OpenFF Gen 2 Opt Set 5 Bayer](https://github.com/openforcefield/qca-dataset-submission/tree/0e6e6da930118e2a2d6402b93c3e3e93830600cc/submissions/2020-03-20-OpenFF-Gen-2-Optimization-Set-5-Bayer)

The following filters were applied to those datasets:

 - `RecordStatusFilter(status=RecordStatusEnum.complete)`
 - `ConnectivityFilter(tolerance=1.2)`
 - `UndefinedStereoFilter()`
 - `ElementFilter(allowed_elements=["H", "C", "N", "O", "S", "P", "F", "Cl", "Br", "I"])`
 - `ConformerRMSDFilter(max_conformers=10)`

### Torsiondrive Datasets

This is the complete TorsionDrive dataset 
used for training OpenFF 2.0.0 Sage, consisting of data drawn from the following datasets: 

- [OpenFF Gen 2 Torsion Set 1 Roche 2](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2020-03-23-OpenFF-Gen-2-Torsion-Set-1-Roche-2)
- [OpenFF Gen 2 Torsion Set 2 Coverage 2](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2020-03-23-OpenFF-Gen-2-Torsion-Set-2-Coverage-2)
- [OpenFF Gen 2 Torsion Set 3 Pfizer Discrepancy 2](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2020-03-23-OpenFF-Gen-2-Torsion-Set-3-Pfizer-Discrepancy-2)
- [OpenFF Gen 2 Torsion Set 4 eMolecules Discrepancy 2](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2020-03-23-OpenFF-Gen-2-Torsion-Set-4-eMolecules-Discrepancy-2)
- [OpenFF Gen 2 Torsion Set 5 Bayer 2](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2020-03-26-OpenFF-Gen-2-Torsion-Set-5-Bayer-2)
- [OpenFF Gen 2 Torsion Set 6 Supplemental 2](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2020-03-26-OpenFF-Gen-2-Torsion-Set-6-supplemental-2)

The `HydrogenBondFilter(method='baker-hubbard')` filter was applied, and the 
following record IDs were dropped due to issues with ForceBalance: 6098580, 2703504, 
2703505, 18045478.

### General Information

- Date: 2025-05-29
- Name: OpenFF SMIRNOFF Sage 2.0.0
- Purpose: Complete set of training data for OpenFF 2.0.0 Sage
- Dataset Submitter: Jennifer A. Clark
- Dataset Curator: Simon Boothroyd
- Dataset Generator: Jessica Maat and Hyesu Jang
---
- Class: OpenFF Optimization Dataset
- Dataset Type: optimization
- Number of unique molecules:   1039
- Number of filtered molecules: 0 
- Number of conformers:         3663
- Number of conformers (min mean max): 1.00, 3.53, 10.00
- Mean molecular weight: 261.37
- Max molecular weight: 544.64
- Set of charges: -2.0, -1.0, 0.0, 1.0
---
* Class: OpenFF TorsionDrive Dataset
- Dataset Type: torsiondrive
* Number of unique molecules: 562
* Number of filtered molecules: 0
* Number of driven torsions: 713
* Number of conformers: 563
* Number of conformers (min, mean, max): 1, 1, 2
* Molecular weight (min, mean, max): 46.07, 224.91, 503.42
* Charges: -1.0, 0.0, 1.0

### QCSubmit generation pipeline

* `opt_main.ipynb`: A jupyter notebook which shows how the optimization dataset was prepared from the input files.
* `td_main.ipynb`: A jupyter notebook which shows how the torsiondrive dataset was prepared from the input files.

## Manifest

* `opt_main.ipynb`: A jupyter notebook which shows how the optimization dataset was prepared from the input files.
* `td_main.ipynb`: A jupyter notebook which shows how the torsiondrive dataset was prepared from the input files.
* `opt_ds_info.json`: Metadata information for optimization dataset record imported by `opt_main.ipynb`
* `td_ds_info.json`: Metadata information for torsiondrive dataset record imported by `td_main.ipynb`
* `scaffold_opt.json.bz2`: The basic optimization dataset information read by qcportal.external.scaffold.
* `scaffold_td.json.bz2`: The basic torsiondrive dataset information read by qcportal.external.scaffold.
* `environment.yaml`: A file to reproduce the conda env used to generate this dataset.
* `environment_full.yaml`: A file to reproduce the conda env used to generate this dataset.

## Optimization Metadata

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

## Torsiondrive Metadata

* Elements: {Br, C, Cl, F, H, I, N, O, P, S}
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