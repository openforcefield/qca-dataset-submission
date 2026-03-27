# OpenFF SMIRNOFF Sage 2.3.0

## Description

A quantum chemical (QC) optimization and torsiondrive datasets generated at the OpenFF default level of theory, B3LYP-D3BJ/DZVP, and curated to train parameters in [OpenFF 2.3.0 Sage](https://github.com/openforcefield/ash-sage-rc2) with NAGL partial charge model AshGC v1.0.

Datasets available on [Zenodo](doi.org/10.5281/zenodo.18436107), DOI: 10.5281/zenodo.18436107.

### Optimization Datasets

The optimization records were selected to maximize chemical diversity using a selection of record IDs listed in the [Sage 2.3.0 repository](https://raw.githubusercontent.com/openforcefield/ash-sage-rc2/32345dddeb6cb249367059fd99607ac2950a5c86/03_fit-valence/02_curate-data/output/optimizations-single-v3.json). These records came from the following datasets:

- [OpenFF Optimization Set 1](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2019-05-16-Roche-Optimization_Set)
- [SMIRNOFF Coverage Set 1](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2019-06-25-smirnoff99Frost-coverage)
- [OpenFF VEHICLe Set 1](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2019-07-02%20VEHICLe%20optimization%20dataset)
- [OpenFF Discrepancy Benchmark 1](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2019-07-05%20eMolecules%20force%20field%20discrepancies%201)
- [OpenFF Ehrman Informative Optimization v0.2](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2019-09-06-OpenFF-Informative-Set)
- [Pfizer discrepancy optimization dataset 1](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2019-09-07-Pfizer-discrepancy-optimization-dataset-1)
- [FDA optimization dataset 1](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2019-09-08-fda-optimization-dataset-1)
- [Kinase Inhibitors: WBO Distributions](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2019-11-27-kinase-inhibitor-optimization)
- [OpenFF Gen 2 Opt Set 1 Roche](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2020-03-20-OpenFF-Gen-2-Optimization-Set-1-Roche)
- [OpenFF Gen 2 Opt Set 2 Coverage](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2020-03-20-OpenFF-Gen-2-Optimization-Set-2-Coverage)
- [OpenFF Gen 2 Opt Set 3 Pfizer Discrepancy](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2020-03-20-OpenFF-Gen-2-Optimization-Set-3-Pfizer-Discrepancy)
- [OpenFF Gen 2 Opt Set 4 eMolecules Discrepancy](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2020-03-20-OpenFF-Gen-2-Optimization-Set-4-eMolecules-Discrepancy)
- [OpenFF Gen 2 Opt Set 5 Bayer](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2020-03-20-OpenFF-Gen-2-Optimization-Set-5-Bayer)
- [OpenFF Sandbox CHO PhAlkEthOH v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2020-09-18-OpenFF-Sandbox-CHO-PhAlkEthOH)
- [OpenFF Aniline Para Opt v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2021-04-02-OpenFF-Aniline-Para-Opt-v1.0)
- [OpenFF Industry Benchmark Season 1 v1.2](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2025-06-24-OpenFF-Industry-Benchmark-Season-1-v1.2)
- [OpenFF Gen2 Optimization Dataset Protomers v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2021-12-21-OpenFF-Gen2-Optimization-Set-Protomers)
- [OpenFF Protein Capped 1-mers 3-mers Optimization Dataset v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2022-05-30-OpenFF-Protein-Capped-1-mers-3-mers-Optimization)
- [OpenFF Iodine Chemistry Optimization Dataset v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2022-07-27-OpenFF-iodine-optimization-set)
- [XtalPi Shared Fragments OptimizationDataset v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2024-01-30-xtalpi-shared-fragments-optimization-v1.0)
- [XtalPi 20-percent Fragments OptimizationDataset v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2024-04-02-xtalpi-20-percent-fragments-optimization-v1.0)
- OpenFF Torsion Benchmark Supplement v1.0
- [OpenFF Torsion Multiplicity Optimization Training Coverage Supplement v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2024-06-20-OpenFF-Torsion-Multiplicity-Optimization-Training-Coverage-Supplement-v1.0)
- [OpenFF Torsion Multiplicity Optimization Benchmarking Coverage Supplement v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2024-06-24-OpenFF-Torsion-Multiplicity-Optimization-Benchmarking-Coverage-Supplement-v1.0)
- [OpenFF Iodine Fragment Opt v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2024-09-10-OpenFF-Iodine-Fragment-Opt-v1.0)
- [OpenFF Sulfur Optimization Training Coverage Supplement v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2024-09-11-OpenFF-Sulfur-Optimization-Training-Coverage-Supplement-v1.0)
- [OpenFF Sulfur Optimization Benchmarking Coverage Supplement v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2024-09-18-OpenFF-Sulfur-Hessian-Training-Coverage-Supplement-v1.0)
- [OpenFF Lipid Optimization Training Supplement v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2024-10-08-OpenFF-Lipid-Optimization-Training-Supplement-v1.0)
- [OpenFF Lipid Optimization Benchmark Supplement v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2024-10-30-OpenFF-Lipid-Optimization-Benchmark-Supplement-v1.0)
- [OpenFF Cresset Additional Coverage Optimizations v4.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2025-03-06-OpenFF-Cresset-Additional-Coverage-Optimizations-v4.0)
- [OpenFF Protein PDB 4-mers v4.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2025-03-05-OpenFF-Protein-PDB-4mer-v4.0)
- [OpenFF Additional Generated ChEMBL Optimizations v4.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2025-04-14-OpenFF-Additional-ChEMBL-Fragment-Optimizations-v4.0)
  
### Torsiondrive Datasets

The optimization records were selected to maximize chemical diversity using a selection of record IDs listed in the [Sage 2.3.0 repository](https://raw.githubusercontent.com/openforcefield/ash-sage-rc2/32345dddeb6cb249367059fd99607ac2950a5c86/03_fit-valence/02_curate-data/output/torsiondrives-v3.json). These records came from the following datasets:

- [OpenFF Group1 Torsions](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2019-05-01-OpenFF-Group1-Torsions)
- [OpenFF Group1 Torsions 2](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2020-01-31-OpenFF-Group1-Torsions-2)
- [OpenFF Group1 Torsions 3](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2020-02-10-OpenFF-Group1-Torsions-3)
- [SMIRNOFF Coverage Torsion Set 1](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2019-07-01-smirnoff99Frost-coverage-torsion)
- [OpenFF Substituted Phenyl Set 1](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2019-07-25-phenyl-set)
- [Pfizer discrepancy torsion dataset 1](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2019-09-07-Pfizer-discrepancy-torsion-dataset-1)
- [OpenFF Primary Benchmark 1 Torsion Set](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2019-12-05-OpenFF-Benchmark-Primary-1-torsion)
- [OpenFF Gen 2 Torsion Set 1 Roche 2](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2020-03-23-OpenFF-Gen-2-Torsion-Set-1-Roche-2)
- [OpenFF Gen 2 Torsion Set 2 Coverage 2](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2020-03-23-OpenFF-Gen-2-Torsion-Set-2-Coverage-2)
- [OpenFF Gen 2 Torsion Set 3 Pfizer Discrepancy 2](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2020-03-23-OpenFF-Gen-2-Torsion-Set-3-Pfizer-Discrepancy-2)
- [OpenFF Gen 2 Torsion Set 4 eMolecules Discrepancy 2](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2020-03-23-OpenFF-Gen-2-Torsion-Set-4-eMolecules-Discrepancy-2)
- [OpenFF Gen 2 Torsion Set 5 Bayer 2](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2020-03-26-OpenFF-Gen-2-Torsion-Set-5-Bayer-2)
- [OpenFF Gen 2 Torsion Set 6 supplemental 2](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2020-03-26-OpenFF-Gen-2-Torsion-Set-6-supplemental-2)
- [OpenFF Fragmenter Validation 1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2020-04-28-Fragmenter-test)
- [OpenFF DANCE 1 eMolecules t142 v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2020-06-01-DANCE-1-eMolecules-t142-selected)
- [OpenFF Rowley Biaryl v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2020-06-17-OpenFF-Biaryl-set)
- [OpenFF-benchmark-ligand-fragments-v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2020-07-27-OpenFF-Benchmark-Ligands)
- [OpenFF Protein Fragments TorsionDrives v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2020-09-16-OpenFF-Protein-Fragments-TorsionDrives)
- [OpenFF WBO Conjugated Series v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2021-01-25-OpenFF-Conjugated-Series)
- [OpenFF Amide Torsion Set v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2021-03-23-OpenFF-Amide-Torsion-Set-v1.0)
- [OpenFF-benchmark-ligand-fragments-v2.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2021-08-10-OpenFF-JACS-Fragments-v2.0)
- [OpenFF multiplicity correction torsion drive data v1.1](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2022-04-29-OpenFF-multiplicity-correction-torsion-drive-data)
- [OpenFF Protein Capped 3-mer Omega v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2023-02-06-OpenFF-Protein-Capped-3-mer-Omega)
- [XtalPi Shared Fragments TorsiondriveDataset v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2024-01-30-xtalpi-shared-fragments-torsiondrive-v1.0)
- [OpenFF Torsion Coverage Supplement v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2024-02-29-OpenFF-Torsion-Coverage-Supplement-v1.0)
- [OpenFF RNA Dinucleoside Monophosphate TorsionDrives v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2024-03-26-OpenFF-RNA-Dinucleoside-Monophosphate-TorsionDrives-v1.0)
- [XtalPi 20-percent Fragments TorsiondriveDataset v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2024-04-02-xtalpi-20-percent-fragments-torsiondrive-v1.0)
- [OpenFF Torsion Drive Supplement v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2024-04-17-OpenFF-Torsion-Drive-Supplement-v1.0)
- [OpenFF Torsion Multiplicity Torsion Drive Coverage Supplement v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2024-06-14-OpenFF-Torsion-Multiplicity-Torsion-Drive-Coverage-Supplement-v1.0)
- [OpenFF Phosphate Torsion Drives v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2024-07-17-OpenFF-Phosphate-Torsion-Drives-v1.0)
- [OpenFF Alkane Torsion Drives v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2024-08-09-OpenFF-Alkane-Torsion-Drives-v1.0)
- [OpenFF Cresset Additional Coverage TorsionDrives v4.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2025-02-12-OpenFF-Cresset-Additional-Coverage-TorsionDrives-v4.0)
- [OpenFF Additional Generated ChEMBL TorsionDrives 4.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2025-04-10-OpenFF-Additional-ChEMBL-Fragment-TorsionDrives-v4.0)
- [OpenFF Additional Generated Guanidine Amidine Derivative and O-Linker TorsionDrives 4.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2025-04-10-OpenFF-Additional-Generated-Guanidine-Amidine-Derivative-O-linkers-TorsionDrives-4.0)
- [OpenFF Gen3 Torsion Set v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2021-04-09-OpenFF-Gen3-Torsion-Set-v1.0)

## General Information

* Date: 2026-01-27
* Name: OpenFF SMIRNOFF Sage 2.3.0
* Dataset submitter: Jennifer A Clark
* Dataset curator: Lily Wang
---
* Class: OpenFF Optimization Dataset
* Dataset Type: optimization
* Purpose: B3LYP-D3BJ/DZVP conformers for training OpenFF 2.3.0 Sage with AshGC v1.0 NAGL partial charge model.
* Number of unique molecules: 4696
* Number of conformers: 4696
* Number of conformers (min, mean, max): 1.00, 1.00, 1.00
* Molecular weight (min, mean, max): 32.05, 207.67, 878.25
* Charges: -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0
---
* Class: OpenFF TorsionDrive Dataset
* Dataset Type: torsiondrive
* Purpose: B3LYP-D3BJ/DZVP conformers for training torsion drives for OpenFF 2.3.0 Sage with AshGC v1.0 NAGL partial charge model.
* Number of unique molecules: 1265
* Number of driven torsions: 1371
* Number of conformers: 1265
* Number of conformers (min, mean, max): 1, 1, 1
* Molecular weight (min, mean, max): 32.05, 165.33, 511.27
* Charges: -3.0, -2.0, -1.0, 0.0, 1.0, 2.0

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

Note that the specifications used to produce the records in this dataset had slight variations, so although each entry corresponds to one specification and therefore one record, there are multiple specifications. Because the original datasets all named their respective specifications "default", here we append that name with the dataset ids of the datasets that share that specification.

* Elements: {N, P, Cl, O, Br, S, I, H, F, C}
---
* Program: geometric
* Program Keywords: {'qccnv': True, 'reset': True, 'enforce': 0.1, 'epsilon': 0, 'coordsys': 'tric'}
* QC Specifications: default-232
  * program: psi4
  * driver: SinglepointDriver.deferred
  * method: b3lyp-d3(bj)
  * basis: dzvp
  * keywords: {'maxiter': 200, 'scf_properties': ['dipole', 'quadrupole', 'wiberg_lowdin_indices', 'mayer_indices']}
  * protocols: {'wavefunction': <WavefunctionProtocolEnum.none: 'none'>, 'stdout': True, 'error_correction': {'default_policy': True, 'policies': None}, 'native_files': <NativeFilesProtocolEnum.none: 'none'>}
  * SCF Properties:
    * dipole
    * quadrupole
    * wiberg_lowdin_indices
    * mayer_indices
---
* Program: geometric
* Program Keywords: {}
* QC Specifications: default-281
  * program: psi4
  * driver: SinglepointDriver.deferred
  * method: b3lyp-d3bj
  * basis: dzvp
  * keywords: {'maxiter': 200, 'scf_properties': ['dipole', 'quadrupole', 'wiberg_lowdin_indices']}
  * protocols: {'wavefunction': <WavefunctionProtocolEnum.none: 'none'>, 'stdout': True, 'error_correction': {'default_policy': True, 'policies': None}, 'native_files': <NativeFilesProtocolEnum.none: 'none'>}
  * SCF Properties:
    * dipole
    * quadrupole
    * wiberg_lowdin_indices
---
* Program: geometric
* Program Keywords: {'tmax': 0.3, 'check': 0, 'qccnv': False, 'reset': False, 'trust': 0.1, 'molcnv': False, 'enforce': 0, 'epsilon': 1e-05, 'maxiter': 300, 'coordsys': 'tric', 'convergence_set': 'gau'}
* QC Specifications: default-296
  * program: psi4
  * driver: SinglepointDriver.deferred
  * method: b3lyp-d3bj
  * basis: dzvp
  * keywords: {'maxiter': 200, 'scf_properties': ['dipole', 'quadrupole', 'wiberg_lowdin_indices', 'mayer_indices']}
  * protocols: {'wavefunction': <WavefunctionProtocolEnum.none: 'none'>, 'stdout': True, 'error_correction': {'default_policy': True, 'policies': None}, 'native_files': <NativeFilesProtocolEnum.none: 'none'>}
  * SCF Properties:
    * dipole
    * quadrupole
    * wiberg_lowdin_indices
    * mayer_indices
---
* Program: geometric
* Program Keywords: {'tmax': 0.3, 'check': 0, 'qccnv': False, 'reset': True, 'trust': 0.1, 'molcnv': False, 'enforce': 0, 'epsilon': 1e-05, 'maxiter': 300, 'coordsys': 'dlc', 'convergence_set': 'gau'}
* QC Specifications: default-315-453-372
  * program: psi4
  * driver: SinglepointDriver.deferred
  * method: b3lyp-d3bj
  * basis: dzvp
  * keywords: {'maxiter': 200, 'scf_properties': ['dipole', 'quadrupole', 'wiberg_lowdin_indices', 'mayer_indices']}
  * protocols: {'wavefunction': <WavefunctionProtocolEnum.none: 'none'>, 'stdout': True, 'error_correction': {'default_policy': True, 'policies': None}, 'native_files': <NativeFilesProtocolEnum.none: 'none'>}
  * SCF Properties:
    * dipole
    * quadrupole
    * wiberg_lowdin_indices
    * mayer_indices
---
* Program: geometric
* Program Keywords: {'tmax': 0.3, 'check': 0, 'qccnv': False, 'reset': True, 'trust': 0.1, 'molcnv': False, 'enforce': 0, 'epsilon': 1e-05, 'maxiter': 300, 'coordsys': 'dlc', 'convergence_set': 'gau'}
* QC Specifications: default-345-365
  * program: psi4
  * driver: SinglepointDriver.deferred
  * method: b3lyp-d3bj
  * basis: dzvp
  * keywords: {'maxiter': 200, 'scf_properties': ['dipole', 'quadrupole', 'wiberg_lowdin_indices', 'mayer_indices', 'mbis_charges']}
  * protocols: {'wavefunction': <WavefunctionProtocolEnum.none: 'none'>, 'stdout': True, 'error_correction': {'default_policy': True, 'policies': None}, 'native_files': <NativeFilesProtocolEnum.none: 'none'>}
  * SCF Properties:
    * dipole
    * quadrupole
    * wiberg_lowdin_indices
    * mayer_indices
    * mbis_charges
---
* Program: geometric
* Program Keywords: {'tmax': 0.3, 'check': 0, 'qccnv': False, 'reset': True, 'trust': 0.1, 'molcnv': False, 'enforce': 0.0, 'epsilon': 1e-05, 'maxiter': 300, 'coordsys': 'dlc', 'convergence_set': 'GAU'}
* QC Specifications: default-379-383-385-387-388-392-393-396-399-412-415-416-426
  * program: psi4
  * driver: SinglepointDriver.deferred
  * method: b3lyp-d3bj
  * basis: dzvp
  * keywords: {'maxiter': 200, 'scf_properties': ['dipole', 'quadrupole', 'wiberg_lowdin_indices', 'mayer_indices']}
  * protocols: {'wavefunction': <WavefunctionProtocolEnum.none: 'none'>, 'stdout': True, 'error_correction': {'default_policy': True, 'policies': None}, 'native_files': <NativeFilesProtocolEnum.none: 'none'>}
  * SCF Properties:
    * dipole
    * quadrupole
    * wiberg_lowdin_indices
    * mayer_indices
---
* Program: geometric
* Program Keywords: {}
* QC Specifications: default-41
  * program: psi4
  * driver: SinglepointDriver.deferred
  * method: b3lyp-d3(bj)
  * basis: dzvp
  * keywords: {}
  * protocols: {'wavefunction': <WavefunctionProtocolEnum.none: 'none'>, 'stdout': True, 'error_correction': {'default_policy': True, 'policies': None}, 'native_files': <NativeFilesProtocolEnum.none: 'none'>}
---
* Program: geometric
* Program Keywords: {}
* QC Specifications: default-43-45-50
  * program: psi4
  * driver: SinglepointDriver.deferred
  * method: b3lyp-d3bj
  * basis: dzvp
  * keywords: {}
  * protocols: {'wavefunction': <WavefunctionProtocolEnum.none: 'none'>, 'stdout': True, 'error_correction': {'default_policy': True, 'policies': None}, 'native_files': <NativeFilesProtocolEnum.none: 'none'>}
---
* Program: geometric
* Program Keywords: {}
* QC Specifications: default-68-69-251-253-255-254-270
  * program: psi4
  * driver: SinglepointDriver.deferred
  * method: b3lyp-d3bj
  * basis: dzvp
  * keywords: {'maxiter': 200, 'scf_properties': ['dipole', 'quadrupole', 'wiberg_lowdin_indices', 'mayer_indices']}
  * protocols: {'wavefunction': <WavefunctionProtocolEnum.none: 'none'>, 'stdout': True, 'error_correction': {'default_policy': True, 'policies': None}, 'native_files': <NativeFilesProtocolEnum.none: 'none'>}
  * SCF Properties:
    * dipole
    * quadrupole
    * wiberg_lowdin_indices
    * mayer_indices

## Torsiondrive Metadata

Note that the specifications used to produce the records in this dataset had slight variations, so although each entry corresponds to one specification and therefore one record, there are multiple specifications. Because the original datasets all named their respective specifications "default", here we append that name with the dataset ids of the datasets that share that specification.

* Elements: {Br, C, Cl, F, H, I, N, O, P, S}
---
* Program: torsiondrive
* Optimization Program: geometric
* Optimization Keywords: {'qccnv': True, 'reset': True, 'enforce': 0.1, 'epsilon': 0, 'coordsys': 'tric'}
* QC Specifications: default-57
  * program: psi4
  * driver: SinglepointDriver.deferred
  * method: b3lyp-d3(bj)
  * basis: dzvp
  * keywords: {'maxiter': 200, 'scf_properties': ['dipole', 'quadrupole', 'wiberg_lowdin_indices', 'mayer_indices']}
  * protocols: {'wavefunction': <WavefunctionProtocolEnum.none: 'none'>, 'stdout': True, 'error_correction': {'default_policy': True, 'policies': None}, 'native_files': <NativeFilesProtocolEnum.none: 'none'>}
  * SCF Properties:
    * dipole
    * quadrupole
    * wiberg_lowdin_indices
    * mayer_indices
---
* Program: torsiondrive
* Optimization Program: geometric
* Optimization Keywords: {'qccnv': True, 'reset': True, 'enforce': 0.1, 'epsilon': 0, 'coordsys': 'tric'}
* QC Specifications: default-242-243-48-70-235-256-257-258-259-265-266-278-282
  * program: psi4
  * driver: SinglepointDriver.deferred
  * method: b3lyp-d3bj
  * basis: dzvp
  * keywords: {'maxiter': 200, 'scf_properties': ['dipole', 'quadrupole', 'wiberg_lowdin_indices', 'mayer_indices']}
  * protocols: {'wavefunction': <WavefunctionProtocolEnum.none: 'none'>, 'stdout': True, 'error_correction': {'default_policy': True, 'policies': None}, 'native_files': <NativeFilesProtocolEnum.none: 'none'>}
  * SCF Properties:
    * dipole
    * quadrupole
    * wiberg_lowdin_indices
    * mayer_indices
---
* Program: torsiondrive
* Optimization Program: geometric
* Optimization Keywords: {'tmax': 0.3, 'check': 0, 'qccnv': True, 'reset': True, 'trust': 0.1, 'molcnv': False, 'enforce': 0.1, 'epsilon': 0, 'maxiter': 300, 'coordsys': 'tric', 'convergence_set': 'gau'}
* QC Specifications: default-283-289-291
  * program: psi4
  * driver: SinglepointDriver.deferred
  * method: b3lyp-d3bj
  * basis: dzvp
  * keywords: {'maxiter': 200, 'scf_properties': ['dipole', 'quadrupole', 'wiberg_lowdin_indices', 'mayer_indices']}
  * protocols: {'wavefunction': <WavefunctionProtocolEnum.none: 'none'>, 'stdout': True, 'error_correction': {'default_policy': True, 'policies': None}, 'native_files': <NativeFilesProtocolEnum.none: 'none'>}
  * SCF Properties:
    * dipole
    * quadrupole
    * wiberg_lowdin_indices
    * mayer_indices
---
* Program: torsiondrive
* Optimization Program: geometric
* Optimization Keywords: {'tmax': 0.3, 'check': 0, 'qccnv': True, 'reset': True, 'trust': 0.1, 'molcnv': False, 'enforce': 0.1, 'epsilon': 0, 'maxiter': 300, 'coordsys': 'dlc', 'convergence_set': 'gau'}
* QC Specifications: default-308-314-324-370-374-317
  * program: psi4
  * driver: SinglepointDriver.deferred
  * method: b3lyp-d3bj
  * basis: dzvp
  * keywords: {'maxiter': 200, 'scf_properties': ['dipole', 'quadrupole', 'wiberg_lowdin_indices', 'mayer_indices']}
  * protocols: {'wavefunction': <WavefunctionProtocolEnum.none: 'none'>, 'stdout': True, 'error_correction': {'default_policy': True, 'policies': None}, 'native_files': <NativeFilesProtocolEnum.none: 'none'>}
  * SCF Properties:
    * dipole
    * quadrupole
    * wiberg_lowdin_indices
    * mayer_indices
---
* Program: torsiondrive
* Optimization Program: geometric
* Optimization Keywords: {'qccnv': True, 'reset': True, 'enforce': 0.1, 'epsilon': 0, 'coordsys': 'tric'}
* QC Specifications: default-36
  * program: psi4
  * driver: SinglepointDriver.deferred
  * method: b3lyp-d3(bj)
  * basis: dzvp
  * keywords: {}
  * protocols: {'wavefunction': <WavefunctionProtocolEnum.none: 'none'>, 'stdout': True, 'error_correction': {'default_policy': True, 'policies': None}, 'native_files': <NativeFilesProtocolEnum.none: 'none'>}
---
* Program: torsiondrive
* Optimization Program: geometric
* Optimization Keywords: {}
* QC Specifications: default-378-380-381-382-384-386-389-390-413-424-423
  * program: psi4
  * driver: SinglepointDriver.deferred
  * method: b3lyp-d3bj
  * basis: dzvp
  * keywords: {'maxiter': 200, 'scf_properties': ['dipole', 'quadrupole', 'wiberg_lowdin_indices', 'mayer_indices']}
  * protocols: {'wavefunction': <WavefunctionProtocolEnum.none: 'none'>, 'stdout': True, 'error_correction': {'default_policy': True, 'policies': None}, 'native_files': <NativeFilesProtocolEnum.none: 'none'>}
  * SCF Properties:
    * dipole
    * quadrupole
    * wiberg_lowdin_indices
    * mayer_indices