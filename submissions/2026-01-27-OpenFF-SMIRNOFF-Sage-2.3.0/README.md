# OpenFF SMIRNOFF Sage 2.3.0

## Description

A quantum chemical (QC) optimization and torsiondrive datasets generated at the OpenFF default level of theory, B3LYP-D3BJ/DZVP, and curated to train parameters in [OpenFF 2.3.0 Sage](https://github.com/openforcefield/ash-sage-rc2) with NAGL partial charge model AshGC.

Datasets available on [Zenodo]()

### Optimization Datasets

The optimization records were selected to maximize chemical diversity using a selection of record IDs listed in the [Sage 2.3.0 repository](https://raw.githubusercontent.com/openforcefield/ash-sage-rc2/32345dddeb6cb249367059fd99607ac2950a5c86/03_fit-valence/02_curate-data/output/optimizations-single-v3.json). These records came from the following datasets:

- "OpenFF Optimization Set 1",
- "SMIRNOFF Coverage Set 1",
- "OpenFF VEHICLe Set 1",
- "OpenFF Discrepancy Benchmark 1",
- "OpenFF Ehrman Informative Optimization v0.2",
- "Pfizer discrepancy optimization dataset 1",
- "FDA optimization dataset 1",
- "Kinase Inhibitors: WBO Distributions",
- "OpenFF Gen 2 Opt Set 1 Roche",
- "OpenFF Gen 2 Opt Set 2 Coverage",
- "OpenFF Gen 2 Opt Set 3 Pfizer Discrepancy",
- "OpenFF Gen 2 Opt Set 4 eMolecules Discrepancy",
- "OpenFF Gen 2 Opt Set 5 Bayer",
- "OpenFF Sandbox CHO PhAlkEthOH v1.0",
- "OpenFF Aniline Para Opt v1.0",
- "OpenFF Industry Benchmark Season 1 v1.2",
- "OpenFF Gen2 Optimization Dataset Protomers v1.0",
- "OpenFF Protein Capped 1-mers 3-mers Optimization Dataset v1.0",
- "OpenFF Iodine Chemistry Optimization Dataset v1.0",
- "XtalPi Shared Fragments OptimizationDataset v1.0",
- "XtalPi 20-percent Fragments OptimizationDataset v1.0",
- "OpenFF Torsion Benchmark Supplement v1.0",
- "OpenFF Torsion Multiplicity Optimization Training Coverage Supplement v1.0",
- "OpenFF Torsion Multiplicity Optimization Benchmarking Coverage Supplement v1.0",
- "OpenFF Iodine Fragment Opt v1.0",
- "OpenFF Sulfur Optimization Training Coverage Supplement v1.0",
- "OpenFF Sulfur Optimization Benchmarking Coverage Supplement v1.0",
- "OpenFF Lipid Optimization Training Supplement v1.0",
- "OpenFF Lipid Optimization Benchmark Supplement v1.0",
- "OpenFF Cresset Additional Coverage Optimizations v4.0",
- "OpenFF Protein PDB 4-mers v4.0",
- "OpenFF Additional Generated ChEMBL Optimizations v4.0"
  
### Torsiondrive Datasets

The optimization records were selected to maximize chemical diversity using a selection of record IDs listed in the [Sage 2.3.0 repository](https://raw.githubusercontent.com/openforcefield/ash-sage-rc2/32345dddeb6cb249367059fd99607ac2950a5c86/03_fit-valence/02_curate-data/output/torsiondrives-v3.json). These records came from the following datasets:

- "OpenFF Group1 Torsions",
- "OpenFF Group1 Torsions 2",
- "OpenFF Group1 Torsions 3",
- "SMIRNOFF Coverage Torsion Set 1",
- "OpenFF Substituted Phenyl Set 1",
- "Pfizer discrepancy torsion dataset 1",
- "OpenFF Primary Benchmark 1 Torsion Set",
- "OpenFF Gen 2 Torsion Set 1 Roche 2",
- "OpenFF Gen 2 Torsion Set 2 Coverage 2",
- "OpenFF Gen 2 Torsion Set 3 Pfizer Discrepancy 2",
- "OpenFF Gen 2 Torsion Set 4 eMolecules Discrepancy 2",
- "OpenFF Gen 2 Torsion Set 5 Bayer 2",
- "OpenFF Gen 2 Torsion Set 6 supplemental 2",
- "OpenFF Fragmenter Validation 1.0",
- "OpenFF DANCE 1 eMolecules t142 v1.0",
- "OpenFF Rowley Biaryl v1.0",
- "OpenFF-benchmark-ligand-fragments-v1.0",
- "OpenFF Protein Fragments TorsionDrives v1.0",
- "OpenFF WBO Conjugated Series v1.0",
- "OpenFF Amide Torsion Set v1.0",
- "OpenFF-benchmark-ligand-fragments-v2.0",
- "OpenFF multiplicity correction torsion drive data v1.1",
- "OpenFF Protein Capped 3-mer Omega v1.0",
- "XtalPi Shared Fragments TorsiondriveDataset v1.0",
- "OpenFF Torsion Coverage Supplement v1.0",
- "OpenFF RNA Dinucleoside Monophosphate TorsionDrives v1.0",
- "XtalPi 20-percent Fragments TorsiondriveDataset v1.0",
- "OpenFF Torsion Drive Supplement v1.0",
- "OpenFF Torsion Multiplicity Torsion Drive Coverage Supplement v1.0",
- "OpenFF Phosphate Torsion Drives v1.0",
- "OpenFF Alkane Torsion Drives v1.0",
- "OpenFF Cresset Additional Coverage TorsionDrives v4.0",
- "OpenFF Additional Generated ChEMBL TorsionDrives 4.0",
- "OpenFF Additional Generated Guanidine Amidine Derivative and O-Linker TorsionDrives 4.0",
- "OpenFF Gen3 Torsion Set v1.0"

## General Information

* Date: 2026-01-27
* Name: OpenFF SMIRNOFF Sage 2.3.0
* Dataset submitter: Jennifer A Clark
* Dataset curator: Lily Wang
---
* Class: OpenFF Optimization Dataset
* Dataset Type: optimization
* Purpose: B3LYP-D3BJ/DZVP conformers for training OpenFF 2.3.0 Sage with AshGC NAGL partial charge model.
* Number of unique molecules: 4696
* Number of conformers: 4696
* Number of conformers (min, mean, max): 1.00, 1.00, 1.00
* Molecular weight (min, mean, max): 32.05, 207.67, 878.25
* Charges: -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0
---
* Class: OpenFF TorsionDrive Dataset
* Dataset Type: torsiondrive
* Purpose: B3LYP-D3BJ/DZVP conformers for training torsion drives for OpenFF 2.3.0 Sage with AshGC NAGL partial charge model.
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
* QC Specifications: default-41
  * program: psi4
  * driver: SinglepointDriver.deferred
  * method: b3lyp-d3(bj)
  * basis: dzvp
  * keywords: {}
  * protocols: {'wavefunction': <WavefunctionProtocolEnum.none: 'none'>, 'stdout': True, 'error_correction': {'default_policy': True, 'policies': None}, 'native_files': <NativeFilesProtocolEnum.none: 'none'>}
---
* Program: geometric
* QC Specifications: default-43-45-50
  * program: psi4
  * driver: SinglepointDriver.deferred
  * method: b3lyp-d3bj
  * basis: dzvp
  * keywords: {}
  * protocols: {'wavefunction': <WavefunctionProtocolEnum.none: 'none'>, 'stdout': True, 'error_correction': {'default_policy': True, 'policies': None}, 'native_files': <NativeFilesProtocolEnum.none: 'none'>}
---
* Program: geometric
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
* Optimization Specification: geometric
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
* Optimization Specification: geometric
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
* Optimization Specification: geometric
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

* Program: torsiondrive
* Optimization Specification: geometric
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
* Optimization Specification: geometric
* QC Specifications: default-36
  * program: psi4
  * driver: SinglepointDriver.deferred
  * method: b3lyp-d3(bj)
  * basis: dzvp
  * keywords: {}
  * protocols: {'wavefunction': <WavefunctionProtocolEnum.none: 'none'>, 'stdout': True, 'error_correction': {'default_policy': True, 'policies': None}, 'native_files': <NativeFilesProtocolEnum.none: 'none'>}
---
* Program: torsiondrive
* Optimization Specification: geometric
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