# OpenFF Industry Benchmark Season 1 v1.2

## Description

This dataset is the public counterpart of the OpenFF Industry Benchmark Season 1.
Each industry partner has selected a range of diverse molecules which represent their current chemical interests.
The dataset will be used in conjunction with private counterparts also designed by each partner to give an unbiased assessment of the progress and current performance of the OpenFF line of force fields in comparison with other contemporary force fields.

The v1.1 dataset features corrected Merck (MRK) molecules with explicit hydrogens.
The original v1.0 dataset did not have explicit hydrogens on these molecules, resulting in poor starting conformers that have largely failed to geometry optimize under QM.
The v1.1 dataset was prepared from the v1.0 dataset, excising the MRK molecules and replacing them with the explicit hydrogen variants prepared using the [Season 1 protocol](https://openforcefield.atlassian.net/wiki/spaces/PS/pages/971898891/Optimization+Benchmarking+Protocol+-+Season+1) via `openff-benchmark`.

This v1.2 dataset has removed a subset of structures from v1.1 that exhibited unrealistic angles, which proved to hinder progress in converging force field parameters. The outlier records are recorded in the [sage-2.2.0 repository](https://github.com/openforcefield/sage-2.2.0/blob/main/05_benchmark_forcefield/process_bm/problem_ids/all_r7_outliers.txt). It also applies our standard post-submission filters.

### General Information

- Date: 2025-06-24
- Name: OpenFF Industry Benchmark Season 1 v1.2
- Purpose: The combination of all publicly chosen compound sets by industry partners from the OpenFF Industry Benchmark Season 1 with standard post-submission filters applied and unrealistic conformers removed.
- Dataset Submitter: Jennifer A. Clark
- Dataset Curator: David Dotson and Joshua Horton
- Dataset Generator: David Dotson and Joshua Horton

- Class: OpenFF Optimization Dataset
- Dataset Type: optimization
- Number of unique molecules:   9835
- Number of filtered molecules: 0
- Number of conformers:         74585
- Number of conformers (min mean max): 1.00, 7.58, 10.00
- Molecular weight (min, mean, max): 16.04, 348.58, 1105.15
- Set of charges: -2.0, -1.0, 0.0, 1.0, 2.0

### QCSubmit generation pipeline

* `main.ipynb`: A jupyter notebook which shows how the optimization dataset was prepared from the input files.

## Manifest

* `main.ipynb`: A jupyter notebook which shows how the optimization dataset was prepared from the input files.
* `ds_info.json`: Metadata information for optimization dataset record imported by `main.ipynb`
* `environment.yaml`: A file to reproduce the conda env used to generate this dataset.
* `environment_full.yaml`: A file to reproduce the conda env used to generate this dataset.
* `scaffold.json.bz2`: The basic dataset information read by qcportal.external.scaffold.

## Optimization Metadata

* Elements: {Br, C, Cl, F, H, N, O, P, S}
* QC Specifications: default
  * method: B3LYP-D3BJ
  * basis: DZVP
  * store_wavefunction: None
  * implicit_solvent: None
  * maxiter: 200
  * keywords: {}
  * SCF Properties:
    * dipole
    * quadrupole
    * wiberg_lowdin_indices
    * mayer_indices
