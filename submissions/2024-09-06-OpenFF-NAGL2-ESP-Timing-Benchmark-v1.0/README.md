# OpenFF NAGL2 ESP Timing Benchmark v1.0

## Description
A dataset of 1005 molecules, sub-sampled from the [OpenFF multi-Br ESP Fragment Conformers v1.1](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2023-11-30-OpenFF-multi-Br-ESP-Fragment-Conformers-v1.1-single-point), the [OpenFF Iodine Fragment Opt v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2024-09-10-OpenFF-Iodine-Fragment-Opt-v1.0), and the [OpenFF ESP Fragment Conformers v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2022-01-16-OpenFF-ESP-Fragment-Conformers-v1.0) datasets. Calculated at the PBE0/def2-TZVPPD level of theory

Sub-sampling was done in `subsample_esp_ds.ipynb` by selecting 945 conformers at random from the ESP Fragment dataset, 53 from the I fragment dataset, and 11 from the multi-Br dataset.

The purpose of this dataset is to estimate timing and storage requirements for this methodology, which is a candidate for NAGL2 dataset generation.

## General information
* Date: 2024-09-06
* Class: OpenFF Basic Dataset
* Purpose: ESP/electric field generation
* Name: OpenFF NAGL2 ESP Timing Benchmark v1.0
* Number of unique molecules: 1005
* Number of conformers: 1009
* Number of conformers (min, mean, max): 1, 1.004, 2
* Molecular weight (min, mean, max): 58.10, 155.15, 329.82
* Charges: -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0
* Dataset submitter: Alexandra McIsaac
* Dataset generator: Alexandra McIsaac

## QCSubmit generation pipeline
* `subsample_esp_ds.ipynb` was used to reduce the parent datasets into a small dataset to use for benchmarking purposes.
* `br_subsample_filtered.json`, `esp_subsample_filtered.json`, `i_subsample_filtered.json` sub-sampled datasets pulled from the multi-Br, ESP50k, and I fragment datasets, produced by `subsample_esp_ds.ipynb`. Used as inputs for `generate-dataset.ipynb`
* `generate-dataset.ipynb` was used to transform the small sub-sampled dataset into a new dataset for submission.

## QCSubmit Manifest
* `dataset.json.bz2`: compressed dataset ready for submission
* `dataset.pdf`: Visualization of dataset molecules
* `dataset.smi`: Smiles strings for dataset molecules
* `subsample_esp_ds.ipynb`: Notebook showing how dataset was pruned from larger parent datasets.
* `generate-dataset.ipynb`: Notebook describing dataset generation and submission
* `input-environment.yaml`: Environment file used to create Python environment for the notebook
* `input-environment-full.yaml`: Fully-resolved environment used to execute the notebook.
* `br_subsample_filtered.json`: subset of multi-Br dataset selected for this benchmark
* `i_subsample_filtered.json`: subset of I fragment dataset selected for this benchmark
* `esp_subsample_filtered.json`: subset of ESP50k dataset selected for this benchmark

## Metadata
* elements: {'H', 'P', 'C', 'N', 'O', 'S', 'Cl', 'I', 'Br', 'F'}
* unique molecules: 1005
* Spec: PBE0/def2-TZVPPD
  * SCF properties:
    * Dipole
    * Quadrupole
    * LowdinCharges
    * MullikenCharges
    * MBISCharges
    * MayerIndices
    * WibergLowdinIndices
    * DipolePolarizabilities
   * QC Spec:
     * name: PBE0/def2-TZVPPD
     * method: PBE0
     * basis: def2-TZVPPD
     * keywords: {'dft_spherical_points': 590, 'dft_radial_points': 99}
