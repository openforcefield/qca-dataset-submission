# OpenFF NAGL2 ESP Timing Benchmark v1.0

## Description
A dataset of 1000 molecules, sub-sampled from the [OpenFF multi-Br ESP Fragment Conformers v1.1](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2023-11-30-OpenFF-multi-Br-ESP-Fragment-Conformers-v1.1-single-point) and the [OpenFF ESP Fragment Conformers v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2022-01-16-OpenFF-ESP-Fragment-Conformers-v1.0) datasets.

Sub-sampling was done in `subsample_esp_ds.ipynb` by selecting 990 conformers at random from the ESP Fragment dataset and 11 from the multi-Br dataset.

The purpose of this dataset is to estimate timing and storage requirements for this methodology, which is a candidate for NAGL2 dataset generation.

## General information
* Date: 2024-09-06
* Class: OpenFF Basic Dataset
* Purpose: ESP/electric field generation
* Name: OpenFF NAGL2 ESP Timing Benchmark v1.0
* Number of unique molecules: 995
* Number of conformers: 1001
* Number of conformers (min, mean, max): 1 1.006 2
* Molecular weight (min, mean, max): 58.10, 149.84, 329.82
* Charges: -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0
* Dataset submitter: Alexandra McIsaac
* Dataset generator: Alexandra McIsaac

## QCSubmit generation pipeline
* `subsample_esp_ds.ipynb` was used to reduce the parent datasets into a small dataset to use for benchmarking purposes.
* `generate-dataset.ipynb` was used to transform the small sub-sampled dataset into a new dataset for submission.

## QCSubmit Manifest
* `dataset.json.bz2`: compressed dataset ready for submission
* `dataset.pdf`: Visualization of dataset molecules
* `dataset.smi`: Smiles strings for dataset molecules
* `subsample_esp_ds.ipynb`: Notebook showing how dataset was pruned from larger parent datasets.
* `generate-dataset.ipynb`: Notebook describing dataset generation and submission
* `input-environment.yaml`: Environment file used to create Python environment for the notebook
* `input-environment-full.yaml`: Fully-resolved environment used to execute the notebook.

## Metadata
* elements: {'P', 'S', 'N', 'C', 'Cl', 'F', 'Br', 'O', 'H'}
* unique molecules: 995
* Spec: wB97X-V/def2-TZVPPD
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
     * name: wB97X-V/def2-TZVPPD
     * method: wB97X-V
     * basis: def2-TZVPPD
     * keywords: {'dft_spherical_points': 590, 'dft_radial_points': 99}
