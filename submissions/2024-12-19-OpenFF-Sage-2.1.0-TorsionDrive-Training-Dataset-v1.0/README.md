# OpenFF Sage 2.1.0 TorsionDrive Training Dataset v1.0

## Description

A quantum chemical (QC) dataset of single point energies and properties was curated to train sage 
[OpenFF 2.1.0 Sage](https://github.com/openforcefield/sage-2.1.0/). This QC dataset with the OpenFF
 default level of theory, B3LYP-D3BJ/DZVP, is used to benchmark Sage geometries and energetics. 
 This is the complete TorsionDrive dataset used for training OpenFF 2.0.0 Sage, consisting of the 
 following datasets: 
 - [OpenFF Gen 2 Torsion Set 1 Roche 2](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2020-03-23-OpenFF-Gen-2-Torsion-Set-1-Roche-2)
 - [OpenFF Gen 2 Torsion Set 2 Coverage 2](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2020-03-23-OpenFF-Gen-2-Torsion-Set-2-Coverage-2)
 - [OpenFF Gen 2 Torsion Set 3 Pfizer Discrepancy 2](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2020-03-23-OpenFF-Gen-2-Torsion-Set-3-Pfizer-Discrepancy-2)
 - [OpenFF Gen 2 Torsion Set 4 eMolecules Discrepancy 2](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2020-03-23-OpenFF-Gen-2-Torsion-Set-4-eMolecules-Discrepancy-2)
 - [OpenFF Gen 2 Torsion Set 5 Bayer 2](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2020-03-26-OpenFF-Gen-2-Torsion-Set-5-Bayer-2)
 - [OpenFF Gen 2 Torsion Set 6 supplemental 2](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2020-03-26-OpenFF-Gen-2-Torsion-Set-6-supplemental-2)
 - [SMIRNOFF Coverage Torsion Set 1](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2019-07-01-smirnoff99Frost-coverage-torsion)
 - [OpenFF Group1 Torsions](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2019-05-01-OpenFF-Group1-Torsions)
 - [OpenFF Group1 Torsions 2](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2020-01-31-OpenFF-Group1-Torsions-2)
 - [OpenFF Group1 Torsions 3](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2020-02-10-OpenFF-Group1-Torsions-3)
 - [Pfizer discrepancy torsion dataset 1](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2019-09-07-Pfizer-discrepancy-torsion-dataset-1)
 - [OpenFF Gen3 Torsion Set v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2021-04-09-OpenFF-Gen3-Torsion-Set-v1.0)
 - [OpenFF Amide Torsion Set v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2021-03-23-OpenFF-Amide-Torsion-Set-v1.0)
 - [OpenFF WBO Conjugated Series v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2021-01-25-OpenFF-Conjugated-Series)
 - [OpenFF DANCE 1 eMolecules t142 v1.0](https://github.com/openforcefield/qca-dataset-submission/tree/7f8ed2ab6c8acc4521c8ca45ff4f587b20f0bcda/submissions/2020-06-01-DANCE-1-eMolecules-t142-selected)
 
 These combined datasets were filtered with:
 
 - `ElementFilter(allowed_elements=['H', 'C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br'])`
 - `HydrogenBondFilter(method='baker-hubbard')`
 - `ChargeCheckFilter()`

## General Information

* Date: 2024-12-19
* Class: OpenFF TorsionDrive Dataset
* Purpose: Complete set of training TorsionDrive data for OpenFF 2.1.0 Sage
* Name: OpenFF Sage 2.1.0 TorsionDrive Training Dataset v1.0
* Number of unique molecules: 953
* Number of filtered molecules: 0
* Number of driven torsions: 1300
* Number of conformers: 974
* Number of conformers (min, mean, max): 1.00, 1.02, 3.00
* Molecular weight (min, mean, max): 32.04, 185.54, 503.42
* Charges: -1.0, 0.0, 1.0
* Dataset submitter: Jennifer A Clark
* Dataset curator: Pavan Behara
* Dataset generator: Simon Boothroyd, John Chodera, Trevor Gokey, Hyesu Jang, Yudong Qiu, Bryon Tjanaka

## QCSubmit Generation Pipeline

* `generate-combined-dataset.py`: A python script which shows how the dataset was prepared from the input files.


## QCSubmit Manifest

* `generate-combined-dataset.py`: A python script which shows how the dataset was prepared from the input files.
* `ds_info.json`: Metadata information for dataset record imported by `generate-combined-dataset.py`
* `output.txt`: Captured output from `generate-combined-dataset.py`
* `dataset.json.bz2`: The basic dataset ready for submission.
* `dataset.pdf`: A pdf file containing molecule 2D structures.
* `dataset.smi`: SMILES for every molecule in the submission.
* `conda_env.yaml`: A file to reproduce the conda env used to generate this dataset.


## Metadata

* Elements: {Br, C, Cl, F, H, N, O, P, S}
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