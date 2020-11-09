# Biphenyl TorsionDriveDataset

### Description

Resubmission of `OpenFF Substituted Phenyl Set 1`, in response to metadata not being [stored in the places we need](https://github.com/openforcefield/qca-dataset-submission/pull/140#issuecomment-705675738) to use the dataset for MM reliably.

Description from original dataset below:

> - The original molecules were provided by Chris Bayly
> - The scripts and data to generate `phenyl_set_inputs.json` is [here](https://github.com/choderalab/fragmenter_data/blob/master/phenyl_benchmark/generate_torsiondrive_inputs.py)
> - Theory: OpenFF high-throughput standard QC reference
> - Additional Properties: Computes and saves Wiberg bond-orders
 
### General Information

- Date: 2020.11.09
- Class: OpenFF TorsionDriveDataset
- Purpose: 
- Collection: TorsionDriveDataset
- Name: OpenFF Substituted Phenyl Set 1 v2.0
- Number of unique molecules: 155
- Number of torsiondrives: 156
- Submitter: David Dotson

### QCSubmit generation pipeline

Extracting original molecules, conformers, and CMILES information; generating QCSubmit dataset:

- `Submission Preparation.ipynb`: This notebook shows how the Optimization dataset was prepared from the input files.

### QCSubmit Manifest

- `biphenyls_set_input.json`: original v1 submission molecules
- `phenyl_set_torsiondrive_inputs.json`: original v1 submission molecules
- `Submission Preparation.ipynb`: Dataset creation notebook with details on decisions made for submission.
- `molecules.pdf`: A pdf file containing molecule 2D structures, with driven torsions highlighted.
- `molecules.smi`: SMILES for every molecule in the submission.
- `dataset.json`: QCSubmit dataset submitted to QCArchive

### Metadata

- elements: {'Br', 'C', 'Cl', 'F', 'H', 'I', 'N', 'O'}
- unique molecules: 155
- number of torsiondrives: 156
- scf properties:
    - `dipole`
    - `quadrupole`
    - `wiberg_lowdin_indices`
    - `mayer_indices`
- qc spec
    - name: default
    - method: B3LYP-D3BJ
    - basis: DZVP
