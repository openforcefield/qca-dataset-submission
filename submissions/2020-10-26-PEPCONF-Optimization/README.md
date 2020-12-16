# PEPCONF OptimizationDataset

### Description

OptimizationDataset of short peptides in various contexts, including disulfide bridges.
The source PEPCONF dataset is documented in [Nature Scientific Data](https://www.nature.com/articles/sdata2018310), and available on [GitHub](https://github.com/aoterodelaroza/pepconf).
This dataset extracts the molecules that were simulated, but uses the QCSubmit infrastructure to generate a new `OptimizationDataset`, so does not use the original conformers.

### General Information

- Date: 2020.10.26
- Class: OpenFF OptimizationDataset
- Purpose: OptimizationDataset of short peptides in various contexts, including disulfide bridges.
- Collection: OptimizationDataset
- Name: OpenFF PEPCONF OptimizationDataset
- Number of unique molecules: 736
- Number of conformers: 7560
- Submitter: John D. Chodera

### QCSubmit generation pipeline

Extracting SMILES:

```bash
git clone https://github.com/aoterodelaroza/pepconf
python extract-pepconf-smiles.py
```

Generating stereoisomers, conformers, and QCSubmit dataset:

- `Dataset Preparation.ipynb`: This notebook shows how the Optimization dataset was prepared from the input files.

### QCSubmit Manifest

- `extract-pepconf-smiles.py` - use the OpenEye Toolkit to perceive unique chemical structures from the PEPCONF repo
- `pepconf.csv` - the unique SMILES with associated names from PEPCONF (absent `_#` conformer suffixes)
- `Dataset Preparation.ipynb`: Dataset creation notebook with details on decisions made for submission.
- `molecules.pdf`: A pdf file containing molecule 2D structures.
- `molecules.smi`: SMILES for every molecule in the submission.

### Metadata

- elements {'C', 'H', 'N', 'O', 'S'}
- unique molecules: 736
- optimizations: 7560
- scf properties:
    - `dipole`
    - `quadrupole`
    - `wiberg_lowdin_indices`
    - `mayer_indices`
- qc spec
    - default
        - method: B3LYP-D3BJ
        - basis: DZVP
        - program: psi4
    - ani2x-qm
        - method: WB97X
        - basis: 6-31G*
        - program: psi4
