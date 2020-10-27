# PEPCONF OptimizationDataset

### Description

OptimizationDataset of short peptides in various contexts, including disulfide bridges
The source PEPCONF dataset is documented in [Nature Scientific Data](https://www.nature.com/articles/sdata2018310), and available on [GitHub](https://github.com/aoterodelaroza/pepconf).
This dataset extracts the molecules that were simulated, but uses the QCSubmit infrastructure to generate a new `OptimizationDataset`, so does not use the original conformers.

### General Information

 - Date: 2020.10.26
 - Class: OpenFF OptimizationDataset
 - Purpose: OptimizationDataset of short peptides in various contexts, including disulfide bridges
 - Collection: OptimizationDataset
 - Name: OpenFF PEPCONF OptimizationDataset
 - Number of unique molecules:
 - Number of torsiondrives:
 - Submitter: John D. Chodera

### QCSubmit generation pipeline

Extracting SMILES:
```bash
git clone https://github.com/aoterodelaroza/pepconf
python extract-pepconf-smiles.py
```

### QCSubmit Manifest
* `extract-pepconf-smiles.py` - use the OpenEye Toolkit to perceive unique chemical structures from the PEPCONF repo
* `pepconf.csv` - the unique SMILES with associated names from PEPCONF (absent `_#` conformer suffixes)

### Metadata
