### Description
Geometry optimization of kinase inhibitor conformers to explore WBO conformation dependency

### General Information
 - Date: 11/27/2019
 - Class: Kinase inhibitors WBO distribution
 - Purpose: Explore WBO conformation dependency
 - Collection: OptimizationDataset
 - Name: Kinase Inhibitors: WBO Distributions
 - Number of Entries: 5609 optimizations (expandable)
 - Submitter: Chaya Stern

### Generation pipeline
1. `01_generate.py`
    This generates the optimization JSON inputs (`torsiondrive_inputs.json.gz`)
    to run:
    `python 01_generate.py -i kinase_inhibitors.smi`

2. `02_create_optimization_dataset.py`
    This script reads `torsiondrive_inputs.json.gz` and create the dataset.

### Manifest
* `kinase_inhibitors.smi` - input molecules
* `01_generate.py` - generate conformers for optimization dataset
* `02_create_optimization_dataset.py` - create optimization dataset
* `kinase_inhibitors.pdf` - visualization of kinase inhibitors in dataset
* `optimization_inputs.json.gz` - optimization inputs JSON
