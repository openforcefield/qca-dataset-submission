### Description
Torsion Drives to explore wavefront propagation for the TorsionDrive paper

### General Information
 - Date: 11/7/2019
 - Class: TorsionDrive Paper
 - Purpose: Explore wavefront propagation
 - Collection: TorsionDriveDataset
 - Name: TorsionDrive paper
 - Number of Entries: 2 1D torsion drives, 2 2D torsiondrives (expandable)
 - Submitter: Yudong Qiu

### Generation pipeline
1. `01_generate.py
    This generates the torsiondrive JSON inputs (`torsiondrive_inputs.json`)

2. `02_create_torsiondrive_dataset.py`
    This script reads `torsiondrive_inputs.json` and create the dataset.
