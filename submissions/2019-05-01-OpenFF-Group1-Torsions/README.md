### Description
This serves as a record of previously submitted TorsionDriveDataset "OpenFF Group1 Torsions"
The files are copied from https://github.com/lpwgroup/forcebalance-qcarchive/releases/tag/v0.0.2

### General Information
 - Date: 5/1/2019
 - Class: OpenFF SMIRNOFF Release 1 Fitting
 - Purpose: Generate data for OpenFF SMIRNOFF fitting
 - Collection: TorsionDriveDataset
 - Name: OpenFF Group1 Torsions
 - Number of Entries: 819 1-D torsion scans
 - Submitter: Yudong Qiu

### Generation pipeline
1. `01_process_molecules.py OpenFF_references.sdf`
    This generates a folder called "processed_molecules", with each molecule in a separate file, like `001_C13H12.sdf`. The file in mol2 and xyz formats are also generated in subfolders.

2. `02_submit_group1_1d.py -s td_scan_configure.yaml -j processed_molecules/*.sdf`
    This goes through all molecules and generate the json file for submission, named `submit_torsion_options.json`.
    It also generate a file `torsion_submit_checkpoint.json` that contains information about jobs submitted.

3. `03_create_torsiondrive_dataset.py`
    This script reads `submit_torsion_options.json` and create the dataset.

4. `04_visualize.py`
    Optionally, this script creates a PDF file with compact 2D structures of the molecules.

### Notes
1. The provided file `OpenFF_references.sdf` contains 468 molecules.

2. Since there are too many dihedrals from all molecules, we apply a few filters to keep only the most important dihedrals.

3. The filters applied in order are 'heavy_atoms', 'no_ring', 'unique_center_bond'.

4. Since this is a record for a previous submission, please do not use this one as a template for submitting new datasets.
