# Fragmenter Testing

### General Information
 - Date: 2020-04-28
 - Class: Fragmentation Research
 - Purpose: Test fragmenter.
 - Collection: TorsionDriveDataset
 - Name: Fragmenter Testing
 - Number of Entries: 22 (expandable)
 - Submitter: Chaya Stern

### Generation pipeline
 - The molecules were taken from the [exhaustive fragmentation benchmark](https://github.com/choderalab/fragmenter_data/tree/master/combinatorial_fragmentation/verify_schemes).

### Quantum Chemistry Information
 - Theory: OpenFF high-throughput standard QC reference
 - Additional Properties: Computes and saves Wiberg bond-orders

### Manifest
* `01_generate.py` - script to generate torsiondrive inputs
* `02_create_torsiondrive_dataset.py` - script to generate torsiondive datase
* `03_visualize.py` - script to visualize input molecules
* `requirements.txt` - conda list
* `selected-torsions.tar.gz` - torsiondrive jobs inputs
* `torsiondrive-inputs.pdf` - visualize input molecules