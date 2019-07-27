# Biphenyl set

### General Information
 - Date: 7/25/2019
 - Class: Fragmentation Research
 - Purpose: Investigate WBO and torsion barrier relationship
 - Collection: TorsionDriveDataset
 - Name: phenyl set
 - Number of Entries: 174 (expandable)
 - Submitter: Chaya Stern

### Generation pipeline
 - The original molecules were provided by Chris Bayly
 - To generate biphenyl input, run `python generate.py`
 - The scripts and data to generate `phenyl_set_inputs.json` is [here](https://github.com/choderalab/fragmenter_data/blob/master/phenyl_benchmark/generate_torsiondrive_inputs.py
 
### Quantum Chemistry Information
 - Theory: OpenFF high-throughput standard QC reference
 - Additional Properties: Computes and saves Wiberg bond-orders
 
### Manifest
* `generate.py` - script to generate torsiondrive inputs
* `biphenyls.smi` - input SMILES 
* `biphenyls.pdf` - visualization of input
* `requirements.txt` - conda list
* `biphenyl_set_inputs.tar.gz` - torsiondrive jobs inputs
* `phenyl_set_inputs.tar.gz` - torsiondrive jobs inputs for the biphenyl set
* `phenyl_set_inputs.tar.gz` - torsiondrive job inputs for phenyl set
* `overlapping_indices.smi` - torsiondrive indices that are already in QCArchive 
