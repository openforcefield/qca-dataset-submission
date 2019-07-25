# Biphenyl set

### General Information
 - Date: 7/25/2019
 - Class: Fragmentation Research
 - Purpose: Initial motivational example
 - Collection: TorsionDriveDataset
 - Name: Biphenyl set
 - Number of Entries: 10 (expandable)
 - Submitter: Chaya Stern

### Generation pipeline
 - The original molecules were provided by Chris Bayely
 - To generate input, run `python generate.py`
 
### Quantum Chemistry Information
 - Theory: OpenFF high-throughput standard QC reference
 - Additional Properties: Computes and saves Wiberg bond-orders
 
### Manifest
* `generate.py` - script to generate torsiondrive inputs
* `biphenyls.smi` - input SMILES 
* `biphenyls.pdf` - visualization of input
* `requirements.txt` - conda list
* `biphenyl_set_inputs.tar.gz` - torsiondrive jobs inputs
