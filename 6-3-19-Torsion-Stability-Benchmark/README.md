# Fragment Stability Benchmark

### General Information
 - Date: 6/9/2019
 - Class: Fragmentation Research
 - Purpose: Further examination of fragmentation schemes to observe x.
 - Collection: TorsionDriveDataset
 - Name: Fragment Stability Benchmark
 - Number of Entries: 4 (expandable)
 - Submitter: Chaya Stern

### Generation pipeline
 - To find the torsion to drive, the bond in the parent molecule had to be mapped to the fragment. The function `get_torsion` in the `generate.py` script does that.
 - `chemi.generate_conformers` was used to generate max 10 conformers. 
 - The original SMILES were chosen based on which rotatable bond in the filtered kinase inhibitor set were the most senstivive to remote substituents (WBO changed by more than 0.1)
 Scripts for choosing the original SMILES live here (https://github.com/choderalab/fragmenter_data/tree/master/frag_opt/compare_growth_path)

### Quantum Chemistry Information
 - Theory: OpenFF high-throughput standard QC reference
 - Additional Properties: Computes and saves Wiberg bond-orders
 
### Manifest
* `generate.py` - script to generate torsiondrive inputs
* `idelalisib_frags.json` - input SMILES with some provenance
* `requirements.txt` - conda list
* `stability_benchmark_inputs.tar.gz` - torsiondrive jobs inputs
