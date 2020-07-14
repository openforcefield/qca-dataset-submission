# OpenFF Group1 Torsions 3

### General Information
 - Date: 2/10/2020
 - Class: Force Field Parameterization
 - Purpose: Generation of additional data for fitting of `t128` and `t129`
 - Collection: TorsiondriveDataset
 - Name: OpenFF Group1 Torsions 3
 - Number of Entries: 6 1-D torsions 
 - Submitter: Hyesu Jang

### Generation pipeline
1. The first step of generating conformers is copyed from "2019-07-01-smirnoff99Frost-coverage-torsion"
2. Group all conformers by their index stored in mdata['cmiles_identifiers']['canonical_isomeric_smiles']
3. For each molecule, analyze the SMIRKs of each torsion, find torsions whose SMIRKs matches to new torsion terms, skip in-ring rotations, check if the SMIRKs already reached a target count (5).
4. If not enough SMIRKs are in the selected_torsions yet, add this torsion to selected_torsions, save as JSON file.
5. Create dataset using the JSON file.

### Notes
 - The details of the torsion selection can be found in `select.log`
 - The SMIRKs target count of 5 is chosen to generate a reasonable number of torsions for scanning.

### Manifest
 - `chosen_supplemented.smi` - input smi file 
 - `01_generate.py` - script to generate OptimizationDataset inputs
 - `optimization_inputs.json` - input molecules
 - `cmiles_failures.json` - Molecules that rdkit failed to generate standardized tautomer
 - `omega_failures.json` - Molecules that Omega failed to generate conformers
 - `skipped_ions.json` - Molecules that were skipped (ions)
 - `requirements.txt` - versions of toolkits used for the python scripts
 - `02_select_torsions.py` - Script to select torsions from generated conformers
 - `select.log` - Stdout of the 02_select_torsions.py script for record
 - `03_create_torsiondrive_dataset.py` - Script for creating the TorsiondriveDataset
 - `04_visualize.py` - create a PDF of all species in the dataset, highlighting central torsion and atoms involved in the torsion (copied from "2019-09-07-Pfizer-discrepancy-torsion-dataset-1")
 - `param_valence.offxml` - force field file 


### Usage
1. Generate the conformers and save in JSON format
    ```
    python 01_generate.py
    ```
    A new file called "optimization_inputs.json" will be generated.

2. Group all conformers, and select torsions based on SMIRNOFF coverage.
    ```
    python 02_select_torsions.py > select.log
    ```
    A new file called "selected_torsions.json" will be generated. Details of the run can be found in "select.log"

3. Create dataset on QCFractal server
    ```
    python 03_create_torsiondrive_dataset.py
    ```
