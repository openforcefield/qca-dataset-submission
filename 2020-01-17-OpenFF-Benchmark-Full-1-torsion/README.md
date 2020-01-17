# OpenFF Benchmark Primary 1 torsion

### General Information
 - Date: 1/17/2020
 - Class: Forcefield Benchmark 
 - Purpose: Validation of optimized force field torsion parameters
 - Collection: TorsiondriveDataset
 - Name: OpenFF Benchmark Full 1 torsion
 - Number of Entries: 595 1-D torsions (max 5 per SMIRKs)
 - Submitter: Hyesu Jang

### Generation pipeline
1. The first step of generating conformers is copyed from "2019-07-01-smirnoff99Frost-coverage-torsion"
2. Group all conformers by their index stored in mdata['cmiles_identifiers']['canonical_isomeric_smiles']
3. For each molecule, analyze the SMIRKs of each torsion, check if the SMIRKs already reached a target count (5).
4. If not enough SMIRKs are in the selected_torsions yet, add this torsion to selected_torsions, save as JSON file.
4. Create dataset using the JSON file.

### Notes
 - The details of the torsion selection can be found in `select.log`
 - The SMIRKs target count of 5 is chosen to generate a reasonable number of torsions for scanning.

### Manifest
 - `primary_set.smi` - input SMILES file containing all molecules in benchmark primary set 1.
 - `01_generate.py` - script to generate OptimizationDataset inputs
 - `optimization_inputs.json` - input molecules
 - `cmiles_failures.json` - Molecules that rdkit failed to generate standardized tautomer
 - `omega_failures.json` - Molecules that Omega failed to generate conformers
 - `skipped_ions.json` - Molecules that were skipped (ions)
 - `requirements.txt` - versions of toolkits used for the python scripts
 - `02_select_torsions.py` - Script to select torsions from generated conformers
 - `select.log` - Stdout of the 02_select_torsions.py script for record
 - `03_create_torsiondrive_dataset.py` - Script for creating the TorsiondriveDataset
 - `04_visualize.py` - create a PDF of all species in the dataset, highlighting central torsion and atoms involved in the torsion
                    (copied from "2019-09-07-Pfizer-discrepancy-torsion-dataset-1")



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