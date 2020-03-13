# OpenFF Gen 2 Torsion Set 4 eMolecules Discrepancy

### General Information
 - Date: 3/12/2020
 - Class: Forcefield Parameterization
 - Purpose: Design 2nd generation torsion dataset for valence parameter fitting
 - Collection: TorsiondriveDataset
 - Name: OpenFF Gen 2 Torsion Set 4 eMolecules Discrepancy
 - Number of Entries: 106 1-D torsions 
 - Submitter: Hyesu Jang

### Generation procedure

1. Read input JSON file and list torsion parameters and effective rotations matched to each parameter.
2. Select torsions (one torsion per parameter) in a way to minimize data degeneracy.
3. Store selected torsions into JSON file.

### Notes
 - The details of the torsion selection can be found in `select.log`;

### Manifest

 - `emolecules_optimization_inputs.json`: input molecules with expanded states using fragmenter and cmiles 
 - `param_valence.offxml`: parameter file 
 - `01_torsion_data_generator.ipynb`: jupyter notebook used in the dataset generation
 - `select.log`: stdout of `01_torsion_data_generator.ipynb` for record
 - `02_create_torsiondrive_dataset.py`: script for creating the TorsiondriveDataset and submitting
 - `requirements.txt`: version of toolkits used in the dataset generation


### Usage

1. Executing jupyter notebook creates `emolecules_selected_torsions.json` which stores the selected torsions.

2. Create dataset on QCFractal server.
    ```
    python 02_create_torsiondrive_dataset.py
    ```
