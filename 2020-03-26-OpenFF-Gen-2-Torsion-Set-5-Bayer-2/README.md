# OpenFF Gen 2 Torsion Set 5 Bayer 2

### General Information
 - Date: 3/26/2020
 - Class: Forcefield Parameterization
 - Purpose: Design 2nd generation torsion dataset for valence parameter fitting
 - Collection: TorsiondriveDataset
 - Name: OpenFF Gen 2 Torsion Set 5 Bayer 2
 - Number of Entries: 219 1-D torsions 
 - Submitter: Hyesu Jang

### Generation procedure

1. Generate `bayer_gen2_torsiondrive_data.pickle` by downloading `OpenFF Gen 2 Torsion Set 5 Bayer`  from qcarchive.
2. Read input JSON file, `bayer_optimization_inputs.json` and list torsion parameters and effective rotations assigned to each parameter.
3. Cluster list of rotations for each parameter so that the list is clustered into at least two clusters. 
4. Using a simple randomized optimization, select one torsion per cluster. Reuse pre-calculated torsion scans if possible. 
5. Store the selected torsions into JSON file, `bayer_2_selected_torsions.json`.

### Notes

 - The details of the torsion selection can be found in `select.log`;
 - 17 torsion scans out of 219 torsions are supposed to be reused. 

### Manifest

 - `bayer_optimization_inputs.json`: input molecules with expanded states using fragmenter and cmiles 
 - `param_valence.offxml`: parameter file 
 - `01_torsion_data_generator.ipynb`: jupyter notebook used in the dataset generation
 - `select.log`: stdout of `01_torsion_data_generator.ipynb` for record
 - `02_create_torsiondrive_dataset.py`: script for creating the TorsiondriveDataset and submitting
 - `requirements.txt`: version of toolkits used in the dataset generation


### Usage

1. Executing jupyter notebook creates `bayer_2_selected_torsions.json` which stores the selected torsions.

2. Create dataset on QCFractal server.
    ```
    python 02_create_torsiondrive_dataset.py
    ```
