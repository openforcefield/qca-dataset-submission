# DANCE 1: eMolecules t142 selected

### General Information

- Date: 6/4/2020
- Class: Forcefield Parameterization
- Purpose: Use molecules selected from the eMolecules database by
  [DANCE](https://github.com/btjnaka/dance) to improve t142 parameterization in
  smirnoff99Frosst.
- Collection: TorsiondriveDataset
- Name: OpenFF DANCE 1 eMolecules t142 v1.0
- Number of Entries: 20 1-D torsions
- Submitter: Bryon Tjanaka (Mobley Lab)

### Generation procedure

1. Output from DANCE is stored in `t142_selected.smi`. Run
   `python 01_generate.py` to turn these molecules into the input JSON file
   `t142_selected.json.gz`. The indices of the `t142` parameter are
   re-calculated while doing this and stored in the `atom_indices` field in the
   JSON file.

### Notes

- The molecules were generated as described
  [here](https://github.com/btjanaka/dance/blob/c339368b398564cecce7e21cd88fdb3f3d2e363e/examples/t142-emolecules/README.md)
  (this link goes to the specific commit).

### Manifest

- `01_generator.py`: Python script used in the dataset generation
- `t142_selected.smi`: Output from DANCE (same order as the OEB)
- `optimization_inputs.json.gz`: Molecules generated from `01_generator.py`
- `02_create_torsiondrive_dataset.py`: script for creating the
  TorsiondriveDataset and submitting
- `requirements.txt`: version of toolkits used in the dataset generation
  (generated with `conda list | tail -n +3 > requirements.txt`)

### Usage

1. Executing `python 01_generate.py` generates `t142_selected.json` which stores
   the selected torsions.
2. Create dataset on QCFractal server.
   ```
   python 02_create_torsiondrive_dataset.py
   ```
