# DANCE 1: eMolecules t142 selected

### General Information

- Date: 6/1/2020
- Class: Forcefield Parameterization
- Purpose: Use molecules selected from the eMolecules database by
  [DANCE](https://github.com/btjnaka/dance) to improve t142 parameterization in
  smirnoff99Frosst.
- Collection: TorsiondriveDataset
- Name: DANCE 1 eMolecules t142 selected
- Number of Entries: 20 1-D torsions
- Submitter: Bryon Tjanaka (Mobley Lab)

### Generation procedure

1. Output from DANCE is stored in `t142_selected.oeb` and `t142_selected.smi`.
   Run `python 01_generate.py` to turn these molecules into the input JSON file
   `t142_selected.json`.

### Notes

- The molecules were generated as described
  [here](https://github.com/btjanaka/dance/blob/c339368b398564cecce7e21cd88fdb3f3d2e363e/examples/t142-emolecules/README.md)
  (this link goes to the specific commit).

### Manifest

- `01_generator.py`: Python script used in the dataset generation
- `02_create_torsiondrive_dataset.py`: script for creating the
  TorsiondriveDataset and submitting
- `t142_selected.oeb`: Output from DANCE
- `t142_selected.smi`: Output from DANCE (same order as the OEB)
- `t142_selected_torsions.json`: Molecules generated from `01_generator.py`
- `requirements.txt`: version of toolkits used in the dataset generation
  (generated with `conda list | tail -n +3 > requirements.txt`)

### Usage

1. Executing `python 01_generate.py` generates `t142_selected.json` which stores
   the selected torsions.
2. Create dataset on QCFractal server.
   ```
   python 02_create_torsiondrive_dataset.py
   ```
