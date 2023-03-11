## Description

This is a single point energy calculations of RNA nucleosides without O5' hydroxyl atom generated from 500K implicit solvent MD and chi torsion scanning. Data generation details can be found at https://github.com/choderalab/create-rna-nucleoside-dataset.

## General Information

 - Date: 2023.03.09
 - Class: Basic dataset 
 - Purpose: Energy calculation
 - Collection: BasicDataset
 - Name: RNA Nucleoside Single Point Dataset v1.0
 - Number of unique molecules:        4
 - Number of filtered molecules:      0
 - Number of conformers:              9555
 - Number of conformers min mean max: 2382 2388.75 2393
 - Mean molecular weight: 243.48
 - Max molecular weight: 267.25
 - Set of Charges: 0.0
 - Dataset Submitter: Kenichiro Takaba
 - Dataset Generator: Kenichiro Takaba
 - Dataset Source: https://github.com/choderalab/create-rna-nucleoside-dataset

## Changelog

Here any information regarding dataset changes are recorded.

## QCSubmit generation pipeline

 - `generate-dataset.ipynb`: This notebook shows how the dataset was prepared from the input files. 
 
## QCSubmit Manifest

- `generate-dataset.ipynb`: Dataset creation notebook with instructions for submission.
- `dataset.json.bz`: The compressed basic dataset ready for submission.
- `dataset.smi`: The smiles file of the RNA nucleosides.
- `dataset.pdf`: A pdf file containing molecule 2D structures.
- `sdf/torsion_scan_a_filtered.sdf.gz`: The compressed sdf for adenosine nucleosides (2393 conformers).
- `sdf/torsion_scan_c_filtered.sdf.gz`: The compressed sdf fof cytidine nucleosides (2382 conformers).
- `sdf/torsion_scan_g_filtered.sdf.gz`: The compressed sdf fof guanosine nucleosides (2392 conformers).
- `sdf/torsion_scan_u_filtered.sdf.gz`: The compressed sdf fof uridine nucleosides (2388 conformers).


## Metadata

- elements {'H', 'N', 'O', 'C''}
- unique molecules 4
- Spec: default
    - scf properties:
        - dipole
        - quadrupole
    - qc spec
        - name: default
        - method: B3LYP-D3BJ
        - basis: DZVP
    
