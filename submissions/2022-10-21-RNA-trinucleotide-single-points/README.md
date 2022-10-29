## Description

This is a single point energy calculations of RNA basepairs and triple bases. Detailed description on how the data is generated can be found at https://github.com/choderalab/rna_bgsu.

## General Information

 - Date: 2022.10.21
 - Class: Basic dataset 
 - Purpose: Energy calculation
 - Collection: BasicDataset
 - Name: RNA Trinucleotide Single Point Dataset v1.0
 - Number of unique molecules:        64
 - Number of filtered molecules:      0
 - Number of conformers:              40835
 - Number of conformers min mean max: 185 638.05 1622
 - Mean molecular weight: 900.35
 - Max molecular weight: 971.64
 - Set of Charges: -2.0
 - Dataset Submitter: Kenichiro Takaba
 - Dataset Generator: Kenichiro Takaba
 - Dataset Source: https://github.com/choderalab/rna-bgsu-trinucleotide

## Changelog

Here any information regarding dataset changes are recorded.

## QCSubmit generation pipeline

 - `generate-dataset.ipynb`: This notebook shows how the dataset was prepared from the input files. 
 
## QCSubmit Manifest

- `generate-dataset.ipynb`: Dataset creation notebook with instructions for submission.
- `rdmols_loopMOTIFS.pkl.tar.gz`: The compressed pickle file containing 40835 rdkit molecules which are extracted from internal/hairpin/junction loop motifs.
- `dataset.json.bz`: The compressed basic dataset ready for submission.
- `dataset.smi`: The smiles file of the nucleic acids.
- `drawimage.py`: Notebook used to save 2D molecule structures used for submission.
- `dataset.pdf`: A pdf file containing molecule 2D structures.
- `pdb/pdb.tar.gz`: The compressed pdb structures that were used to generate the rdkit molecule objects which are saved as pickle files.

## Metadata

- elements {'O', 'N', 'C', 'H', 'P'}
- unique molecules 64
- Spec: default
    - scf properties:
        - dipole
        - quadrupole
    - qc spec
        - name: default
        - method: B3LYP-D3BJ
        - basis: DZVP
    
