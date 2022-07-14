## Description

This is a single point energy calculations of RNA basepairs and triple bases. Detailed description on how the data is generated can be found at https://github.com/choderalab/rna_bgsu.

## General Information

 - Date: 2022.07.13
 - Class: Basic dataset 
 - Purpose: Energy calculation
 - Collection: BasicDataset
 - Name: RNA Single Point Dataset v1.0
 - Number of unique molecules:        94
 - Number of filtered molecules:      0
 - Number of conformers:              4489
 - Number of conformers min mean max: 1  47.76 174
 - Mean molecular weight: 833.83
 - Max molecular weight: 971.64
 - Set of Charges: [-2.0, 0.0]
 - Dataset Submitter: Kenichiro Takaba
 - Dataset Generator: Kenichiro Takaba
 - Dataset Source: https://github.com/choderalab/rna_bgsu

## Changelog

Here any information regarding dataset changes are recorded.

## QCSubmit generation pipeline

 - `generate-dataset.ipynb`: This notebook shows how the dataset was prepared from the input files. 
 
## QCSubmit Manifest

- `generate-datasetn.ipynb`: Dataset creation notebook with instructions for submission.
- `rdmols_loopMOTIFS.pkl`: The pickle file which contains 4056 rdkit molecules which are internal/hairpin/junction loop motifs.
- `rdmols_triplebaseDB_exemplar.pkl`: The pickle file which contains 295 rdkit molecules which are experimental base triples.
- `rdmols_basepairCATALOG.pkl`: The pickle file which contains 138 rdkit molecules which are RNA base pairs.
- `dataset.json.bz`: The compressed basic dataset ready for submission.
- `dataset.smi`: The smiles file of the nucleic acids.
- `drawimage.py`: Notebook used to save 2D molecule structures used for submission.
- `dataset.pdf`: A pdf file containing molecule 2D structures.
- `pdb/pdb.tar.gz`: The compressed pdb structures that were used to generate the rdkit molecule objects which are saved as pickle files.
- `mols.sdf.tar.gz`: A compressed sdf file converted from the pdb structures.

## Metadata

- elements {'P', 'N', 'O', 'C', 'H'}
- unique molecules 94
- Spec: default
    - scf properties:
        - dipole
        - quadrupole
        - wiberg_lowdin_indices
        - mayer_indices
    - qc spec
        - name: default
        - method: B3LYP-D3BJ
        - basis: DZVP
- Spec: wb97m-d3bj/def2-tzvppd
    - scf properties:
        - dipole
        - quadrupole
        - wiberg_lowdin_indices
        - mayer_indices
    - qc spec
        - name: wB97M-D3BJ/def2-TZVPPD
        - method: wB97M-D3BJ
        - basis: def2-TZVPPD
    - keywords
        - wcombine: False
    
