### Description

This is a TorsionDrive dataset consisting of 36 1-D torsions selected for benchmarking different QM levels.


### Torsion selection scheme 

1. Consideration of the coverage of chemical diversity
    - Variations in central bonds, formal charges, element compositions, intramolecular interactions;
    - Inclusion of molecules (1) w/ non-zero formal charges, (2) w/ strong internal interactions, (3) w/ central bond conjugated( < 10 kcal/mol rotational barrier) or (4) w/ halogen
    
2. Consistency in molecular size


### General Information 

- Date: 2020-10-29
- Class: OpenFF torsiondrive 
- Purpose: a torsiondrive dataset calculated with B3LYP-D3BJ/6-31+G** for benchmarking different QM levels 

- Collection: TorsionDriveDataset
- Name: OpenFF Theory Benchmarking Set B3LYP-D3BJ 6-31+Gss
- Number of Entries: 36
- Submitter: Hyesu Jang
 

### Manifest

- `QCSubmit workflow.ipynb`: curation of hand-picked torsions using qcsubmit package.
- `input_torsions.json`: QCArchive entry of hand-picked torsions
- `dataset.json`: The torsiondrive dataset ready for submission.
- `theory-bm-set_settings.yaml`: The QCSubmit settings used to generate the torsiondrive dataset.
- `theory-bm-set-curated.smi`:  The smiles file of the dataset.
- `theory-bm-set-curated.pdf`: A pdf file containing molecule 2D structures with targeted torsions for driving highlighted.


### Metadata

- elements {'C', 'N', 'Cl', 'H', 'P', 'O', 'S', 'F'}
- unique molecules 31
- torsiondrives 36
- scf properties:
    - dipole
    - quadrupole
    - wiberg_lowdin_indices
    - mayer_indices
- qc spec
    - name: default
    - method: B3LYP-D3BJ
    - basis: 6-31+G**
