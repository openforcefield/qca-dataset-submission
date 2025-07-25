# OpenFF Lipid Torsion Drives v1.0

### Description

A torsion drive data set created to improve the coverage of lipid-like headgroup and backbone parameters in Sage.

### General Information

- Date: 2025-07-25
- Class: OpenFF TorsionDrive Dataset
- Purpose: Improve lipid headgroup and backbone coverage in Sage
- Name: 2025-07-25-OpenFF-Lipid-Torsion-Drives-v1.0
- Number of unique molecules: 16
- Number of driven torsions: 78
- Number of filtered molecules: 0
- Number of conformers: 321
- Number of conformers min mean max: 1, 4.12, 5
- Set of charges: [-1.0, 0.0, 1.0]
- Mean molecular weight: 194.30
- Max molecular weight: 297.22
- Dataset Submitter: Julianne Hoeflich
- Dataset Generator: Julianne Hoeflich 


### QCSubmit generation pipeline

- `generate-dataset.py`: This script shows how the dataset wsa prepared from the input gile `input.smi`.
- The list of SMILES in `input.smi` was curated by hand to correspond to relevant lipid model compounds. 

### QCSubmit Manifest

- `input.smi`: Input SMILES strings for dataset molecules
- `generate-dataset.py`: Script describing dataset generation and submission
- `input-environment.yaml`: Environment file used to create Python environment for the script
- `full-environment.yaml`: Fully-resolved environment used to execute the script
- `td.toml`: Experimental input file for defining variables used throughout the QCA submission process

- `dataset.json.bz2`: The compressed dataset ready for submission.
- `dataset.pdf`: A pdf file containing molecule 2D structures.
- `dataset.smi`: SMILES for every molecule in the submission.
 
### Metadata
* Elements: {C, P, N, H, O}
* Spec: default
     * basis: DZVP
     * implicit_solvent: None
     * keywords: {}
     * maxiter: 200
     * method: B3LYP-D3BJ
     * program: psi4
    * SCF properties:
        * dipole
        * quadrupole
        * wiberg_lowdin_indices
        * mayer_indices