# OpenFF Lipid Torsion Drives v4.0

### Description

A torsiondrive dataset was created to improve coverage of lipid-like headgroup, tail, and backbone parameters in the Sage force field. The model compounds were selected to enhance parameterization relevant to lipid force fields and include: methylacetate, dimethylphosphate, ethyltrimethylammonium, choline, 2-hexene, isopropylbutyrate, propylbutyrate, propylmethylphosphate, an esterified glycerol analogue ((M)EGLY), an esterified glycerol-phosphate analogue ((M)PGLY), cis-2-hexene, cis-5-decene, cis-7-pentadecene, and a lauroyl (LA) tail residue. These molecules were manually selected based on literature precedent and were chosen to mimic the training molecules used in the CHARMM36 and Amber lipid force fields.

The intended use of this dataset is to improve Sage parameters relevant to lipid simulations. The dataset provides minimal coverage for glycerol backbone parameters, ester parameters, alkene parameters, and amine and phospahte headgroup parameters. This effort builds on previous datasets to curate lipid-relevant QM data, including: 
- 2024-07-17-OpenFF-Phosphate-Torsion-Drives-v1.0
- 2024-08-09-OpenFF-Alkane-Torsion-Drives-v1.0
- 2024-10-08-OpenFF-Lipid-Optimization-Training-Supplement-v1.0
- 2024-10-30-OpenFF-Lipid-Optimization-Benchmark-Supplement-v1.0

This dataset was generated at the B3LYP-D3BJ/DZVP level of theory using Psi4 with no implicit solvent and default specs. It includes molecules containing the elements C, H, O, N, and P, spanning formal charges of -1, 0, and +1. The molecular weights range up to 297.22 amu, with a mean of 194.30 amu.

### General Information

- Date: 2025-07-25
- Class: OpenFF TorsionDrive Dataset
- Purpose: Improve lipid headgroup and backbone coverage in Sage
- Name: OpenFF Lipid Torsion Drives v4.0
- Number of unique molecules: 16
- Number of driven torsions: 78
- Number of filtered molecules: 0
- Number of conformers: 321
- Number of conformers min mean max: 1, 4.12, 5
- Set of charges: [-1.0, 0.0, 1.0]
- Mean molecular weight: 194.30
- Set of charges: [-1.0, 0.0, 1.0]
- Mean molecular weight: 194.30
- Min molecular weight: 74.08
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

- `dataset.json.bz2`: The compressed dataset ready for submission.
- `dataset.pdf`: A pdf file containing molecule 2D structures.
- `output.smi`: SMILES for every molecule in the submission.
 
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
