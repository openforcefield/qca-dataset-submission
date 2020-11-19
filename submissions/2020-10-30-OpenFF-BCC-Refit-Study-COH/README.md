# OpenFF BCC Refit Study COH

### Description

A data set curated for the initial stage of the on-going OpenFF study which aims to co-optimize the AM1BCC bond charge correction (BCC) parameters against an experimental training set of density and enthalpy of mixing data points and a QM training set of electric field data.

The initial data set is limited to only molecules composed of C, O, H. This limited scope significantly reduces the number of BCC parameters which must be retrained, thus allowing for easier convergence of the initial optimizations.

The included molecules are those included in the experimental data set as well as an additional set chosen to ensure that each BCC parameter to train has been sufficiently (at least five instances) represented and exercised.

The conformers included in the set where generated using version 0.0.1a4 of the openff-recharge package. The exact conformer generation settings are attached as provenance.

### General Information

 - Date: 2020.10.30
 - Class: OpenFF Basic Dataset
 - Purpose: C,H,O training data for BCC refits.
 - Collection: BasicDataset
 - Name: OpenFF BCC Refit Study COH v1.0
 - Number of unique molecules: 94
 - Number of records: 215
 - Submitter: Simon Boothroyd
 
### QCSubmit generation pipeline

 - `generate.py`: A script which shows how the dataset was prepared from the input files. 
 
### QCSubmit Manifest

- `generate.py`: The script used to prepare the dataset.
- `dataset.json.xz`: The basic dataset ready for submission.
- `dataset.pdf`: A pdf file containing molecule 2D structures.
- `molecules.smi`: SMILES for every molecule in the submission.
 
### Metadata

- elements {'C', 'H', 'O'}
- unique molecules: 94
- optimizations: 215
- scf properties:
    - `dipole`
    - `quadrupole`
    - `wiberg_lowdin_indices`
    - `mayer_indices`
- qc spec
    - name: resp-2-vacuum
    - method: pw6b95
    - basis: aug-cc-pV(D+d)Z
    - store_wavefunction: orbitals_and_eigenvalues
- qc spec
    - name: resp-2-water
    - method: pw6b95
    - basis: aug-cc-pV(D+d)Z
    - store_wavefunction: orbitals_and_eigenvalues
    - implicit_solvent:
        - cavity_Type: GePol
        - cavity_Area: 0.3
        - cavity_Scaling: true
        - cavity_RadiiSet: Bondi
        - medium_SolverType: CPCM
        - medium_Solvent: H2O
