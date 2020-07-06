### Description

This is the initial test of running constrained optimizations on various protein fragments prepared by David Cerutti.
Here we just have ALA as the central residue, once this is confirmed to run we will expand to all combinations.
Details from David 
The way it's set up, we are scanning phi and psi of the central residue with a random selection of ALA, GLY, SER, or VAL to the N- or C-terminus, and the customary ACE and NME blocking groups outside of that. We are targeting B3LYP, def2-dzvpp with Becke-Johnson D3 corrections for the optimization and thereafter to use B3LYP, ma-def2-tzvpp for the single point calculation getting the energy and gradient. That is "minimally augmented" def2-triple zeta / PP. 

### General Information
 - Date: 06/7/2020
 - Class: OpenFF constrained optimization 
 - Purpose: Optimizations 
 - Collection: OptimizationDataset
 - Name: OpenFF Protein Fragments v1.0
 - Number of unique molecules: 16
 - Submitter: Josh Horton
 
 ### QCSubmit generation pipeline
 - `Qcsubmit protein prep.ipynb`: This notebook shows how the optimization dataset was prepared from the input files and constraint files. 
 
 ### QCSubmit Manifest
- `Qcsubmit protein prep.ipynb`: Dataset creation notebook with instructions for submission.
- `protein_dataset.json`: The optimization dataset ready for submission.
- `protein_fragments.smi`: The smiles file of the biaryl molecules.
- `protein_dataset.pdf`: A pdf file containing molecule 2D structures.
- `*_ALA_*`: Folders containing the input molecule conformations and the corresponding dihedral restraints.
 
 ### Metadata
- elements {'O', 'C', 'H', 'N'}
- unique molecules 16
- optimizations 567
- scf properties:
    - dipole
    - quadrupole
    - wiberg_lowdin_indices
    - mayer_indices
- qc spec
    - name: default
    - method: B3LYP-D3BJ
    - basis: def2-dzvpp (to be confirmed)
    