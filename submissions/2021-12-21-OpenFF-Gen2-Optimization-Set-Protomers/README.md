# OpenFF Gen 2 Optimization Set Protomers

## General Information
 - Date: 12/21/2021
 - Class: Forcefield Parameterization
 - Purpose: Fill in Gen 2 optimization datasets with missing reasonable protomers
 - Collection: OptimizationDataset
 - Name: OpenFF Gen2 Optimization Dataset Protomers v1.0
 - Number of unique molecules        110
 - Number of filtered molecules      0
 - Number of conformers              610
 - Number of conformers min mean max 1   5.55 10
 - Mean molecular weight: 280.17
 - Max molecular weight: 542.59
 - Charges: [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0]
 - Submitter: Pavan Behara

## Changelog

Here any information regarding dataset changes are recorded.

## QCSubmit Manifest

1. check_protonation_states.py - code that checks the reasonable protonation states of molecules in Gen 2 optimization datasets and writes down the missing states to a file
2. confs_with_different_protonation_states.smi - smiles files that has the smiles patterns of missing reasonable protonation states
3. Dataset_Generation.ipynb - notebook to read the smiles file and generate optimization dataset
4. dataset.json.bz2 - generated dataset file
5. protomers_not_present_in_gen2_optimization_sets.pdf - visualizing the molecules in dataset.json.bz2

## Metadata

- elements {'O', 'F', 'S', 'Br', 'Cl', 'C', 'P', 'H', 'I', 'N'}
- unique molecules 110
- scf properties:
    - dipole
    - quadrupole
    - wiberg_lowdin_indices
    - mayer_indices
    - mbis_charges
- qc spec
    - name: default
    - method: b3lyp-d3bj
    - basis: dzvp