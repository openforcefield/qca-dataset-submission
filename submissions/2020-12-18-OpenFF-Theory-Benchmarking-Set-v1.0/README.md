### Description

This is a TorsionDrive dataset containing 1-D torsions for theory benchmarking.

### Dataset Information

- Date: 2020-12-18
- Class: OpenFF TorsionDrive Dataset
- Name: OpenFF Theory Benchmarking Set v1.0
- Version: 1.0 
- Description: A torsiondrive dataset for theory benchmarking
- Long Description URL: https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2020-12-14-theory-bm-torsiondrive-set
- Changelog: { 1.0 :{ 'submitter' : 'Hyesu Jang', 'date': '2020-12-18', 'description': 'A torsiondrive dataset for theory benchmarking with additional dataset entries'}}
- Dataset generator: Hyesu Jang
- Dataset Submitter: Hyesu Jang
- Set of charges: {-1, 0, 1}
- Mean molecular weight: 168.00
- Max molecular weight: 233.29
- Enumerate stereoisomers: False
- Enumerate tautomers: False

### Torsion selection scheme
1. 36 1-D torsions have been initially selected from Roche and Coverage set based on two criteria:
    1. Consideration of the coverage of chemical diversity
        - Variations in central bonds, formal charges, element compositions, intramolecular interactions;
        - Inclusion of molecules (1) w/ non-zero formal charges, (2) w/ strong internal interactions, (3) w/ central bond conjugated( < 10 kcal/mol rotational barrier) or (4) w/ halogen
    2. Consistency in molecular size

2. Additional 1-D torsions have been selected to include more diverse charged molecules.
    1. From the recently isomerized Roche set by Pavan, extracted charged molecules.
    2. Separated the molecules into groups based on the type of charged functional group.
        - -1 charged functional groups: `c[O-]`, `C(=O)[N-]`, `c[N-]c`, `S(=O)(=O)[N-]`, `S(=[N-])(=O)`
        - +1 charged functional groups: `[NH+,nH+](=,:[C,c])[C,c]`, `[NH+]([*])[*]`, `[NH2+]([*])[*]`, `[NH3+][*]`
    3. Then selected one molecule per each group (by picking a center molecule using MACCS keys fingerprint.)
    4. Using qcsubmit tool, generated 1-D torsions from the selected molecules. 


## Manifest

- `QCSubmit workflow.ipynb`: curation of torsiondrive dataset using qcsubmit package.
- `dataset.json.bz2`: The torsiondrive dataset ready for submission.
- `compute*.json`: additional specs added after the first submission.
## Metadata

```
{'submitter': 'hyejang',
 'creation_date': '2020-12-14',
 'collection_type': 'TorsiondriveDataset',
 'dataset_name': 'OpenFF Theory Benchmarking Set v1.0',
 'short_description': 'Torsiondrives for theory benchmarking',
 'long_description_url': 'https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2020-12-14-theory-bm-torsiondrive-set-v1.0',
 'long_description': 'A torsiondrive dataset for theory benchmarking',
 'elements': ['C', 'N', 'Cl', 'F', 'O', 'S', 'P', 'H']}
```

## QC specification 

```
{'B3LYP-D3BJ/DZVP ': {'method': 'B3LYP-D3BJ',
  'basis': 'dzvp',
  'program': 'psi4',
  'spec_name': 'B3LYP-D3BJ/DZVP ',
  'spec_description': 'A torsiondrive dataset for benchmarking B3LYP-D3BJ/DZVP ',
  'store_wavefunction': 'none',
  'implicit_solvent': None},
 'B3LYP-D3BJ/DEF2-TZVP ': {'method': 'B3LYP-D3BJ',
  'basis': 'def2-tzvp',
  'program': 'psi4',
  'spec_name': 'B3LYP-D3BJ/DEF2-TZVP ',
  'spec_description': 'A torsiondrive dataset for benchmarking B3LYP-D3BJ/DEF2-TZVP ',
  'store_wavefunction': 'none',
  'implicit_solvent': None},
 'B3LYP-D3BJ/DEF2-TZVPD ': {'method': 'B3LYP-D3BJ',
  'basis': 'def2-tzvpd',
  'program': 'psi4',
  'spec_name': 'B3LYP-D3BJ/DEF2-TZVPD ',
  'spec_description': 'A torsiondrive dataset for benchmarking B3LYP-D3BJ/DEF2-TZVPD ',
  'store_wavefunction': 'none',
  'implicit_solvent': None},
 'B3LYP-D3BJ/DEF2-TZVPP ': {'method': 'B3LYP-D3BJ',
  'basis': 'def2-tzvpp',
  'program': 'psi4',
  'spec_name': 'B3LYP-D3BJ/DEF2-TZVPP ',
  'spec_description': 'A torsiondrive dataset for benchmarking B3LYP-D3BJ/DEF2-TZVPP ',
  'store_wavefunction': 'none',
  'implicit_solvent': None},
 'B3LYP-D3BJ/DEF2-TZVPPD ': {'method': 'B3LYP-D3BJ',
  'basis': 'def2-tzvppd',
  'program': 'psi4',
  'spec_name': 'B3LYP-D3BJ/DEF2-TZVPPD ',
  'spec_description': 'A torsiondrive dataset for benchmarking B3LYP-D3BJ/DEF2-TZVPPD ',
  'store_wavefunction': 'none',
  'implicit_solvent': None},
 'B3LYP-D3BJ/DEF2-QZVP ': {'method': 'B3LYP-D3BJ',
  'basis': 'def2-qzvp',
  'program': 'psi4',
  'spec_name': 'B3LYP-D3BJ/DEF2-QZVP ',
  'spec_description': 'A torsiondrive dataset for benchmarking B3LYP-D3BJ/DEF2-QZVP ',
  'store_wavefunction': 'none',
  'implicit_solvent': None},
 'B3LYP-D3BJ/6-31+GSS ': {'method': 'B3LYP-D3BJ',
  'basis': '6-31+gss',
  'program': 'psi4',
  'spec_name': 'B3LYP-D3BJ/6-31+GSS ',
  'spec_description': 'A torsiondrive dataset for benchmarking B3LYP-D3BJ/6-31+GSS ',
  'store_wavefunction': 'none',
  'implicit_solvent': None},
 'B3LYP-D3BJ/6-311+GSS ': {'method': 'B3LYP-D3BJ',
  'basis': '6-311+gss',
  'program': 'psi4',
  'spec_name': 'B3LYP-D3BJ/6-311+GSS ',
  'spec_description': 'A torsiondrive dataset for benchmarking B3LYP-D3BJ/6-311+GSS ',
  'store_wavefunction': 'none',
  'implicit_solvent': None}}
```

## Provenance 

```
{'qcsubmit': '0.1.1', 'openforcefield': '0.8.1', 'openeye': '2020.1.1'}
```
