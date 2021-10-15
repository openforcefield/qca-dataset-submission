# OpenFF BCC Refit Study COH v2.0

### Description

A data set used for training ESP-fitting based typed atomic polarizabilities with a direct approximation.
Seven optimizations will be performed on each conformer: an unpolarized ESP and six polarized ESPs with imposed external electric fields.

This relatively small data set consisted of 39 small molecules that have commonly occurred functional groups in drug-like molecule, such as hydroxyl group, aromatic ring, carbonyl group, methyl group, and amines.
One positively charged and one negatively charged molecule are fitted simultaneously to capture the electronic polarization.

The initial data set is limited to only molecules composed of C, O, H, and N to enable an element-based and relatively fast test.

The v1.1 version of this dataset fixes incorrect values for the `PERTURB_DIPOLE` keyword; these are meant to be a list of floats, not a string.


### General Information

- Date: 2021.10.05
- Class: OpenFF Optimization Dataset
- Purpose: C,H,O,N training data for training ESP-fitting based atomic polarizabilities and partial charges
- Collection: OptimizationDataset
- Name: OpenFF RESP Polarizability Optimizations v1.1 
- Number of unique molecules        39
- Number of filtered molecules      0
- Number of conformers              105
- Number of conformers min mean max 1    2.69  5
- Dataset Submitter: Willa Wang 
- Dataset Generator: Willa Wang 
- Set of charges: [-1.0, 0.0, 1.0]
- Mean molecular weight: 90.54
- Max molecular weight: 158.24
- Enumerate stereoisomers: False
- Enumerate tautomers: False
- Enumerate protomers: False

### QCSubmit generation pipeline

- `generate-dataset.ipynb`: A notebook which shows how the dataset was prepared from the input files.
- `generate-dataset-v1.1.ipynb`: A notebook which shows how the dataset was prepared from the input files, for v1.1
- `molecules.smi`: A text file which includes the SMILES strings of input molecules.

### QCSubmit Manifest

- `generate-dataset.ipynb`
- `dataset.json.bz2`: The basic dataset ready for submission.
- `dataset.pdf`: A pdf file containing molecule 2D structures.
- `dataset.smi`: SMILES for every molecule in the submission.

Corresponding files for v1.1 also included.
 
### Metadata

```
{'collection_type': 'OptimizationDataset',
 'creation_date': datetime.date(2021, 10, 14),
 'dataset_name': 'OpenFF RESP Polarizability Optimizations v1.1',
 'elements': {'O', 'C', 'N', 'H'},
 'long_description': 'A data set used for training typed polarizabilities '
                     'using direct polarization.\n'
                     'This data set only includes element C, H, N, and O.',
 'long_description_url': HttpUrl('https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2021-10-01-OpenFF-resppol-mp2-single-point', scheme='https', host='github.com', tld='com', host_type='domain', path='/openforcefield/qca-dataset-submission/tree/master/submissions/2021-10-01-OpenFF-resppol-mp2-single-point'),
 'short_description': 'Optimizations of ESP-fitting based direct '
                      'polarizabilities.',
 'submitter': 'willawang'}
```

## QC specifications:

```
Spec: MP2/aug-cc-pVTZ/X-
{'basis': 'aug-cc-pVTZ',
 'implicit_solvent': None,
 'keywords': {'E_CONVERGENCE': '1.0e-8',
              'PERTURB_DIPOLE': [-0.01, 0.0, 0.0],
              'PERTURB_H': True,
              'PERTURB_WITH': 'DIPOLE',
              'mp2_type': 'df',
              'scf_type': 'df'},
 'maxiter': 200,
 'method': 'MP2',
 'program': 'psi4',
 'scf_properties': ['dipole',
                    'quadrupole',
                    'wiberg_lowdin_indices',
                    'mayer_indices'],
 'spec_description': 'The quantum chemistry specification used to generate '
                     'data for typed polarizabilities training.',
 'spec_name': 'MP2/aug-cc-pVTZ/X-',
 'store_wavefunction': 'orbitals_and_eigenvalues'}
Spec: MP2/aug-cc-pVTZ/X+
{'basis': 'aug-cc-pVTZ',
 'implicit_solvent': None,
 'keywords': {'E_CONVERGENCE': '1.0e-8',
              'PERTURB_DIPOLE': [0.01, 0.0, 0.0],
              'PERTURB_H': True,
              'PERTURB_WITH': 'DIPOLE',
              'mp2_type': 'df',
              'scf_type': 'df'},
 'maxiter': 200,
 'method': 'MP2',
 'program': 'psi4',
 'scf_properties': ['dipole',
                    'quadrupole',
                    'wiberg_lowdin_indices',
                    'mayer_indices'],
 'spec_description': 'The quantum chemistry specification used to generate '
                     'data for typed polarizabilities training.',
 'spec_name': 'MP2/aug-cc-pVTZ/X+',
 'store_wavefunction': 'orbitals_and_eigenvalues'}
Spec: MP2/aug-cc-pVTZ/Y-
{'basis': 'aug-cc-pVTZ',
 'implicit_solvent': None,
 'keywords': {'E_CONVERGENCE': '1.0e-8',
              'PERTURB_DIPOLE': [0.0, -0.01, 0.0],
              'PERTURB_H': True,
              'PERTURB_WITH': 'DIPOLE',
              'mp2_type': 'df',
              'scf_type': 'df'},
 'maxiter': 200,
 'method': 'MP2',
 'program': 'psi4',
 'scf_properties': ['dipole',
                    'quadrupole',
                    'wiberg_lowdin_indices',
                    'mayer_indices'],
 'spec_description': 'The quantum chemistry specification used to generate '
                     'data for typed polarizabilities training.',
 'spec_name': 'MP2/aug-cc-pVTZ/Y-',
 'store_wavefunction': 'orbitals_and_eigenvalues'}
Spec: MP2/aug-cc-pVTZ/Y+
{'basis': 'aug-cc-pVTZ',
 'implicit_solvent': None,
 'keywords': {'E_CONVERGENCE': '1.0e-8',
              'PERTURB_DIPOLE': [0.0, 0.01, 0.0],
              'PERTURB_H': True,
              'PERTURB_WITH': 'DIPOLE',
              'mp2_type': 'df',
              'scf_type': 'df'},
 'maxiter': 200,
 'method': 'MP2',
 'program': 'psi4',
 'scf_properties': ['dipole',
                    'quadrupole',
                    'wiberg_lowdin_indices',
                    'mayer_indices'],
 'spec_description': 'The quantum chemistry specification used to generate '
                     'data for typed polarizabilities training.',
 'spec_name': 'MP2/aug-cc-pVTZ/Y+',
 'store_wavefunction': 'orbitals_and_eigenvalues'}
Spec: MP2/aug-cc-pVTZ/Z-
{'basis': 'aug-cc-pVTZ',
 'implicit_solvent': None,
 'keywords': {'E_CONVERGENCE': '1.0e-8',
              'PERTURB_DIPOLE': [0.0, 0.0, -0.01],
              'PERTURB_H': True,
              'PERTURB_WITH': 'DIPOLE',
              'mp2_type': 'df',
              'scf_type': 'df'},
 'maxiter': 200,
 'method': 'MP2',
 'program': 'psi4',
 'scf_properties': ['dipole',
                    'quadrupole',
                    'wiberg_lowdin_indices',
                    'mayer_indices'],
 'spec_description': 'The quantum chemistry specification used to generate '
                     'data for typed polarizabilities training.',
 'spec_name': 'MP2/aug-cc-pVTZ/Z-',
 'store_wavefunction': 'orbitals_and_eigenvalues'}
Spec: MP2/aug-cc-pVTZ/Z+
{'basis': 'aug-cc-pVTZ',
 'implicit_solvent': None,
 'keywords': {'E_CONVERGENCE': '1.0e-8',
              'PERTURB_DIPOLE': [0.0, 0.0, 0.01],
              'PERTURB_H': True,
              'PERTURB_WITH': 'DIPOLE',
              'mp2_type': 'df',
              'scf_type': 'df'},
 'maxiter': 200,
 'method': 'MP2',
 'program': 'psi4',
 'scf_properties': ['dipole',
                    'quadrupole',
                    'wiberg_lowdin_indices',
                    'mayer_indices'],
 'spec_description': 'The quantum chemistry specification used to generate '
                     'data for typed polarizabilities training.',
 'spec_name': 'MP2/aug-cc-pVTZ/Z+',
 'store_wavefunction': 'orbitals_and_eigenvalues'}
Spec: MP2/aug-cc-pVTZ
{'basis': 'aug-cc-pVTZ',
 'implicit_solvent': None,
 'keywords': None,
 'maxiter': 200,
 'method': 'MP2',
 'program': 'psi4',
 'scf_properties': ['dipole',
                    'quadrupole',
                    'wiberg_lowdin_indices',
                    'mayer_indices'],
 'spec_description': 'The quantum chemistry specification used to generate '
                     'data for typed polarizabilities training.',
 'spec_name': 'MP2/aug-cc-pVTZ',
 'store_wavefunction': 'orbitals_and_eigenvalues'}
```
