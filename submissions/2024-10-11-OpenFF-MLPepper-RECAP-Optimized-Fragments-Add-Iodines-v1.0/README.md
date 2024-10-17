# OpenFF BCC Refit Study COH v2.0

### Description

A single point dataset created by combining the [50k ESP from Simon](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2022-01-16-OpenFF-ESP-Fragment-Conformers-v1.0) and 
[Br substituted set from Lily](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2023-11-30-OpenFF-multi-Br-ESP-Fragment-Conformers-v1.1-single-point), filtereing by Cl and Br and replacing them successively with iodines:

```
with open("./all_dataset.smi", "r") as input:
    for mols in input:
        smiles.append(mols.strip('\n'))

for molecule in smiles:
    if 'Cl' in molecule:
        iodine_smiles.append(molecule.replace('Cl','I'))
    if 'Br' in molecule:
        iodine_smiles.append(molecule.replace('Br','I'))
        
with open('./iodine_smiles.smi','w+') as file:
    for iodine in iodine_smiles:
        file.write(iodine + '\n')

```

Each fragment had 5 conformations generated which were optimised locally using an AIMNET2 model trained to `wb97m-d3`. 
This adds an extension of iodines molecules to the [original mlpepper dataset](https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2024-07-26-MLPepper-RECAP-Optimized-Fragments-v1.0).

The aim of the dataset is to provide polarised and gas phase electrostatic properties which can be used to generate ML models 
for partial charge prediction. Unlike past datasets the wavefunction will not be saved to recompute the ESP instead we recommend building the ESP 
from the MBIS atomic multipoles which save substantial amount of space. 

An off equilibrium data set will also be generated to enable conformation dependent prediction of charges. 

### General Information

* Number of conformers: 6041
* Number of conformers per molecule (min, mean, max): 1, 1.08, 4
* Mean molecular weight: 280.08
* Max molecular weight: 701.59
* Charges: [-4.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0]
## Metadata
* Elements: {N, Si, O, I, Br, F, Cl, B, C, P, S, H}
* Spec: wb97x-d/def2-tzvpp
	* basis: def2-tzvpp
	* implicit_solvent: None
	* keywords: {'dft_spherical_points': 590, 'dft_radial_points': 99, 'debug': 1}
	* maxiter: 200
	* method: wb97x-d
	* program: psi4
	* SCF properties:
		* dipole
		* quadrupole
		* lowdin_charges
		* mulliken_charges
		* mbis_charges
		* mayer_indices
		* wiberg_lowdin_indices
		* dipole_polarizabilities
* Spec: wb97x-d/def2-tzvpp/ddx-water
	* basis: def2-tzvpp
	* implicit_solvent: {'ddx_model': 'pcm', 'ddx_radii_scaling': 1.1, 'ddx_radii_set': 'uff', 'ddx_solvent_epsilon': 78.4, 'ddx_solvent': 'water'}
	* keywords: {'dft_spherical_points': 590, 'dft_radial_points': 99, 'debug': 1}
	* maxiter: 200
	* method: wb97x-d
	* program: psi4
	* SCF properties:
		* dipole
		* quadrupole
		* lowdin_charges
		* mulliken_charges
		* mbis_charges
		* mayer_indices
		* wiberg_lowdin_indices
		* dipole_polarizabilities

### QCSubmit generation pipeline

- `generate-dataset.ipynb`: A notebook which shows how the dataset was prepared from the input files.

### QCSubmit Manifest

- `generate-dataset.ipynb`
- `dataset.json.bz2`: The basic dataset ready for submission.
- `dataset.pdf`: A pdf file containing molecule 2D structures.
- `dataset.smi`: SMILES for every molecule in the submission.
 
### Metadata

* Elements: {Si, F, Cl, C, H, S, P, N, I, O, B, Br}
* Spec: wb97x-d/def2-tzvpp
        * basis: def2-tzvpp
        * implicit_solvent: None
        * keywords: {'dft_spherical_points': 590, 'dft_radial_points': 99, 'debug': 1}
        * maxiter: 200
        * method: wb97x-d
        * program: psi4
        * SCF properties:
                * dipole
                * quadrupole
                * lowdin_charges
                * mulliken_charges
                * mbis_charges
                * mayer_indices
                * wiberg_lowdin_indices
                * dipole_polarizabilities
* Spec: wb97x-d/def2-tzvpp/ddx-water
        * basis: def2-tzvpp
        * implicit_solvent: {'ddx_model': 'pcm', 'ddx_radii_scaling': 1.1, 'ddx_radii_set': 'uff', 'ddx_solvent_epsilon': 78.4, 'ddx_solvent': 'water'}
        * keywords: {'dft_spherical_points': 590, 'dft_radial_points': 99, 'debug': 1}
        * maxiter: 200
        * method: wb97x-d
        * program: psi4
        * SCF properties:
                * dipole
                * quadrupole
                * lowdin_charges
                * mulliken_charges
                * mbis_charges
                * mayer_indices
                * wiberg_lowdin_indices
                * dipole_polarizabilities
```
{'collection_type': 'OptimizationDataset',
 'creation_date': datetime.date(2021, 7, 2),
 'dataset_name': 'OpenFF BCC Refit Study COH v2.0',
 'elements': {'H', 'O', 'C'},
 'long_description': 'A data set curated for the initial stage of the on-going '
                     'OpenFF study which aims to co-optimize the AM1BCC bond '
                     'charge correction (BCC) parameters against an '
                     'experimental training set of density and enthalpy of '
                     'mixing data points and a QM training set of electric '
                     'field data.\n'
                     '\n'
                     'The initial data set is limited to only molecules '
                     'composed of C, O, H. This limited scope significantly '
                     'reduces the number of BCC parameters which must be '
                     'retrained, thus allowing for easier convergence of the '
                     'initial optimizations.\n'
                     '\n'
                     'The included molecules were combinatorially generated to '
                     'cover a range of alcohol, ether, and carbonyl containing '
                     'molecules.',
 'long_description_url': HttpUrl('https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2021-06-22-OpenFF-BCC-Refit-Study-COH-v2.0', scheme='https', host='github.com', tld='com', host_type='domain', path='/openforcefield/qca-dataset-submission/tree/master/submissions/2021-06-22-OpenFF-BCC-Refit-Study-COH-v2.0'),
 'short_description': 'Optimizations of diverse, para-substituted aniline '
                      'derivatives.',
 'submitter': 'simonboothroyd'}
```

## QC specifications:

<!-- copy spec outputs from generation notebook here -->

```

```

## SCF Properties

<!-- copy scf property outputs from generation notebook here -->

```
[<SCFProperties.Dipole: 'dipole'>,
 <SCFProperties.Quadrupole: 'quadrupole'>,
 <SCFProperties.WibergLowdinIndices: 'wiberg_lowdin_indices'>,
 <SCFProperties.MayerIndices: 'mayer_indices'>]
```
