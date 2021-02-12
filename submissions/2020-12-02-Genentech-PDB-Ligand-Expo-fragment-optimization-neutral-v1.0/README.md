### Description

This is the second of the Genentech PDB Ligand Expo Dataset referenced here (https://github.com/openforcefield/qca-dataset-submission/pull/48). This contains neutral molecules and with number of rotors more than three. These molecules are fragmented.
Conformers were generated using a RMS cutoff of 3 Ångstroms.

### General Information
 - Date: 2020.12.02
 - Class: OpenFF Optimization dataset
 - Purpose: Improve coverage of FF
 - Collection: OptimizationDataset
 - Name: Genentech PDB Ligand Expo fragment optimization neutral v1.0
 - Number of unique molecules: 2319
 - Number of unique conformers: 2366
 - Number of tasks submitted: 2363
 - Submitter: Pavan Behara
 
### QCSubmit generation pipeline

```
./01.run-qcsubmit.ipynb
```
`pubLigsNeutralGoodDensity.sdf` contains the sdf file  
QCSubmit was used to filter duplicate SMILES, filter rotors < 3, enumerate stereoisomers, and to generate conformers. 

### QCSubmit Manifest
 
- `pubLigsNeutralGoodDensity.sdf`: The initial molecule inputs containing the 3D structures
- `01.run-qcsubmit.ipynb`: The commands needed setup and run QCSubmit to create the submission.

### Metadata

```
{
	"submitter": "pavankum",
	"creation_date": "2020-12-02",
	"collection_type": "OptimizationDataset",
	"dataset_name": "Genentech PDB Ligand Expo fragment optimization neutral v1.0", 
	"short_description": "",
	"long_description_url": "https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2020-12-02-Genentech-PDB-Ligand-Expo-fragment-optimization-neutral-v1.0",
	"long_description": "This dataset contains Genentech PDB Ligand Expo neutral set of molecules that fragments molecules with greater than 3 rotors",
	"elements": [C ,Cl ,F ,N ,O ,I ,S ,H ,Br]
}
```

