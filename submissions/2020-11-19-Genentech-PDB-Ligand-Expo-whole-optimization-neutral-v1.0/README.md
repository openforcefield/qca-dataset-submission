### Description

This is the first of the Genentech PDB Ligand Expo Dataset referenced here (https://github.com/openforcefield/qca-dataset-submission/pull/48). This contains neutral molecules and with number of rotors less than or equal to three. Filtering for rotors brings down the original set to nearly 20% of its size (considering the initial number of molecules before conformer generation ).
Conformers were generated using a RMS cutoff of 3 Ångstroms.

### General Information
 - Date: 2020.11.19
 - Class: OpenFF Optimization dataset
 - Purpose: Improve coverage of FF
 - Collection: OptimizationDataset
 - Name: Genentech PDB Ligand Expo whole optimization neutral v1.0
 - Number of unique molecules: 445 
 - Number of unique conformers: 445
 - Number of tasks submitted: 445
 - Submitter: Pavan Behara
 
### QCSubmit generation pipeline

```
./01.run-qcsubmit.ipynb
```
`pubLigsNeutralGoodDensity.sdf` contains the sdf file  
QCSubmit was used to filter duplicate SMILES, filter rotors > 3, enumerate stereoisomers, and to generate conformers. 

### QCSubmit Manifest
 
- `pubLigsNeutralGoodDensity.sdf`: The initial molecule inputs containing the 3D structures
- `01.run-qcsubmit.ipynb`: The commands needed setup and run QCSubmit to create the submission.

### Metadata

```
{
	"submitter": "pavankum",
	"creation_date": "2020-11-19",
	"collection_type": "OptimizationDataset",
	"dataset_name": "Genentech PDB Ligand Expo whole optimization neutral v1.0", 
	"short_description": "",
	"long_description_url": "https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2020-11-19-Genentech-PDB-Ligand-Expo-whole-optimization-neutral-v1.0",
	"long_description": "This contains neutral molecules from Genentech PDB Ligand Expo dataset filtering molecules with rotors > 3.",
	"elements": [
	]
}
```

