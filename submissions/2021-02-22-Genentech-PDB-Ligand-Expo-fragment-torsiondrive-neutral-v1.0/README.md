### Description

This is the third of the Genentech PDB Ligand Expo Datasets referenced here (https://github.com/openforcefield/qca-dataset-submission/pull/48).
This contains torsiondrives for fragmented neutral molecules.
Conformers were generated using a RMS cutoff of 3 Ångstroms.

### General Information
 - Date: 2021.02.22
 - Class: OpenFF TorsionDrive dataset
 - Purpose: Improve coverage of FF
 - Collection: TorsionDriveDataset
 - Name: Genentech PDB Ligand Expo fragment torsiondrive neutral v1.0
 - Number of unique molecules: --
 - Number of unique conformers: --
 - Number of tasks submitted: --
 - Submitter: Pavan Behara
 
### QCSubmit generation pipeline

```
./01.run-qcsubmit.ipynb
```
`pubLigsNeutralGoodDensity.sdf` contains the sdf file  
QCSubmit was used to filter duplicate SMILES, enumerate stereoisomers, and to generate conformers. 

### QCSubmit Manifest
 
- `pubLigsNeutralGoodDensity.sdf`: The initial molecule inputs containing the 3D structures
- `01.run-qcsubmit.ipynb`: The commands needed setup and run QCSubmit to create the submission.

### Metadata

```
{
	"submitter": "pavankum",
	"creation_date": "2021-02-22",
	"collection_type": "TorsionDriveDataset",
	"dataset_name": "Genentech PDB Ligand Expo fragment torsiondrive neutral v1.0", 
	"short_description": "",
	"long_description_url": "https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2021-02-22-Genentech-PDB-Ligand-Expo-fragment-torsiondrive-neutral-v1.0",
	"long_description": "This dataset contains torsiondrives for Genentech PDB Ligand Expo neutral set of fragmented molecules",
	"elements": [C ,Cl ,F ,N ,O ,I ,S ,H ,Br]
}
```

### Compute specs

```
{
  "qc_specifications": {
    "default-dlc": {
      "method": "B3LYP-D3BJ",
      "basis": "DZVP",
      "program": "psi4",
      "spec_name": "default-dlc",
      "spec_description": "Standard OpenFF optimization quantum chemistry specification, with DLC internal coordinates for geomeTRIC.",
      "store_wavefunction": "none",
      "implicit_solvent": null
    }
```
