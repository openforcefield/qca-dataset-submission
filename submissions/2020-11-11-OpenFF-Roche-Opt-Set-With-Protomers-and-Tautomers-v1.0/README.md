### Description

This is an optimization dataset from Gen 1 fitting [https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2019-05-16-Roche-Optimization_Set]() with the protomers and tautomers enumerated.
Conformers were generated using a RMS cutoff of 3 Ångstroms.

Note: The 1.0 version of this dataset contained a fundamental error in indexing molecules, and so some molecules were excluded from the submission.
      The 1.1 version given here is the corrected version.

### General Information
 - Date: 2020.11.11
 - Class: OpenFF Optimization dataset
 - Purpose: Improve coverage of FF with addition of protonation states and tautomers
 - Collection: OptimizationDataset
 - Name: OpenFF Roche Opt Set With Protomers and Tautomers v1.1
 - Number of unique molecules: 1376
 - Number of unique conformers: 1376
 - Number of tasks submitted: 1376
 - Submitter: Pavan Behara
 
### QCSubmit generation pipeline

```
./01.run-qcsubmit.ipynb
```
`OpenFF_references.sdf.tar` contains the compressed sdf file with 468 molecules 
QCSubmit was used to filter duplicate SMILES, enumerate protomers, tautomers, stereoisomers, and to generate conformers. 

### QCSubmit Manifest
 
- `OpenFF_references.sdf.tar`: The initial molecule inputs containing the 3D structures
- `01.run-qcsubmit.ipynb`: The commands needed setup and run QCSubmit to create the submission.

### Metadata

```
{
	"submitter": "pavankum",
	"creation_date": "2020-11-11",
	"collection_type": "OptimizationDataset",
	"dataset_name": "OpenFF Roche Opt Set With Protomers and Tautomers v1.1", 
	"short_description": "",
	"long_description_url": "https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2020-11-11-OpenFF Roche Opt-Set-With-Protomers-and-Tautomers-v1.1",
	"long_description": "This dataset contains a tuatomer and protomer expanded version of the Roche Optimization set, which was used in Gen 1 fitting of FF.",
	"elements": [
		"Cl",
		"S",
		"C",
		"F",
		"H",
		"O",
		"N"
	]
}
```

