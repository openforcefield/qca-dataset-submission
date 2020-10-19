### Description

This dataset is a collection of CHO molecules.
The molecules are from the AlkEthOH and PhEthOH datasets originally used to build the smirnoff99Frosst parameters. The AlkEthOH was taken from [https://github.com/openforcefield/open-forcefield-data/tree/master/Model-Systems/AlkEthOH_distrib]()

The AlkEthOH SMILES contained some undefined stereo centers.
When this was encountered, all stereoisomers were generated and submitted for calculation.
Conformers were generated using a RMS cutoff of 3 Ångstroms.

### General Information
 - Date: 2020.09.18
 - Class: OpenFF Optimization dataset
 - Purpose: A set of CHO molecules to test force field typing schemes.
 - Collection: OptimizationDataset
 - Name: OpenFF Sandbox CHO PhAlkEthOH v1.0
 - Number of unique molecules: 7408
 - Number of unique stereoisomers: 10505
 - Number of unique conformations: 12592
 - Number of tasks submitted: 12271
 - Submitter: Trevor Gokey
 
### QCSubmit generation pipeline

```
./00.preprocess_input_data.sh
./01.run-qcsubmit.ipynb
```

The pipeline starts with preprocessing PhEthOH zip file to generate the SMILES files from the contained mol2 files.
Then, QCSubmit was used to filter duplicate SMILES, enumerate stereoisomers, and generate conformers. 

### QCSubmit Manifest
 
- `smi/AlkEthOH_chain.smi.bz2`
- `smi/AlkEthOH_rings.smi.bz2`
- `smi/PhEthOH.smi.bz2`: The initial molecule inputs containing the SMILES used for submission.
- `smi/PhEthOH/PhEthOH_mol2files_frosst-20200919T083119Z-001.zip`: The mol2 files containing the PhEthOh molecules.
- `00.preprocess_input_data.sh`: The commands needed to convert the extract and convert the mol2 files to create `PhEthOH.smi`
- `01.run-qcsubmit.ipynb`: The commands needed setup and run QCSubmit to create the submission.
- `lib/process_PheEthOH_mol2_to_smi.sh`: The helper command to extract and convert the mol2 files to SMILES
- `lib/to_smi.sh`: Converts a mol2 using `tleap` and `obabel`.
- `conda-env.yml`: The conda environment to run the preparation steps

### Metadata

```
{
	"submitter": "trevorgokey",
	"creation_date": "2020-09-18",
	"collection_type": "OptimizationDataset",
	"dataset_name": "OpenFF Sandbox CHO PhAlkEthOH v1.0", 
	"short_description": "A diverse set of CHO molecules",
	"long_description_url": "https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2020-09-18-OpenFF-Sandbox-CHO-PhAlkEthOH",
	"long_description": "This dataset contains a stereo-expanded version of the AlkEthOH dataset, and the original PhEthOH dataset, which were used in the original derivation of the smirnoff99Frosst parameters.",
	"elements": [
		"C",
		"H",
		"O"
	]
}
```
