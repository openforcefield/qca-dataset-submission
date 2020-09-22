### Description

This dataset is a collection of CHO molecules. The molecules are the AlkEthOH and PhEthOH datasets originally used to build the smirnoff99Frosst parameters.

### General Information
 - Date: 09/18/2020
 - Class: OpenFF Optimization dataset
 - Purpose: A set of CHO molecules to test force field typing schemes.
 - Collection: OptimizationDataset
 - Name: OpenFF Sandbox CHO PhAlkEthOH v1.0
 - Number of unique molecules: 7541
 - Number of unique stereoisomers: 10638
 - Number of unique conformations: 40133
 - Submitter: Trevor Gokey
 
### QCSubmit generation pipeline

```
./00.preprocess_input_data.sh
./01.make_input_molecules.sh
./02.submit_dataset.sh
```

The pipeline starts with preprocessing PhEthOH zip file to generate the SMILES files from the contained mol2 files. Then, the next step enumerates and generates the SMILES to generate conformers. Finally, it takes the `metadata.json`, `compute.json`, and the molecule JSON files and submits the dataset.
 
### QCSubmit Manifest
 
- `compute.json`: The file containing the compute specifications.
- `metadata.json`: The file containing the dataset metadata.
- `smi/AlkEthOH_chain.smi`
- `smi/AlkEthOH_rings.smi`
- `smi/PhEthOH.smi`: The initial molecule inputs containing the SMILES used for submission.
- `smi/AlkEthOH_rings.smi.json.lzma`
- `smi/AlkEthOH_chain.smi.json.lzma`
- `smi/PhEthOH.smi.json.lzma`: The files containing a list of QCSchema molecules converted from the SMILES.
- `smi/AlkEthOH_chain.smi.log`
- `smi/AlkEthOH_rings.smi.log`
- `smi/PhEthOH.smi.log`: A record of how many stereoisomers and conformers in the JSON file that were generated from the SMILES.
- `smi/PhEthOH/PhEthOH_mol2files_frosst-20200919T083119Z-001.zip`: The mol2 files containing the PhEthOh molecules.
- `00.preprocess_input_data.sh`: The commands needed to convert the extract and convert the mol2 files to create `PhEthOH.smi`
- `01.make_input_molecules.sh`: The commands needed to enumerate stereoisomers and generate conformers from `smi` files to `json` files.
- `02.submit_dataset.sh`: The commands needed to submit the JSON files in this submission to a QCArchive Optimization dataset.
- `00.preprocess_input_data.sh.log`
- `01.make_input_molecules.sh.log`
- `02.submit_dataset.sh.log`: The output logs of each of the three steps.
- `lib/process_PheEthOH_mol2_to_smi.sh`: The helper command to extract and convert the mol2 files to SMILES
- `lib/smi_to_json.sh`: The helper command to convert a file of SMILES to QCSchema.
- `lib/to_smi.sh`: Converts a mol2 using `tleap` and `obabel`.
- `conda-env.yml`: The conda environment to submit the calculations

 ### Metadata
```
{
	"submitter": "trevorgokey",
	"creation_date": "2020-09-18",
	"collection_type": "OptimizationDataset",
	"dataset_name": "OpenFF Sandbox CHO PhAlkEthOH v1.0", 
	"short_description": "A diverse set of CHO molecules",
	"long_description_url": "https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/2020-09-18-OpenFF-Sandbox-CHO-PhAlkEthOH",
	"long_description": "This dataset contains an expanded set of the AlkEthOH and PhEthOH datasets, which were used in the original derivation of the smirnoff99Frosst parameters.",
	"elements": [
		"C",
		"H",
		"O"
	],
	"change_log": [
		{"author": "trevorgokey",
		 "date": "2020-09-18",
		 "version": "1.0",
		 "description": "A diverse set of CHO molecules. The molecules in this set were generated to include all stereoisomers if chirality was ambiguous from the SMILES input. Conformations were generated which had an RMSD of at least 4 Angstroms from all other conformers"
		}
	]
}
```
