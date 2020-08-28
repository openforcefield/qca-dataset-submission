# Pfizer discrepancy optimization dataset 1 (`OptimizationDataset`)

Prepared by John Chodera `<john.chodera@choderalab.org>` with molecules from Xinjun Hou `<Xinjun.Hou@pfizer.com>` and Brajesh Rai `<Brajesh.Rai@pfizer.com>`

This database is a subset of 100 challenging small molecule fragments where HF/minix followed by B3LYP/6-31G*//B3LYP/6-31G** differed substantially from OPLS3e.
Diffuse functions for sulfur-containing compounds and compounds with negative formal charge were added.

Molecules were submitted as-is with no expansion of stereochemistry, protonation, or tautomers.

Pfizer code for torsion drive exploration has already been made public as part of the manuscript (currently under revision).
https://github.com/PfizerRD/torsional-strain

### General Information

 - Date: 2019-09-07
 - Class: Forcefield Parametrization
 - Purpose: Explore discrepancies between QM and OPLS3e
 - Collection: OptimizationDataset
 - Name: Pfizer discrepancy optimization dataset 1
 - Number of Entries: 100 unique molecules, 352 conformers
 - Submitter: John Chodera

### Generation pipeline

 - Used `fragmenter.chemi.generate_conformers` to generate input conformers
 - Used `cmiles` to generate cmiles identifiers

### Quantum Chemistry Information

 - Theory: OpenFF high-throughput standard QC reference
 - Additional Properties: Computes and saves Wiberg bond-orders

### Notes

 - Some molecules could not be Kekulized: `Can't kekulize mol.`

### Manifest

 - `01_generate.py` - script to generate OptimizationDataset inputs
 - `02_create_optimization_dataset.py` - create the dataset
 - `03_visualize.py` - create a PDF of all species in the dataset
 - `PFE-OFF-100Frags.smi.gz` - Pfizer dataset of 100 SMILES strings
 - `optimization_inputs.json.gz` - QC Schema JSON input molecules for QCFractal
 - `optimization_inputs.smi` - SMILES input molecules for QCFractal (for inspection)
 - `optimization_inputs.pdf` - PDF 2D structures of input molecules for QCFractal (for inspection)
 - `optimization_inputs.sdf` - SDF file of generated conformers
 - `duplicates.smi` - rejected duplicate canonical SMILES strings
 - `cmiles_failures.smi` - Molecules that RDKit failed to generate standardized tautomer
 - `omega_failures.smi` - Molecules that Omega failed to generate conformers
 - `skipped_ions.smi` - Molecules that were skipped (ions)
 - `environment.yml` - versions of toolkits used (conda environment file)

### Usage

0. Generate SMILES files
```bash
python 00_extract_smiles.py
```
1. Generate the conformers and save in JSON format
```bash
python 01_generate.py
```
A new file called `optimization_inputs.json` will be generated.
2. Create dataset on QCFractal server
```bash
python 02_create_optimization_dataset.py optimization_inputs.json smirnoff_coverage client_config.yaml
```
3. Create PDF visualization
```bash
python 03_visualize.py
```
