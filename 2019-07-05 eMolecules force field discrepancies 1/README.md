# eMolecules force field discrepancy optimization dataset 1

Prepared by John Chodera `<john.chodera@choderalab.org>` and David Mobley `<dmobley@uci.edu>` based on work done by Jordan Ehrman `<ehrmanj@uci.edu>` supervised by Mobley, Caitlin Bannan and Victoria Lim (UCI).

### Provenance Information

The Mobley lab has been conducting a project to energy minimize all of eMolecules under the GAFF, GAFF2, smirnoff99Frosst and MMFF94 and 94s force fields and compare geometries, looking for significant differences. Results are compared by RMSD and torsion fingerprint difference (code being made available at https://github.com/MobleyLab/off-ffcompare and a paper is in preparation). Here, we have selected a subset of the resulting minimized molecules based on picking those with geometries that are substantially different in `smirnoff99Frosst` relative to the other force fields. We further filtered the results based on picking molecules with three or fewer rotatable bonds, since at present we are not planning on fragmenting the molecules before running geometry optimizations. The molecule set was selected from the minimized geometries and RMSD/TFD results using Ehrman's [dataproliferator](https://github.com/MobleyLab/off-ffcompare/blob/721b4c8a67276c381314a2287b4c1cfe76356f01/dataproliferator.py) script. Subsequent filtering by rotatable bonds and visualization was done in the attached Jupyter notebook, `filter_by_bonds.ipynb`, resulting in the final set `SMIRNOFF_sub3rot.tar.gz` provided here.

Following this, the set was further prepared by expanding undefined stereochemistry (but not protonation or tautomeric states) and submitted for optimization and Hessian computations.

### General Information

 - Date: 2019-07-05
 - Class: Forcefield Parametrization
 - Purpose: Gather quantum chemical valence data on small molecules that show discrepancies between force fields
 - Collection: OptimizationDataset
 - Name: OpenFF Discrepancy Benchmark 1
 - Number of Entries: 2802 unique molecules, 19712 total conformers
 - Submitter: John Chodera

### Generation pipeline

 - Used `fragmenter.fragment.expand_states` to expand stereoisomers
 - Used `fragmenter.chemi.generate_conformers` to generate input conformers
 - Used `cmiles` to generate cmiles identifiers

### Quantum Chemistry Information

 - Theory: OpenFF high-throughput standard QC reference
 - Additional Properties: Computes and saves Wiberg bond-orders

### Notes

 - Some molecules could not be Kekulized: `Can't kekulize mol.`

### Manifest

 - `00_extract_smiles.py` - script to extract SMILES from `SMIRNOFF_sub3rot.tar.gz` and generate `input.smi`
 - `01_generate.py` - script to generate OptimizationDataset inputs
 - `02_create_optimization_dataset.py` - create the dataset
 - `03_visualize.py` - create a PDF of all species in the dataset
 - `SMIRNOFF_sub3rot.tar.gz` - original file containing molecules of interest
 - `input.smi` - initial SMILES
 - `optimization_inputs.json.gz` - QC Schema JSON input molecules for QCFractal
 - `optimization_inputs.smi` - SMILES input molecules for QCFractal (for inspection)
 - `optimization_inputs.pdf` - PDF 2D structures of input molecules for QCFractal (for inspection)
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
