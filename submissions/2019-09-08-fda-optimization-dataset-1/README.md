# FDA optimization dataset 1 (`OptimizationDataset`)

Prepared by John Chodera `<john.chodera@choderalab.org>`

The [ZINC15 FDA dataset](http://zinc.docking.org/substances/subsets/fda/) was retrieve in `mol2` format on Sun Sep  8 20:44:34 EDT 2019 via:
http://zinc.docking.org/substances/subsets/fda.mol2?count=all

Some molecules appear in multiple protonation/tautomeric states.
The following [oechem filterfile](https://docs.eyesopen.com/toolkits/python/molproptk/filter_files.html) was applied:
```
MIN_NUM_HVY 3 "Minimum number of heavy atoms"
MAX_NUM_HVY 55 "Maximum number of heavy atoms"
MIN_ROT_BONDS 0 "Minimum number of rotatable bonds"
MAX_ROT_BONDS 5 "Maximum number of rotatable bonds"
ADJUST_ROT_FOR_RING true "BOOLEAN for whether to estimate degrees of freedom in rings"
ELIMINATE_METALS Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd
ALLOWED_ELEMENTS H,C,N,O,F,P,S,Cl,Br,I,B
MAX_RING_SIZE 20 "Maximum atoms in any ring system"
```
Up to 20 conformers/molecule were enumerated.

### General Information

 - Date: 2019-09-08
 - Class: Forcefield Parametrization
 - Purpose: Assess quality of conformations and energetics for FDA-approved inhibitors
 - Collection: OptimizationDataset
 - Name: FDA optimization dataset 1
 - Number of Entries: 1038 unique molecules, 6675 conformers
 - Submitter: John Chodera

### Generation pipeline

 - Filtered according to `oechem-filterfile`
 - Used `fragmenter.chemi.generate_conformers` to generate input conformers (max_confs=20)
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
 - `fda.mol2.gz` - FDA molecules retrieved from ZINC15 in mol2 format
 - `optimization_inputs.json.gz` - QC Schema JSON input molecules for QCFractal
 - `optimization_inputs.smi.gz` - SMILES input molecules for QCFractal (for inspection)
 - `optimization_inputs.pdf` - PDF 2D structures of input molecules for QCFractal (for inspection)
 - `optimization_inputs.sdf.gz` - SDF file of generated conformers
 - `duplicates.smi.gz` - rejected duplicate canonical SMILES strings
 - `cmiles_failures.smi.gz` - Molecules that RDKit failed to generate standardized tautomer
 - `omega_failures.smi.gz` - Molecules that Omega failed to generate conformers
 - `skipped.smi.gz` - Molecules that were skipped or filtered out
 - `environment.yml` - versions of toolkits used (conda environment file)

### Usage

1. Generate the conformers and save in JSON format
```bash
python 01_generate.py
```
A new file called `optimization_inputs.json` will be generated.
2. Create dataset on QCFractal server
```bash
python 02_create_optimization_dataset.py
```
3. Create PDF visualization
```bash
python 03_visualize.py
```
