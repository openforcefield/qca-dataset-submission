# smirnoff99Frost parameter coverage

### General Information
 - Date: 6/25/2019
 - Class: smirnoff99Frost parameter coverage
 - Purpose: Set of small molecules that use all smirnoff99Frost parameters
 - Collection: OptimizationDataset
 - Name: 
 - Number of Entries: 156 (expandable)
 - Submitter: Chaya Stern

### Generation pipeline
 - Used `fragmenter.fragment.expand_states`` to expand protonation states and stereoisomers
 - Used `fragmenter.chemi.generate_conformers`` to generate input conformers
 - Used `cmiles`` to generate cmiles identifiers
 - Original SMILES were taken from [openforcefiles/open-forcefield-data/Utilize-All_Parameters](https://github.com/openforcefield/open-forcefield-data/tree/master/Utilize-All-Parameters/selected)

### Quantum Chemistry Information
 - Theory: OpenFF high-throughput standard QC reference
 - Additional Properties: Computes and saves Wiberg bond-orders
 
### Notes
 - Removed `CSSCCN=C=S` from `chosen.smi` because it causes `oequacpac.OEEnumerateFormalCharges()` 
   to segfault. 
 - Filtered out ions. This set includes molecules that are missing both valence and nonbonded parameters. 
 - Omega failed to generate conformers for some molecules. The errors were all: `Warning: force field setup failed due to missing parameters for molecule` 
 - RDkit failed to generate standard smiles for some molecules. All failures were due to wrong Explicit valence. 
  
### Manifest
 - `chosen.smi` - input SMILES file. 
 - `generate.py` - script to generate OptimizationDataset inputs
 - `optimization_inputs.json` - input molecules
 - `cmiles_failures.json` - Molecules that rdkit failed to generate standardized tautomer
 - `omega_failures.json` - Molecules that Omega failed to generate conformers
 - `skipped_ions.json` - Molecules that were skipped (ions)
 - `requirements.txt` - versions of toolkits used
 

