# Heteroaromatic rings of the future (optimization dataset)

Prepared by John Chodera `<john.chodera@choderalab.org>`

VEHICLe (virtual exploratory heterocyclic library) dataset of 24,847 aromatic heterocyclic rings from [1].

Prepared by enumerating protonation and tautomeric states and submitted for optimization and Hessian computations.

### General Information
 - Date: 2019-07-02
 - Class: Forcefield Parametrization
 - Purpose: VEHICLe (virtual exploratory heterocyclic library)
 - Collection: OptimizationDataset
 - Name: Open Force Field VEHICLe 1.0 optimization dataset
 - Number of Entries: 84761
 - Submitter: John Chodera

### Generation pipeline
 - Used `fragmenter.fragment.expand_states` to expand protonation states and stereoisomers
 - Used `fragmenter.chemi.generate_conformers` to generate input conformers
 - Used `cmiles` to generate cmiles identifiers
 - SMILES downloaded from ftp://ftp.ebi.ac.uk/pub/databases/chembl/VEHICLe/

### Quantum Chemistry Information
 - Theory: OpenFF high-throughput standard QC reference
 - Additional Properties: Computes and saves Wiberg bond-orders

### Notes
 - Some molecules could not be Kekulized: `Can't kekulize mol.`

### Manifest
 - `01_generate.py` - script to generate OptimizationDataset inputs
 - `02_create_optimization_dataset.py` - create the dataset
 - `03_visualize.py` - create a PDF of all species in the dataset
 - `VEHICLe.csv` - downloaded VEHICLe database containing SMILES
 - `VEHICLe.smi` - input VEHICLe database containing only SMILES
 - `optimization_inputs.json.gz` - input molecules for QCFractal
 - `cmiles_failures.smi` - Molecules that RDKit failed to generate standardized tautomer
 - `omega_failures.smi` - Molecules that Omega failed to generate conformers
 - `skipped_ions.smi` - Molecules that were skipped (ions)
 - `environment.yml` - versions of toolkits used (conda environment file)

### Usage
1. Generate the conformers and save in JSON format
   ```bash
   python generate.py
   ```
   A new file called `optimization_inputs.json` will be generated.

2. Create dataset on QCFractal server
    ```bash
    python create_optimization_dataset.py optimization_inputs.json smirnoff_coverage client_config.yaml
    ```

## References

[1] William R. Pitt, David M. Parry, Benjamin G. Perry, and Colin R. Groom.
Heteroaromatic rings of the future.
J Med Chem 52:2952, 2009.
[[DOI]](http://doi.org/10.1021/jm801513z) [[FTP]](ftp://ftp.ebi.ac.uk/pub/databases/chembl/VEHICLe/)
