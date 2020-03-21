# OpenFF Gen 2 Opt Set 5 Bayer

### General Information
 - Date: 3/20/2020
 - Class: Force Field Parameterization
 - Purpose: 2nd generation optimization dataset for bond and valence parameter fitting
 - Collection: OptimizationDataset
 - Name: OpenFF Gen 2 Opt Set 5 Bayer
 - Number of Entries: - 16424 optimizations
 - Number of Molecules: - 524 molecules
 - Submitter: Jessica Maat (jmaat@uci.edu)

### Generation procedure

1. load_DS : Loads the ds from .txt file or from QCA DS
2. paramUsage: Creates dictionary of parameters (keys) and molecules which use these parameters
3. selectDiverseMols: Returns a list of diverse molecules after clustering with DBSCAN based off graph similarity scores
4. create_QCSChema: Expands tautomeric and isomeric state using cmiles and fragmenter and creates .json file from final list of molecules for QCA input
5. Submit the optimizations to QCFractal server

### Manifest

Cluster/filtering:
 - `generateOptDS.py`: Generates a .json of filtered optimization inputs from QCA data set
 - `openff_unconstrained-1.0.0-RC1.offxml`: .offxml file that contain parameters molecules are assigned to

Final .json:
 - `optimization_inputs.json.gz`: .json file that contains optimization inputs

Logs for filtering/clustering process:
 - `anglebond.p`: Dictionary of the parameters as keys and smiles that utilize the parameter key.
 - `optimization_inputs.sdf`: sdf file of the molecules that are optimization inputs
 - `duplicates.smi`: Duplicate smiles from makeJson function
 - `cmiles_failures.smi`: Molecules that failed in cmiles generation
 - `omega_failures.smi`: Molecules that failed during omega generation
 - `skipped_ions.smi`: List of the ions that were skipped as optimization inputs
 - `optimization_inputs.sdf`: .sdf file of the optimization inputs, including conformers
 - `optimization_inputs.smi`: Smiles inputs of the optimization inputs
 - `bond.txt`: A list of the bond parameters that are covered in the data set
 - `angleCoverage.txt`: A list of the angle parameters that are covered in the data set
 - `finalCounts.txt`: Gives the total molecules in the final data set and the total number of conformers for optimization

Requirements:
 - `requirements.txt`: Packages and versions used

### Usage
1. Run generateOptDS.py with a string containing the OptimizationDataset name.
    ```
    python 01_generateOptDS.py '20191216_BayerPatentSpace_forOpenFF_SmilesOnly.sdf'
    ```
2. Submit the .json file that is created in step 1 to QCFractal server.
   ```
   python 02_create_optimization_dataset.py 'OpenFF Gen 2 Opt Set 5 Bayer'
   ```


### Contributors:
 - Jessica Maat, Hyesu Jang, David Mobley, Lee-Ping Wang, Jeffrey Wagner, Josh Horton, Chaya Stern

