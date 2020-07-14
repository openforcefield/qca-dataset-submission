# 2-D Nitrogen Grid Optimizations

### General Information
 - Date: 12/11/2019
 - Class: Forcefield Parameterization
 - Purpose: Set of diverse trivalent nitrogen molecules for 2-D grid optimization
 - Collection: GridOptimization
 - Name: Nitrogen Grid Optimization
 - Number of Entries: 323 2-D Grid optimizations
 - Submitter: Jessica Maat

### Generation
Details regarding the generation of this set of molecule is outline in the DANCE repository by Bryon Tjanaka - https://github.com/btjanaka/dance


### Scripts in File
 - `Molecules_to_run.tar.gz` - A folder that contains the .sdf files of the trivalent nitrogen molecules for the grid optimization calculations
 - `sendmols_to_server.py` - Contains script to set up the json input for the 2-d grid optimization for a given coordinate file, such as an `*.sdf`
 - `requirements.txt` - Contains the installments and versions used
 - `nitrogen_Jobs_2dscans.json` - Contains the qcschema for the molecules to run

### Notes
 - `sendmols_to_servers.py` should be modified to read in the .sdf files in `Molecules_to_run`

