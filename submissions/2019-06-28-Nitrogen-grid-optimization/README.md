# Nitrogen Grid Optimizations

### General Information
 - Date: 7/2/2019
 - Class: Forcefield Parameterization
 - Purpose: Set of diverse trivalent nitrogen molecules for 1-D grid optimization
 - Collection: GridOptimization
 - Name: OpenFF Trivalent Nitrogen Set 1
 - Number of Entries: 323 1-D Grid optimizations
 - Submitter: Jessica Maat

### Generation
Details regarding the generation of this set of molecule is outline in the DANCE repository by Bryon Tjanaka - https://github.com/btjanaka/dance


### Scripts in File
 - `Molecules_to_run.tar.gz` - A folder that contains the .sdf files of the trivalent nitrogen molecules for the grid optimization calculations
 - `sendmols_to_server.py` - Contains script to set up the json input for the 1-d grid optimization for a given coordinate file, such as an `*.sdf`
 - `requirements.txt` - Contains the installments and versions used

### Notes
 - `sendmols_to_servers.py` should be modified to read in the .sdf files in `Molecules_to_run`

