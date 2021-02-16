# Nitrogen Grid Optimizations

### General Information
 - Date: Jan 16, 2020
 - Class: Forcefield Parameterization
 - Purpose: Set of diverse trivalent nitrogen molecules for 1-D grid optimization, this is a secondary dataset
 - Collection: GridOptimization
 - Name: OpenFF Trivalent Nitrogen Set 3
 - Number of Entries: 116 1-D Grid optimizations
 - Submitter: Jessica Maat

### Generation
These molecules were selected by Jessica Maat to sample the chemical enviornment around a trivalent nitrogen and look at combinatorial effects of various functional groups on the pyramidalization of a trivalent nitrogen center.

### Scripts in File
 - `sendmols_to_server.py` - Contains script to set up the json input for the 1-d grid optimization for a given coordinate file, such as an `*.sdf`
 - `requirements.txt` - Contains the installments and versions used
 - `1d_scan_jobs.json` - Contains the QC schema

