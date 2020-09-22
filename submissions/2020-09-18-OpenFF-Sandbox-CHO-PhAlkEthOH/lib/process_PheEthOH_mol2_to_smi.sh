#!/bin/bash

# This file takes a mol2 and uses tleap to convert to PDB, then OpenBabel to 
# convert to SMILES. This is done because these mol2 files use AMBER types and 
# don't have any element names, so e.g. Os is parsed incorrectly as Osmium.
to_smi=$(realpath lib/to_smi.sh)

cd smi/PhEthOH

echo "Unzipping PhEthOH_mol2files_frosst-20200919T083119Z-001.zip"
unzip -q -f PhEthOH_mol2files_frosst-20200919T083119Z-001.zip
cd PhEthOH_mol2files_frosst

# gather the mol2 files into a sorted list, convert to SMILES using 16 jobs
xargs -a <(ls -1 *.mol2 | sort -k 1.9n)  -L 1 -P 16 bash $to_smi 

# Take the individual SMILES files and combine into a single dataset
for i in $(ls -1 *.smi | sort -k 1.9n ) ; do
	cat $i
done | tqdm --total=$(ls -1 *.smi | wc -l) --ncols=80 > ../../PhEthOH.smi 

echo "Wrote PhEthOH.smi"
