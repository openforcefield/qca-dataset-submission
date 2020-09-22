#!/bin/bash

mol2=$1
pdb=${1/mol2/pdb}
smi=${1/mol2/smi}

cat << EOF | tleap -f - &> /dev/null
x = loadmol2 $mol2
savepdb x $pdb
EOF

# obabel makes a lot of noise
obabel -i pdb $pdb -o smi -O $smi &> /dev/null

rm -f $pdb
