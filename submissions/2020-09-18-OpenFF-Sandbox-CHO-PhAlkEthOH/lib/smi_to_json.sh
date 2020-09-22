#!/bin/bash
smi=$2
cd $1
python3 -u -m offsb.ui.smiles.load $smi -c 3 -n 100 -o ${smi}.log -f ${smi}.json --json | sed '/No principle axes found during inertial alignment/d'
lzma -9v --force ${smi}.json
