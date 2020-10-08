#!/bin/bash

# This command requires tleap and obabel

log=$(basename $0).log

lib/process_PheEthOH_mol2_to_smi.sh |& tee ${log}
sed -i 's/^.*\r//' ${log}
