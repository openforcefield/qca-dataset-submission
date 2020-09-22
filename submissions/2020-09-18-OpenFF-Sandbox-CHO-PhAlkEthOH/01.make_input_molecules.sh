#!/bin/bash

log=$(basename $0).log
smi_dir="./smi"

# Converts the SMILES into a QCArchive compatible JSON format
echo | tee $log
date | tee -a $log

run() {
	local smi_dir=$1
	local smi=$2
	local log=$3
	echo "Converting ${smi} to JSON conformers. Writing results to ${smi_dir}/${smi}.log" | tee -a ${log}
	lib/smi_to_json.sh ${smi_dir} ${smi} | tee -a ${log}
	date | tee -a $log
	tail -n 4 ${smi_dir}/${smi}.log | tee -a ${log}
	echo | tee -a ${log}

	sed -i 's/^.*\r//' ${log}
}

run ${smi_dir} AlkEthOH_chain.smi ${log}
run ${smi_dir} AlkEthOH_rings.smi ${log}
run ${smi_dir} PhEthOH.smi ${log}
