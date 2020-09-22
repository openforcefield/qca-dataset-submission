#!/bin/bash

log=$(basename $0).log
for molecules in \
	AlkEthOH_chain.smi.json.lzma \
	AlkEthOH_rings.smi.json.lzma \
	PhEthOH.smi.json.lzma
do 
	date
	python3 -m offsb.ui.qca.submit-optimizations \
		--verbose \
		--server localhost:7777 \
		--compute-tag openff \
		--priority normal \
		--metadata metadata.json \
		--compute-spec compute.json \
		--input-molecules smi/${molecules}
done |& tee ${log}

# Since the progress bar updates use \r, strip intermediate progress bars such 
# that only the final progress bar is saved
sed -i 's/^.*\r//g' ${log}
