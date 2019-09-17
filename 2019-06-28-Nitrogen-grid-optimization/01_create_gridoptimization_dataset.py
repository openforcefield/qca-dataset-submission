#!/usr/bin/env python

import os
import json
import tarfile

#import qcportal as ptl
import qcfractal.interface as ptl
import qcelemental as qcel

dataset_name = "OpenFF Trivalent Nitrogen Set 1" 
input_json = 'gridoptimization_inputs.json'

with open(input_json, "r") as handle:
    inputs = json.load(handle)

print(f"Found {len(inputs)} gridoptimizations")

print("Initializing dataset...")
client = ptl.FractalClient.from_file()
#client = ptl.FractalClient("localhost:7777", verify=False)

# create a new dataset with specified name
ds = ptl.collections.GridOptimizationDataset(dataset_name, client=client)

# create specification for this dataset
opt_spec = {
    "program": "geometric",
    "keywords": {
        "coordsys": "tric",
        "enforce": 0.1,
        "reset": True,
        "qccnv": True,
        "epsilon": 0.0,
    }
}

kw = ptl.models.KeywordSet(values={'maxiter': 200,
 'scf_properties': ['dipole',
  'quadrupole',
  'wiberg_lowdin_indices',
  'mayer_indices']})
kw_id = client.add_keywords([kw])[0]

#qc_spec = {"driver": "gradient", "method": "uff", "program": "rdkit"}
qc_spec = {"driver": "gradient", "method": "B3LYP-d3bj", "basis": "dzvp", "program": "psi4", "keywords": kw_id}
ds.add_specification("default", opt_spec, qc_spec, description="Standard OpenFF gridoptimization specification.")

#qc_spec = {"driver": "gradient", "method": "MP2", "basis": "def2-SV(P)", "program": "psi4", "keywords": kw_id}
#ds.add_specification("mp2", opt_spec, qc_spec, description="MP2/def2-SV(P) .")

# add molecules
print(f"Adding {len(inputs)} torsions")
i = 0
for index, data in inputs.items():
    attributes = data['cmiles_ids']
    scans = data['keywords']['scans']
    initial_molecule = data['initial_molecule']
    preopt = data['keywords']['preoptimization']
    ds.add_entry(index, initial_molecule, scans, preoptimization=preopt, attributes=attributes)
    i += 1

print("Submitting tasks...")
comp = ds.compute("default", tag="openff")
print(comp)

print("Complete!")
