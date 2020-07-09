#!/usr/bin/env python

import os
import json
import tarfile
import tqdm

import qcportal as ptl

# import qcfractal.interface as ptl
import qcelemental as qcel

dataset_name = "OpenFF Trivalent Nitrogen Set 2"
input_json = "nitrogen_Jobs_2dscans.json"
UPDATE = True

with open(input_json, "r") as handle:
    inputs = json.load(handle)

print(f"Found {len(inputs)} gridoptimizations")

print("Initializing dataset...")
client = ptl.FractalClient.from_file()
# client = ptl.FractalClient("localhost:7777", verify=False)

# create a new dataset with specified name
if UPDATE:
    print("Pulling dataset from server...")
    ds = client.get_collection("GridOptimizationDataset", dataset_name)
else:
    ds = ptl.collections.GridOptimizationDataset(dataset_name, client=client)
    print("Creating new dataset...")

    # create specification for this dataset
    opt_spec = {
        "program": "geometric",
        "keywords": {
            "coordsys": "tric",
            "enforce": 0.1,
            "reset": True,
            "qccnv": True,
            "epsilon": 0.0,
        },
    }

    kw = ptl.models.KeywordSet(
        values={
            "maxiter": 200,
            "scf_properties": [
                "dipole",
                "quadrupole",
                "wiberg_lowdin_indices",
                "mayer_indices",
            ],
        }
    )
    kw_id = client.add_keywords([kw])[0]

    # qc_spec = {"driver": "gradient", "method": "uff", "program": "rdkit"}
    qc_spec = {
        "driver": "gradient",
        "method": "B3LYP-d3bj",
        "basis": "dzvp",
        "program": "psi4",
        "keywords": kw_id,
    }
    ds.add_specification(
        "default",
        opt_spec,
        qc_spec,
        description="Standard OpenFF gridoptimization specification.",
    )

# qc_spec = {"driver": "gradient", "method": "MP2", "basis": "def2-SV(P)", "program": "psi4", "keywords": kw_id}
# ds.add_specification("mp2", opt_spec, qc_spec, description="MP2/def2-SV(P) .")

inputs = list(inputs.items())

# add molecules
print(f"Adding {len(inputs)} gridoptimizations")
i = 0
for index, data in tqdm.tqdm(inputs):
    attributes = data["cmiles_ids"]
    scans = data["keywords"]["scans"]
    initial_molecule = data["initial_molecule"]
    preopt = data["keywords"]["preoptimization"]

    # Check to make sure their is a decent amount of connectivity to watch bohr/angstrom issues
    conn = qcel.molutil.guess_connectivity(
        initial_molecule["symbols"], initial_molecule["geometry"]
    )
    assert (len(conn) + 3) > len(initial_molecule["symbols"]), conn

    try:
        ds.add_entry(
            index,
            initial_molecule,
            scans,
            preoptimization=preopt,
            attributes=attributes,
        )
    except:
        continue
    i += 1

print("Submitting tasks...")
comp = ds.compute("default", tag="openff")
print(comp)

print("Complete!")
