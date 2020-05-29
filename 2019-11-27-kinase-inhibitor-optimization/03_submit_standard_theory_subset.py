#!/usr/bin/env python

import copy
import json
import tqdm
import tarfile
from collections import Counter
import qcportal as ptl
import qcelemental as qcel
from qcelemental.models import Molecule

dataset_name = "Kinase Inhibitors: WBO Distributions"
update = True

client = ptl.FractalClient.from_file()
print(client)

with open("subset_1.json", "r") as handle:
    subset = json.load(handle)

ds = client.get_collection("OptimizationDataset", dataset_name)

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

qc_spec = {
    "driver": "gradient",
    "method": "b3lyp-d3(bj)",
    "program": "psi4",
    "basis": "dzvp",
    "keywords": kw_id,
}  # Keywords 2 id map compute wiberg bond orders, dipoles, and quadrupoles.

# create specification for this dataset
ds.add_specification(
    "default",
    opt_spec,
    qc_spec,
    description="Standard OpenFF torsiondrive specification.",
)

print("Submitting tasks...")
comp = ds.compute("default", tag="openff", priority="low", subset=subset)
print(comp)

print("Complete!")
