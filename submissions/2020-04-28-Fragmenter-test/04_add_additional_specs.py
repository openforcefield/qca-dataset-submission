#!/usr/bin/env python

import os
import qcportal as ptl

#collection_name = "Fragmenter paper"
collection_name = "OpenFF Fragmenter Validation 1.0"
UPDATE = True
local_run = False

print("Connecting to server...")
if local_run:
    client = ptl.FractalClient("localhost:7777", verify=False)
else:
    client = ptl.FractalClient.from_file()

if UPDATE:
    print("Retreiving dataset:", collection_name)
    ds = client.get_collection("TorsionDriveDataset", collection_name)

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

    qc_spec = {
        "driver": "gradient",
        "method": "",
        "basis": "",
        "program": "psi4",
        "keywords": 2
    }  # Keywords 2 id map compute wiberg bond orders, dipoles, and quadrupoles.

    methods = [
            "B3LYP-d3bj",
            "BLYP-d3bj",
            "LRC-WPBEH",
            "MN15-d3bJ",
            "TPSSH-d3bJ"
    ]
    bases   = [
        "dzvp",
        "aug-cc-pvdz",
        "aug-cc-pvtz",
        "def2-tzvppd",
    ]
    for method in methods:
        for basis in bases:
            qc_spec['method'] = method
            qc_spec['basis']  = basis
            model = "/".join((method,basis))
            # already the "default" spec, so skip
            if (method.lower() == "b3lyp-d3bj") and (basis.lower() == "dzvp"):
                continue
            print("Submitting model:", model)
            ds.add_specification(
                model,
                opt_spec,
                qc_spec,
                description="OpenFF model {:s} TorsionDrive exploration of charged molecules".format(model)
            )
            ds.save()
            ds.compute(model, tag="openff", priority="LOW")

    print("Complete!")
