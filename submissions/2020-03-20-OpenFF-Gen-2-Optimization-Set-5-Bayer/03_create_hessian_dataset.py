#!/usr/bin/env python

import copy
import json
import tqdm
import tarfile
import pandas as pd
from collections import Counter
#import qcportal as ptl
import qcfractal.interface as ptl
from qcelemental.models import Molecule

# Updates the dataset rather than recomputing
UPDATE = True
name = "OpenFF Gen 2 Opt Set 5 Bayer"

print("Initializing dataset...")
client = ptl.FractalClient.from_file()
#client = ptl.FractalClient("localhost:7777", verify=False)

# create or pull a previous dataset with specified name

# Pull the optimization dataset
opt_ds = client.get_collection("OptimizationDataset", name)
opt_ds.query("default", force=True)

if UPDATE:
    hess_ds = client.get_collection("Dataset", name)
else:
    hess_ds = ptl.collections.Dataset(name, client=client, default_program="psi4", default_driver="hessian")

    kw = ptl.models.KeywordSet(values={'maxiter': 200,
     'scf_properties': ['dipole',
      'quadrupole',
      'wiberg_lowdin_indices',
      'mayer_indices']})
    hess_ds.add_keywords("default", "psi4", kw, default=True)
    hess_ds.save()

print(f"Generating computation records")

molecule_ids = []
molecule_idx = []
for idx, opt in opt_ds.df["default"].iteritems():

    if pd.isnull(opt):
        continue

    if opt.status != "COMPLETE":
        continue

    if idx in hess_ds.df.index:
        continue

    molecule_ids.append(opt.final_molecule)
    molecule_idx.append(idx)

print(f'Querying info for {len(molecule_ids)} molecules')

# Query/submit 250 at a time
for i in range(0, len(molecule_ids), 250):
    query_ids = molecule_ids[i:i+250]
    query_idx = molecule_idx[i:i+250]
    molecule_objs = client.query_molecules(id=query_ids)
    assert len(molecule_objs) == len(query_ids)
    assert len(molecule_objs) == len(query_idx)

    for idx,mobj in zip(query_idx, molecule_objs):
        hess_ds.add_entry(idx, mobj)

    hess_ds.save()

print("Submitting tasks...")
r = hess_ds.compute("B3LYP-d3bj", "DZVP", keywords="default", program="psi4", tag="openff", priority="high")
print(r)

hess_ds.save()

print("Complete!")
