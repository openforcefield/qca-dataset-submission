#!/usr/bin/env python

import copy
import json
import tqdm
import tarfile
from collections import Counter
#import qcportal as ptl
import qcfractal.interface as ptl
from qcelemental.models import Molecule

# Updates the dataset rather than recomputing
UPDATE = False


print("Initializing dataset...")
client = ptl.FractalClient.from_file()
#client = ptl.FractalClient("localhost:7777", verify=False)

# create or pull a previous dataset with specified name

# Pull hte optimization dataset
opt_ds = client.get_collection("OptimizationDataset", "OpenFF Optimization Set 1")
opt_ds.query("default", force=True)

if UPDATE:
    hess_ds = client.get_collection("Dataset", "OpenFF Optimization Set 1")
else:
    hess_ds = ptl.collections.Dataset("OpenFF Optimization Set 1", client=client)
    hess_ds.data.default_program = "psi4"
    hess_ds.data.default_driver = "hessian"
    
    kw = ptl.models.KeywordSet(values={'maxiter': 200,
     'scf_properties': ['dipole',
      'quadrupole',
      'wiberg_lowdin_indices',
      'mayer_indices']})
    hess_ds.add_keywords("default", "psi4", kw, default=True)
    hess_ds.save()

# add molecprint(f"Adding {len(molecules_dict)} molecules")
print(f"Generating computation records")
known_jobs = set(hess_ds.df.index)

comp = []
for idx, opt in opt_ds.df["default"].iteritems():
    if opt.status == "INCOMPLETE":
        continue

    if idx in hess_ds.df.index:
        continue
        
    record = {"name": idx, "molecule_id": opt.final_molecule}
    comp.append(record)


print(f"Adding {len(comp)} computations")
hess_ds.data.records.extend([ptl.collections.dataset.MoleculeRecord(**x) for x in comp])

print("Submitting tasks...")
r = hess_ds.compute("B3LYP-d3bj", "DZVP", keywords="default", program="psi4", tag="openff", priority="normal")
print(r)

hess_ds.save()
#print(comp)

print("Complete!")
