#!/usr/bin/env python

import os
import json
import tarfile

import qcportal as ptl
from qcelemental.models import Molecule

def read_submit_torsion_options(input_json):
    """ Read data from submit_torsion_options.json

    Returns
    -------
    selected_torsions: dict
        Dictionary for selected torsions, has this structure:
        {
            canonical_torsion_index1: {
                'initial_molecules': [ Molecule1a, Molecule1b, .. ],
                'atom_indices': [ (0,1,2,3) ],
                'attributes': {'canonical_explicit_hydrogen_smiles': .., 'canonical_isomeric_smiles': .., ..}
            },
            ..
        }
    """
    with open(input_json) as infile:
        submit_torsion_options = json.load(infile)
    selected_torsions = {}
    for data in submit_torsion_options:
        canonical_torsion_index = data['attributes']['canonical_torsion_label']
        selected_torsions[canonical_torsion_index] = {
            'initial_molecules': [data['initial_molecule']],
            'atom_indices': data['keywords']['dihedrals'],
            'attributes': data['attributes']['cmiles_id'],
        }
    return selected_torsions

print("Reading selected_torsions...")

if not os.path.exists('submit_torsion_options.json'):
    with tarfile.open('submit_torsion_options.json.tar.gz') as f:
        f.extractfile('submit_torsion_options.json')
selected_torsions = read_submit_torsion_options('submit_torsion_options.json')

print(f"Found {len(selected_torsions)} torsions")

print("Initializing dataset...")
# client = ptl.FractalClient.from_file()
client = ptl.FractalClient("localhost:7777", verify=False)

# create a new dataset with specified name
ds = ptl.collections.TorsionDriveDataset("OpenFF Group1 Torsions", client=client)

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
qc_spec = {"driver": "gradient", "method": "B3LYP-d3bj", "basis": "dzvp", "program": "psi4"}
ds.add_specification("default", opt_spec, qc_spec, description="Standard OpenFF torsiondrive specification.")

# add molecules
print(f"Adding {len(selected_torsions)} torsions")
i = 0
for canonical_torsion_index, torsion_data in selected_torsions.items():
    attributes = torsion_data['attributes']
    torsion_atom_indices = torsion_data['atom_indices']
    grid_spacings = [15] * len(torsion_atom_indices)
    initial_molecules = torsion_data['initial_molecules']
    print(i, canonical_torsion_index, len(initial_molecules))
    ds.add_entry(canonical_torsion_index, initial_molecules, torsion_atom_indices, grid_spacings, energy_upper_limit=0.05, attributes=attributes)
    i += 1
print("Submitting tasks...")
comp = ds.compute("default", tag="openff", priority="normal")
print(comp)

print("Complete!")
