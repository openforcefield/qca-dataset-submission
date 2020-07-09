#!/usr/bin/env python

import copy
import json
import tqdm
import tarfile
from collections import Counter
import qcportal as ptl
from qcelemental.models import Molecule


def read_molecules(input_json):
    """ Extract the molecules and the index of them from the input json file

    Parameters
    ----------
    input_json: str,
        JSON file name to the output json of generate.py
        The data in the json file should be a list of {'initial_molecules': [..], 'cmiles_identifiers':{}}.

    Returns
    -------
    molecules_dict: dict
        The dictionary maps the index of a molecule to a Molecule object. e.g.
        {
            index1: Molecule1,
            index2: Molecule2,
        }

    molecule_attributes: dict
        The dicitonary maps the index of a molecule to the attributes of the molecule, e.g.
        {
            index1: {'canonical_explicit_hydrogen_smiles': .., 'canonical_isomeric_smiles': .., ..}
        }

    Note
    ----
    1. The mdata['cmiles_identifiers']['canonical_isomeric_smiles'] is selected as the index.
    2. For molecules have the same "canonical_isomeric_smiles", we use index-1, index-2 to distinguish them.
    """
    molecules_dict = {}
    molecule_attributes = {}

    if input_json.endswith(".tar") or input_json.endswith(".tar.gz"):

        extract_file = input_json.replace(".gz", "").replace(".tar", ".json")
        with tarfile.open(input_json, 'r') as infile:

            molecule_data_list = json.load(infile.extractfile(extract_file))

    if input_json.endswith(".gz"):

        import gzip
        with gzip.open(input_json, 'r') as infile:

            molecule_data_list = json.loads(infile.read().decode('utf-8'))

    else:
        with open(input_json) as infile:
            molecule_data_list = json.load(infile)

    index_counter = Counter()
    for mdata in molecule_data_list:
        initial_molecules = mdata['initial_molecules']
        cmiles_ids = mdata['cmiles_identifiers']
        index = cmiles_ids['canonical_isomeric_smiles']
        for i_conformer, initial_molecule in enumerate(initial_molecules):
            qcel_molecule = Molecule.from_data(initial_molecule)
            # use count to generate unique index
            index_count = index_counter[index]
            this_index = f'{index}-{index_count}'
            index_counter[index] += 1
            assert this_index not in molecules_dict, f"Multiple molecules have the same index, please check {mdata}"
            molecules_dict[this_index] = qcel_molecule
            molecule_attributes[this_index] = cmiles_ids
    return molecules_dict, molecule_attributes


print("Extracting molecules...")
molecules_dict, molecule_attributes = read_molecules("optimization_inputs.json.gz")

print("Initializing dataset...")
#client = ptl.FractalClient("localhost:7777", verify=False) # TODO: Should this be changed to remote address?
client = ptl.FractalClient().from_file()

# create a new dataset with specified name
ds = ptl.collections.OptimizationDataset("Pfizer Discrepancy Optimization Dataset 1", client=client)

kw = ptl.models.KeywordSet(values={'maxiter': 200,
 'scf_properties': ['dipole',
  'quadrupole',
  'wiberg_lowdin_indices',
  'mayer_indices']})
kw_id = client.add_keywords([kw])[0]

# create specification for this dataset
opt_spec = {"program": "geometric"}
qc_spec = {"driver": "gradient", "method": "B3LYP-d3bj", "basis": "dzvp", "program": "psi4", "keywords": kw_id}
ds.add_specification("default", opt_spec, qc_spec, description="Standard OpenFF optimization quantum chemistry specification.")

# add molecules
print(f"Adding {len(molecules_dict)} molecules")
for molecule_index, molecule in tqdm.tqdm(molecules_dict.items()):
    attributes = molecule_attributes[molecule_index]
    ds.add_entry(molecule_index, molecule, attributes=attributes)

print("Submitting tasks...")
comp = ds.compute("default", tag="openff", priority="normal")
print(comp)

print("Complete!")
