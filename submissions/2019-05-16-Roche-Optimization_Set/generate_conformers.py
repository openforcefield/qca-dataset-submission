"""
This script skips the fragmentation step to find torsions and generate multiple conformations of molecules
"""

from fragmenter import workflow_api, chemi
from qcfractal import interface as ptl
import cmiles
from openeye import oechem
import json


client = ptl.FractalClient('https://localhost:7777/', verify=False)
fragmenter_wf = workflow_api.WorkFlow(client=client, workflow_id='torsiondrive_input', workflow_json='example_workflow.json')


# Load Roche molecules
roche_mols = chemi.file_to_oemols('OpenFF_references.sdf')
smiles = [cmiles.utils.mol_to_smiles(mol, mapped=False, explicit_hydrogen=False) for mol in roche_mols]

# Put smiles in format for generating torsion input files
frags = {}
for sm in smiles:
    identifiers = cmiles.get_molecule_ids(sm, toolkit='openeye', strict=False)
    frags[sm] = {'identifiers':identifiers, 'provenance': {'routine': {}}}

# Generate torsiondrive input
torsiondrive_inputs = {}
for frag in frags:
    td_input = fragmenter_wf.generate_torsiondrive_input(frags[frag])
    torsiondrive_inputs.update(td_input)

# Save
with open('torsiondrive_input.json', 'w') as f:
    json.dump(torsiondrive_inputs, f, indent=2, sort_keys=True)

# Added comment 2019-09-04
# Since this script was written using an older verison of fragmenter before there was a function to generate optimizaiton
# inputs, some manipulation of the data was needed. In addition, molecules that do not have torsions needed conformers


unrestrained = {}
no_conformers = []
for frag in torsiondrive_inputs:
    jobs = list(torsiondrive_inputs[frag]['torsiondrive_input'].keys())
    if not jobs:
        print(frag)
        unrestrained[frag] = {}
        no_conformers.append(frag)
    else:
        # First remove ids from initial molecules
        for mol in torsiondrive_inputs[frag]['torsiondrive_input'][jobs[0]]['initial_molecule']:
            ids = mol.pop('identifiers')
        unrestrained[frag] = {'cmiles_ids': ids,
                             'initial_molecules': torsiondrive_inputs[frag]['torsiondrive_input'][jobs[0]]['initial_molecule']}

# Add conformers manually
for smile in no_conformers:
    mol_id = cmiles.get_molecule_ids(smile, strict=False)
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smile)
    conformers = chemi.generate_conformers(mol)
    initial_molecules = [cmiles.utils.mol_to_map_ordered_qcschema(conf, mol_id) for conf in conformers.GetConfs()]
    print(len(initial_molecules))
    for m in initial_molecules:
        m.pop('identifiers')
    unrestrained[smile] = {'cmiles_ids': mol_id,
                          'initial_molecules': initial_molecules}

with open('geom_opt.json', 'w') as f:
    json.dump(unrestrained, f, indent=2, sort_keys=True)
