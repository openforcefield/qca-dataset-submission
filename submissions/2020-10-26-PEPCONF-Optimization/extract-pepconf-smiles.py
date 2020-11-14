#!/bin/env python

"""
Extract SMILES strings
"""

datasets = [
    'dipeptide',
    'tripeptide',
    'disulfide',
    'bioactive',
    'cyclic',
]

def xyz_to_smiles(filename):
    """
    Parse an .xyz file, heuristically perceive chemistry, and return canonical isomeric SMILES
    """
    from openeye import oechem
    ifs = oechem.oemolistream(filename)
    for oemol in ifs.GetOEGraphMols():
        return oechem.OEMolToSmiles(oemol)

all_smiles = dict()
import glob, os, re
for dataset in datasets:
    filenames = glob.glob(os.path.join('pepconf', dataset, 'xyz', '*.xyz'))
    for filename in filenames:
        head, tail = os.path.split(filename)
        root, ext = os.path.splitext(tail)
        match = re.match('^(?P<name>\S+)_(\d+)$', root)
        name = match.group('name')
        smiles = xyz_to_smiles(filename)
        all_smiles[smiles] = name

# Write CSV
print('Writing...')
with open('pepconf.csv', 'w') as outfile:
    for smiles, name in all_smiles.items():
        outfile.write(f'{smiles},{name}\n')
