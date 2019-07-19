#!/usr/bin/env python

"""
Extract SMILES from SMIRNOFF_sub3rot.tar.gz to generate input.smi
"""

from openeye import oechem
from openeye import oemolprop
import gzip

filterfile = oechem.oeifstream('oechem-filterfile')
filter = oemolprop.OEFilter(filterfile)

ifs = oechem.oemolistream('nci-250k.smi.gz')
ofs = oechem.oemolostream('input.smi')
for mol in ifs.GetOEMols():
    if filter(mol):
        smiles = oechem.OEMolToSmiles(mol)
        oechem.OEWriteMolecule(ofs, mol)

ifs.close()
ofs.close()
