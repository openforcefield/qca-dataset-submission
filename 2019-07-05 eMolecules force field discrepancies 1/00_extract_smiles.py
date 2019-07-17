#!/usr/bin/env python

"""
Extract SMILES from SMIRNOFF_sub3rot.tar.gz to generate input.smi
"""

import os
import shutil
from openeye import oechem
import tarfile

oemol = oechem.OEMol()
with open('input.smi', 'w') as outfile:
    with tarfile.open(name='SMIRNOFF_sub3rot.tar.gz', mode='r:gz') as tar:
        for tarinfo in tar:
            if '.mol2' in tarinfo.name:
                print(tarinfo.name)
                # Extract the file
                tar.extract(tarinfo.name, path='extract')
                # Read mol2
                ifs = oechem.oemolistream(os.path.join('extract', tarinfo.name))
                oechem.OEReadMolecule(ifs, oemol)
                ifs.close()
                # Generate SMILES
                smiles = oechem.OEMolToSmiles(oemol)
                outfile.write(smiles + '\n')

shutil.rmtree('extract')
