#!/usr/bin/env python

"""
Generate a PDF of all small molecules in the JSON dataset.
"""

from fragmenter import fragment, chemi
import json
import gzip

# Read the compressed dataset
with gzip.open('optimization_inputs.json.gz', 'r') as f:
    data = json.loads(f.read().decode('utf-8'))

# Extract SMILES
smiles_list = [ item['cmiles_identifiers']['canonical_isomeric_smiles'] for item in data ]

# Build OEMols
from openeye import oechem
oemols = list()
for smiles in smiles_list:
    oemol = oechem.OEMol()
    oechem.OESmilesToMol(oemol, smiles)
    oemols.append(oemol)

# Generate a PDF of all molecules in the set
pdf_filename = 'VEHICLe.pdf'
chemi.to_pdf(oemols, pdf_filename)
