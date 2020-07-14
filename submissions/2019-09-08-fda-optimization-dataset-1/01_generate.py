import fragmenter
import cmiles
import json
import warnings
import logging
from openeye import oechem
from openeye import oemolprop

def save_smiles(smiles, filename):
    """Write smiles str to smi file"""
    if filename.endswith('.gz'):
        import gzip
        with gzip.open(filename, 'w') as f:
            for smi in smiles:
                f.write( (smi + '\n').encode('utf-8') )
    else:
        with open(filename, 'w') as f:
            for smi in smiles:
                f.write(smi + '\n')

def read_smiles(filename):
    with open(filename, 'r') as f:
        smiles = f.read().split('\n')
    print(smiles)
    return smiles

logging.getLogger().setLevel(logging.INFO)

filterfile = oechem.oeifstream('oechem-filterfile')
filter = oemolprop.OEFilter(filterfile)

MAX_CONFS = 20

# Read SMILES
oemols = fragmenter.chemi.file_to_oemols('fda.mol2.gz')

optimization_input = []
processed_smiles = []
skipped = []
duplicates = [] # duplicate protonation/tautomeric states
omega_failures = []
cmiles_failures = []

# Write out SDF file of all conformations
ofs = oechem.oemolostream('optimization_inputs.sdf.gz')

# Drop duplicate input molecules
print('De-duplicating input molecules')
smiles_set = set()
new_oemols = list()
for oemol in oemols:
    smiles = oechem.OEMolToSmiles(oemol)
    if oemol not in smiles_set:
        new_oemols.append(oemol)
        smiles_set.add(smiles)
    else:
        duplicates.append(smiles)
print(f'Retained {len(new_oemols)} unique molecules out of {len(oemols)}')
oemols = new_oemols

for index, input_oemol in enumerate(oemols):
    # Apply filter criteria
    if not filter(input_oemol):
        skipped.append(smiles)
        continue

    smiles = oechem.OEMolToSmiles(input_oemol)
    print(f'input molecule {index} / {len(oemols)} : {smiles}')

    # Generate SMILES to use for CMILES
    try:
        smiles = cmiles.utils.mol_to_smiles(input_oemol, isomeric=True, mapped=False, explicit_hydrogen=True)
    except Exception as e:
        cmiles_failures.append(smiles)
        print(e)
        continue

    # Generate mapped CMILES for molecules with all explicit hydrogens
    try:
        cmiles_ids = cmiles.get_molecule_ids(smiles)
    except Exception as e:
        cmiles_failures.append(smiles)
        continue

    # Generate molecule using mapped SMILES
    mapped_smiles = cmiles_ids['canonical_isomeric_explicit_hydrogen_mapped_smiles']
    oemol = cmiles.utils.load_molecule(mapped_smiles)

    # Generate conformers
    try:
        # Omega fails for some molecules.
        conformers = fragmenter.chemi.generate_conformers(oemol, max_confs=MAX_CONFS)
    except RuntimeError:
        logging.info('Omega failed to generate conformers for {}'.format(cmiles_ids['canonical_isomeric_smiles']))
        # Omega failed
        omega_failures.append(cmiles_ids['canonical_isomeric_smiles'])
        continue
    print(f'    {conformers.NumConfs()} confomers generated')

    # Convert to QCSchema
    qcschema_molecules = [cmiles.utils.mol_to_map_ordered_qcschema(conf, mapped_smiles) for conf in conformers.GetConfs()]

    # Append to QCSchema-based optimization input
    optimization_input.append({'initial_molecules': qcschema_molecules,
                               'cmiles_identifiers': cmiles_ids})

    # Write to SDF
    oechem.OEWriteMolecule(ofs, conformers)
    processed_smiles.append(mapped_smiles)

import gzip
with gzip.open('optimization_inputs.json.gz', 'w') as f:
    f.write(json.dumps(optimization_input, indent=2, sort_keys=True).encode('utf-8'))

ofs.close()

save_smiles(processed_smiles, 'optimization_inputs.smi.gz')
save_smiles(duplicates, 'duplicates.smi.gz')
save_smiles(omega_failures, 'omega_failures.smi.gz')
save_smiles(cmiles_failures, 'cmiles_failures.smi.gz')
save_smiles(skipped, 'skipped.smi.gz')
