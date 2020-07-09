import fragmenter
import cmiles
import json
import warnings
import logging
from openeye import oechem

def save_smiles(smiles, filename):
    """Write smiles str to smi file"""
    with open(filename, 'w') as f:
        for smi in smiles:
            f.write(smi + '\n')

def read_smiles(filename):
    with open(filename, 'r') as f:
        smiles = f.read().split('\n')
    print(smiles)
    return smiles

logging.getLogger().setLevel(logging.INFO)

# Read SMILES
oemols = fragmenter.chemi.file_to_oemols('PFE-OFF-100Frags.smi.gz')

optimization_input = []
processed_canonical_smiles = []
skipped = []
duplicates = [] # duplicate states
omega_failures = []
cmiles_failures = []

# Write out SDF file of all conformations
ofs = oechem.oemolostream('optimization_inputs.sdf')

for mol in oemols:
    # Filter out single atom molecules
    if mol.GetMaxAtomIdx() == 1:
        skipped.append(cmiles.utils.mol_to_smiles(mol, mapped=False))
        continue

    # Expand protonation states and stereoisomers
    states = fragmenter.states.enumerate_states(mol, stereoisomers=False, tautomers=False)
    for s in states:
        # Some states have valences that rdkit does not accept.
        try:
            cmiles_ids = cmiles.get_molecule_ids(s)
        except:
            cmiles_failures.append(s)
            continue

        # Drop duplicates
        canonical_smiles = cmiles_ids['canonical_smiles']
        if canonical_smiles in processed_canonical_smiles:
            logging.info('Found duplicate canonical SMILES {}'.format(canonical_smiles))
            duplicates.append(canonical_smiles)
            continue
        else:
            processed_canonical_smiles.append(canonical_smiles)

        # Generate molecule using mapped SMILES
        mapped_smiles = cmiles_ids['canonical_isomeric_explicit_hydrogen_mapped_smiles']
        m = cmiles.utils.load_molecule(s)
        try:
            # Omega fails for some molecules.
            conformers = fragmenter.chemi.generate_conformers(m)
        except RuntimeError:
            logging.info('Omega failed to generate conformers for {}'.format(cmiles_ids['canonical_isomeric_smiles']))
            # Omega failed
            omega_failures.append(cmiles_ids['canonical_isomeric_smiles'])
            continue
        qcschema_molecules = [cmiles.utils.mol_to_map_ordered_qcschema(conf, mapped_smiles) for conf in conformers.GetConfs()]
        optimization_input.append({'initial_molecules': qcschema_molecules,
                                   'cmiles_identifiers': cmiles_ids})

        # Write to SDF
        oechem.OEWriteMolecule(ofs, conformers)

import gzip
with gzip.open('optimization_inputs.json.gz', 'w') as f:
    f.write(json.dumps(optimization_input, indent=2, sort_keys=True).encode('utf-8'))

ofs.close()

save_smiles(processed_canonical_smiles, 'optimization_inputs.smi')
save_smiles(duplicates, 'duplicates.smi')
save_smiles(omega_failures, 'omega_failures.smi')
save_smiles(cmiles_failures, 'cmiles_failures.smi')
save_smiles(skipped, 'skipped_ions.smi')
