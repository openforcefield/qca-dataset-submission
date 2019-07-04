from fragmenter import fragment, chemi
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
oemols = chemi.file_to_oemols('chosen_supplemented.smi')

optimization_input = []
skipped = []
omega_failures = []
cmiles_failures = []

for mol in oemols:
    # Filter out single atom molecules
    if mol.GetMaxAtomIdx() == 1:
        skipped.append(cmiles.utils.mol_to_smiles(mol, mapped=False))
        continue

    # Expand protonation states and stereoisomers
    states = fragment.expand_states(mol, stereoisomers=True, protonation=True, tautomers=False)
    for s in states:
        # Some states have valences that rdkit does not accept.
        try:
            cmiles_ids = cmiles.get_molecule_ids(s)
        except:
            cmiles_failures.append(s)
            continue
        mapped_smiles = cmiles_ids['canonical_isomeric_explicit_hydrogen_mapped_smiles']
        m = cmiles.utils.load_molecule(s)
        try:
            # Omega fails for some molecules.
            conformers = chemi.generate_conformers(m)
        except RuntimeError:
            logging.info('Omega failed to generate conformers for {}'.format(cmiles_ids['canonical_isomeric_smiles']))
            # Omega failed
            omega_failures.append(cmiles_ids['canonical_isomeric_smiles'])
            continue
        qcschema_molecules = [cmiles.utils.mol_to_map_ordered_qcschema(conf, mapped_smiles) for conf in conformers.GetConfs()]
        optimization_input.append({'initial_molecules': qcschema_molecules,
                                   'cmiles_identifiers': cmiles_ids})
# Add quacpac failure
s = read_smiles('quacpac_failures.smi')[0]
cmiles_ids = cmiles.get_molecule_ids(s, strict=False)
mapped_smiles = cmiles_ids['canonical_isomeric_explicit_hydrogen_mapped_smiles']
m = cmiles.utils.load_molecule(s)
conformers = chemi.generate_conformers(m)
qcschema_molecules = [cmiles.utils.mol_to_map_ordered_qcschema(conf, mapped_smiles) for conf in conformers.GetConfs()]
optimization_input.append({'initial_molecules': qcschema_molecules,
                                   'cmiles_identifiers': cmiles_ids})

with open('optimization_inputs.json', 'w') as f:
    json.dump(optimization_input, f, indent=2, sort_keys=True)

save_smiles(omega_failures, 'omega_failures.smi')
save_smiles(cmiles_failures, 'cmiles_failures.smi')
save_smiles(skipped, 'skipped_ions.smi')
