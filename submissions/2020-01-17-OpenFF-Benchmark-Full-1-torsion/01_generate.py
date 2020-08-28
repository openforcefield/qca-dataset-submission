import fragmenter
import cmiles
import json
import warnings
import logging
from openeye import oechem

logging.getLogger().setLevel(logging.INFO)

def save_smiles(smiles, filename):
    """Write smiles str to smi file"""
    with open(filename, 'w') as f:
        for smi in smiles:
            f.write(smi + '\n')

# Read smi file containing SMILES of benchmark primary set
oemols =  fragmenter.chemi.file_to_oemols('full_set_filtered.smi')

optimization_input = []
processed_canonical_smiles = []
skipped = []
duplicates = [] # duplicate states
omega_failures = []
cmiles_failures = []

for mol in oemols:
    # Filter out single atom molecules
    if mol.GetMaxAtomIdx() == 1: 
        skipped.append(cmiles.utils.mol_to_smiles(mol, mapped=False))
        continue
    
    # Expand protonation states and stereoisomers (fragment.expand_states -> states.enumerate_states, 'protonation=True' removed)
    states = fragmenter.states.enumerate_states(mol, stereoisomers=True, tautomers=False)
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
            # Omega failes for some molecules.
            conformers = fragmenter.chemi.generate_conformers(m)
        except RuntimeError:
            logging.info('Omega failed to generate conformers for {}'.format(cmiles_ids['canonical_isomeric_smiles']))
            # Omega failed
            omega_failures.append(cmiles_ids['canonical_isomeric_smiles'])
            continue
        qcschema_molecules = [cmiles.utils.mol_to_map_ordered_qcschema(conf, mapped_smiles) for conf in conformers.GetConfs()]
        optimization_input.append({'initial_molecules': qcschema_molecules,
                                   'cmiles_identifiers': cmiles_ids})

with open('optimization_inputs.json', 'w') as f:
    json.dump(optimization_input, f, indent=2, sort_keys=True)

save_smiles(processed_canonical_smiles, 'optimization_inputs.smi')
save_smiles(duplicates, 'duplicates.smi')
save_smiles(omega_failures, 'omega_failures.smi')
save_smiles(cmiles_failures, 'cmiles_failures.smi')
save_smiles(skipped, 'skipped_ions.smi')
