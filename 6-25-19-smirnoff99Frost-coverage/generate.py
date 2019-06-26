from fragmenter import fragment, chemi
import cmiles
import json
import warnings

# Read SMILES
oemols = chemi.file_to_oemols('chosen.smi')

optimization_input = []
skipped = []
omega_failures = []
cmiles_failures = []
for mol in oemols:
    # Filter out single atom molecules
    if mol.GetMaxAtomIdx() == 1:
        skipped.append(cmiles.utils.mol_to_smiles(mol))
        continue

    # Expand protonation states and stereoisomers
    states = fragment.expand_states(mol, stereoisomers=True, protonation=True, tautomers=False)
    for s in states:
        try:
            cmiles_ids = cmiles.get_molecule_ids(s)
        except:
            cmiles_failures.append(s)
            continue
        mapped_smiles = cmiles_ids['canonical_isomeric_explicit_hydrogen_mapped_smiles']
        m = cmiles.utils.load_molecule(s)
        try:
            conformers = chemi.generate_conformers(m)
        except RuntimeError:
            warnings.warn('Omega failed to generate conformers for {}'.format(cmiles_ids['canonical_isomeric_smiles']))
            # Omega failed
            omega_failures.append(cmiles_ids['canonical_isomeric_smiles'])
            continue
        qcschema_molecules = [cmiles.utils.mol_to_map_ordered_qcschema(conf, mapped_smiles) for conf in conformers.GetConfs()]
        optimization_input.append({'initial_molecules': qcschema_molecules,
                                   'cmiles_identifiers': cmiles_ids})


with open('optimization_inputs.json', 'w') as f:
    json.dump(optimization_input, f, indent=2, sort_keys=True)

with open('omega_failures.json', 'w') as f:
    json.dump(omega_failures, f, indent=2, sort_keys=True)

with open('cmiles_failures.json', 'w') as f:
    json.dump(cmiles_failures, f, indent=2, sort_keys=True)

with open('skipped_ions.json', 'w') as f:
    json.dump(skipped, f, indent=2, sort_keys=True)