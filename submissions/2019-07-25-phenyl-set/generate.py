import cmiles
import fragmenter
import openeye
import json
import getpass

oemols = fragmenter.chemi.file_to_oemols('biphenyls.smi')
bond_maps = []
wbos = []
td_json = {}
for mol in oemols:
    # Add atom map
    cmiles.utils.add_atom_map(mol)
    charged_mol = fragmenter.chemi.get_charges(mol, keep_confs=-1)
    # Find rotatable bond
    for bond in charged_mol.GetBonds():
        if bond.IsRotor():
            map_idx = (bond.GetBgn().GetMapIdx(), bond.GetEnd().GetMapIdx())
            bond_maps.append(map_idx)
            wbos.append(bond.GetData('WibergBondOrder'))
    # Find torsion around this bond
    mapped_smiles = openeye.oechem.OEMolToSmiles(mol)
    torsion = fragmenter.torsions.find_torsion_around_bond(charged_mol, map_idx)
    qcarchive_mols = cmiles.utils.mol_to_map_ordered_qcschema(charged_mol, mapped_smiles)
    job_idx = cmiles.utils.to_canonical_label(mapped_smiles, torsion)
    cmiles_identifiers = cmiles.get_molecule_ids(mapped_smiles)
    # Sanity check that the mapped SMILES are the same
    if not mapped_smiles == cmiles_identifiers['canonical_isomeric_explicit_hydrogen_mapped_smiles']:
        print('mapped SMILES do not match. {}, {}'.format(mapped_smiles,
                                                          cmiles_identifiers['canonical_isomeric_explicit_hydrogen_smiles']))
    td_json[job_idx] = {'dihedral': [torsion],
                        'grid': [15],
                        'cmiles_identifiers': cmiles_identifiers,
                        'input_molecules': qcarchive_mols,
                        'provenance': {'fragmenter_version': fragmenter.__version__,
                                       'openeye_version': openeye.__version__,
                                       'username': getpass.getuser(),}}


# Generate figure to visualize input molecules with OE WBO
fragmenter.chemi.to_pdf(oemols, fname='biphenyls.pdf', cols=2, rows=3, bond_map_idx=bond_maps, bo=wbos)

# Save td inputs
with open('biphenyls_set_input.json', 'w') as f:
    json.dump(td_json, f, indent=2, sort_keys=True)