import cmiles
import fragmenter
import json

mol_id = cmiles.get_molecule_ids('OCCO', strict=False)
mapped_smiles = (mol_id['canonical_isomeric_explicit_hydrogen_mapped_smiles'])
mol = cmiles.utils.load_molecule(mapped_smiles)

torsions = fragmenter.torsions.find_torsions(mol)

dihedrals_list = [torsions['internal']['torsion_0'], torsions['terminal']['torsion_0']]
single_conformer = fragmenter.chemi.generate_conformers(mol, max_confs=1)
mult_conformers_grid = fragmenter.chemi.generate_grid_conformers(mol, dihedrals=dihedrals_list, intervals=[90, 120])

qm_mol_single_conf = cmiles.utils.mol_to_map_ordered_qcschema(single_conformer, mapped_smiles)
qm_mol_mult_conf = [cmiles.utils.mol_to_map_ordered_qcschema(conf, mapped_smiles) for conf in mult_conformers_grid.GetConfs()]

job_index_1d = cmiles.utils.to_canonical_label(mapped_smiles, dihedrals_list[0])
job_index_2d = job_index_1d + ',' + cmiles.utils.to_canonical_label(mapped_smiles , dihedrals_list[1])
job_index_1d_mult = job_index_1d + '_' + str(len(qm_mol_mult_conf))
job_index_2d_mult = job_index_2d + '_' + str(len(qm_mol_mult_conf))


torsion_drive_inputs = {
    job_index_1d: {
        'dihedral': [dihedrals_list[0]],
        'grid': [15],
        'input_molecules': qm_mol_single_conf,
        'cmiles_identifiers': mol_id
    },
    job_index_2d: {
        'dihedral': dihedrals_list,
        'grid': [15, 15],
        'input_molecules': qm_mol_single_conf,
        'smiles_identifiers': mol_id

    },
    job_index_1d_mult: {
        'dihedral': [dihedrals_list[0]],
        'grid': [15],
        'input_molecules': qm_mol_mult_conf,
        'cmiles_identifiers': mol_id
    },
    job_index_2d_mult: {
        'dihedral': dihedrals_list,
        'grid': [15, 15],
        'input_molecules': qm_mol_mult_conf,
        'cmiles_identifiers': mol_id
    }
}

with open('torsiondrive_inputs.json', 'w') as f:
    json.dump(torsion_drive_inputs, f, indent=2, sort_keys=True)