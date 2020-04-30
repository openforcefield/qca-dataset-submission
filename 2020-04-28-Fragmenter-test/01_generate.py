import json
import cmiles
import fragmenter

with open('example-molecules.json') as f:
    frags_mols = json.load(f)

torsion_dict = {}
keep_track = {}
for mol in frags_mols:
    print(mol)
    keep_track[mol] = {}
    with open('../benchmark_fragmentation_schemes/{}/{}_wbo_dists.json'.format(mol, mol)) as f:
        data = json.load(f)
    with open('../benchmark_fragmentation_schemes/{}/{}_pfizer_wbo_dists.json'.format(mol, mol)) as f:
        data_pfizer = json.load(f)

    for bond in frags_mols[mol]:
        bond_idx = fragmenter.utils.deserialize_bond(bond)
        for frag in ('parent', '0.03'):
            print(frag)
            smiles = data[bond][frag]['frag']
            oemol = cmiles.utils.load_molecule(smiles)
            torsion = fragmenter.torsions.find_torsion_around_bond(oemol, bond=bond_idx)
            conformers = fragmenter.chemi.generate_conformers(oemol, max_confs=1, strict_stereo=False, strict_types=False)
            cmiles_ids = cmiles.get_molecule_ids(smiles, strict=False)
            mapped_smiles = cmiles_ids['canonical_isomeric_explicit_hydrogen_mapped_smiles']
            qcarchive_mols = [cmiles.utils.mol_to_map_ordered_qcschema(conf, mapped_smiles) for conf in conformers.GetConfs()]
            job_idx = cmiles.utils.to_canonical_label(smiles, torsion)
            torsion_dict[job_idx] = {
                'initial_molecules': qcarchive_mols,
                'dihedral': [torsion],
                'grid': [15],
                'cmiles_identifiers': cmiles_ids
            }

            keep_track[mol][frag] = job_idx
        frag = 'pfizer'
        print(frag)
        smiles = data_pfizer[bond]['frag']
        oemol = cmiles.utils.load_molecule(smiles)
        torsion = fragmenter.torsions.find_torsion_around_bond(oemol, bond=bond_idx)
        conformers = fragmenter.chemi.generate_conformers(oemol, max_confs=1, strict_stereo=False, strict_types=False)
        cmiles_ids = cmiles.get_molecule_ids(smiles, strict=False)
        mapped_smiles = cmiles_ids['canonical_isomeric_explicit_hydrogen_mapped_smiles']
        qcarchive_mols = [cmiles.utils.mol_to_map_ordered_qcschema(conf, mapped_smiles) for conf in
                          conformers.GetConfs()]
        job_idx = cmiles.utils.to_canonical_label(smiles, torsion)
        torsion_dict[job_idx] = {
            'initial_molecules': qcarchive_mols,
            'dihedral': [torsion],
            'grid': [15],
            'cmiles_identifiers': cmiles_ids
        }

    keep_track[mol][frag] = job_idx

with open('selected-torsions.json', 'w') as jsonfile:
    json.dump(torsion_dict, jsonfile, indent=2)

with open('torsiondrive-job-idx.json', 'w') as jsonfile:
    json.dump(keep_track, jsonfile, indent=2)


