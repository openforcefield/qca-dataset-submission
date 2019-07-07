import json
from openeye import oechem
import cmiles
from fragmenter import torsions, chemi


def get_torsion(mol, rot_bond):
    """
    Find torsion in fragment that corresponds to same rotatable bond in parent molecule
    parameters:
    ----------
    mol: OEMol with map indices to parent molecule
    rot_bond: tuple
        (mapidx_1, mapidx_2)
    return:
    -------
    dihdral: list of dihedral atom indices
    """
    tors = [[tor.a, tor.b, tor.c, tor.d] for tor in oechem.OEGetTorsions(mol)]
    filtered_torsions = torsions.one_torsion_per_rotatable_bond(tors)
    mapped_tors = [[i.GetMapIdx() for i in t] for t in filtered_torsions]
    cbs = [(t[1], t[2]) for t in mapped_tors]
    try:
        dihedral = mapped_tors[cbs.index(rot_bond)]
    except ValueError:
        dihedral = mapped_tors[cbs.index(tuple(reversed(rot_bond)))]
    dihedral = [d-1 for d in dihedral]
    return dihedral


with open('Idelalisib_frags.json', 'r') as f:
    filtered_frags = json.load(f)


# Prepare torsiondrive input
torsiondrive_jobs = {}
keep_track = {}
for frag in filtered_frags:
    can_smiles = filtered_frags[frag]['identifiers']['canonical_isomeric_explicit_hydrogen_mapped_smiles']

    # Find the bond in fragment that corresponds to bond in parent molecule to run torison scan
    for i, bond in enumerate(filtered_frags[frag]['provenance']['central_rot_bond']):
        map_to_parent = filtered_frags[frag]['provenance']['map_to_parent'][i]
        map_to_parent_mol = oechem.OEMol()
        oechem.OESmilesToMol(map_to_parent_mol, map_to_parent)
        des_bond = tuple(int(i) for i in bond.split('[')[-1].split(']')[0].split(','))
        dihedral_in_parent = get_torsion(map_to_parent_mol, des_bond)
        job_label = cmiles.utils.to_canonical_label(mapped_smiles=map_to_parent, labeled_atoms=dihedral_in_parent)

        # Keep track of map idx in parent of labeled fragment
        if job_label not in keep_track:
            keep_track[job_label] = [{'parent_molecule': filtered_frags[frag]['provenance']['parent_molecule'],
                                'central_bond': filtered_frags[frag]['provenance']['central_rot_bond'][i],
                                'map_to_parent': filtered_frags[frag]['provenance']['map_to_parent'][i]}]
        else:
            keep_track[job_label].append({'parent_molecule': filtered_frags[frag]['provenance']['parent_molecule'],
                                'central_bond': filtered_frags[frag]['provenance']['central_rot_bond'][i],
                                'map_to_parent': filtered_frags[frag]['provenance']['map_to_parent'][i]})

        if job_label in torsiondrive_jobs:
            continue
        # Map dihedral onto new atom map
        # Add explicit hydrogen
        oechem.OEAddExplicitHydrogens(map_to_parent_mol)
        cmiles._cmiles_oe.canonical_order_atoms(map_to_parent_mol)
        map_to_parent = cmiles.utils.get_atom_map(map_to_parent_mol, map_to_parent)
        dih = [map_to_parent[d+1] for d in dihedral_in_parent]

        # Generate starting conformations. First add hydrogens
        oechem.OEAddExplicitHydrogens(map_to_parent_mol)
        conformers = chemi.generate_conformers(map_to_parent_mol, strict_types=False, max_confs=10)
        # Genereate QCSchema molecules
        qcschema_mols = [cmiles.utils.mol_to_map_ordered_qcschema(conf, can_smiles) for conf in conformers.GetConfs()]

        torsiondrive_jobs[job_label] = {'initial_molecule': qcschema_mols, 'dihedral': [dih], 'grid': [15],
              'identifiers': filtered_frags[frag]['identifiers'],  'provenance':
                                        {'canonicalization': filtered_frags[frag]['provenance']['canonicalization']}}

with open('stability_benchmark_inputs.json', 'w') as f:
    json.dump(torsiondrive_jobs, f, indent=2, sort_keys=True)

with open('job_label_dict.json', 'w') as f:
    json.dump(keep_track, f, indent=2, sort_keys=True)
