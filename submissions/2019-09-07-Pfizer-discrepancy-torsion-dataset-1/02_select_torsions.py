#!/usr/bin/env python

import json
from collections import Counter, defaultdict
import tarfile

import cmiles
from openforcefield.topology import Molecule

def enumerate_torsions(oemol):
    """Enumerate torsions to drive using the OpenEye toolkit.

    Enumerates torsions that are
    * Marked as rotatable
    * Not in rings
    * Prioritizes heaviest atoms involved in the torsion

    Returns
    -------
    torsion_idx_list : list of (int,int,int,int)
        List of rotatable torsions in molecule, if any
    """
    torsion_idx_list = list()

    for oebond in oemol.GetBonds():
        if not oebond.IsRotor(): continue
        if oebond.IsInRing(): continue
        # Find heaviest (i,j,k,l) tuple
        max_mass = -1
        optimal_torsion_idx = None
        jatom = oebond.GetBgn()
        katom = oebond.GetEnd()
        for iatom in jatom.GetAtoms():
            for latom in katom.GetAtoms():
                torsion_idx = (iatom.GetIdx(), jatom.GetIdx(), katom.GetIdx(), latom.GetIdx())
                if len(set(torsion_idx)) != 4: continue
                mass = iatom.GetAtomicNum() + latom.GetAtomicNum()
                if mass > max_mass:
                    max_mass = mass
                    optimal_torsion_idx = torsion_idx
        torsion_idx_list.append(optimal_torsion_idx)

    return torsion_idx_list

def read_molecules(input_json):
    """
    Parameters
    ----------
    input_json: str,
        JSON file name to the output json of generate.py (prepared as if for an OptimizationDataset)
        The data in the json file should be a list of {'initial_molecules': [..], 'cmiles_identifiers':{}}.

    Returns
    -------
    molecules_list_dict: dict
        The dictionary maps the index of a molecule to a Molecule object. e.g.
        {
            index1: [Molecule_json1a, Molecule_json1b, ..],
            index2: [Molecule_json2a, Molecule_json2b, ..],
        }

    """
    if input_json.endswith(".tar") or input_json.endswith(".tar.gz"):

        extract_file = input_json.replace(".gz", "").replace(".tar", ".json")
        with tarfile.open(input_json, 'r') as infile:

            molecule_data_list = json.load(infile.extractfile(extract_file))

    if input_json.endswith(".gz"):

        import gzip
        with gzip.open(input_json, 'r') as infile:

            molecule_data_list = json.loads(infile.read().decode('utf-8'))

    else:
        with open(input_json) as infile:
            molecule_data_list = json.load(infile)

    print(f'{len(molecule_data_list)} molecules read')
    return molecule_data_list

def generate_selected_torsions(input_json):
    """Identify torsions that can be driven.

    Parameters
    ----------
    input_json: str,
        JSON file name to the output json of generate.py (prepared as if for an OptimizationDataset)
        The data in the json file should be a list of {'initial_molecules': [..], 'cmiles_identifiers':{}}.

    Returns
    -------
    torsions_dict: dict
        Dictionary for selected torsions, has this structure:
        {
            canonical_torsion_index1: {
                'initial_molecules': [ Molecule1a, Molecule1b, .. ],
                'atom_indices': [ (0,1,2,3) ],
                'attributes': {'canonical_explicit_hydrogen_smiles': .., 'canonical_isomeric_smiles': .., ..}
            },
            ..
        }

    Note
    ----
    The 'atom_indices' in return dict value is a list with only one item, because we select only 1-D torsion for now.

    """
    molecule_data_list = read_molecules(input_json)

    # generate torsion_dict
    torsions_dict = {}
    ntorsions = 0
    for mol_index, json_mol in enumerate(molecule_data_list):
        mapped_smiles = json_mol['cmiles_identifiers']['canonical_isomeric_explicit_hydrogen_mapped_smiles']
        print(f'{mol_index} : {mapped_smiles}')
        # round trip from QCFractal molecule to OpenEye molecule then to Off Molecule
        # this is needed for now to ensure atom indices are consistent
        qcjson_mol = json_mol['initial_molecules'][0]
        oemol = cmiles.utils.load_molecule(qcjson_mol)
        off_mol = Molecule.from_openeye(oemol, allow_undefined_stereo=True)
        torsion_idx_list = enumerate_torsions(oemol)
        for atom_indices in torsion_idx_list:
            torsions_dict[ntorsions] = {
                'initial_molecules': [ qcjson_mol ],
                'atom_indices': [ atom_indices ],
                'attributes': json_mol['cmiles_identifiers'],
            }
            print(f"  - torsion {atom_indices} added")
            ntorsions += 1

    print(f'{ntorsions} torsions added')
    return torsions_dict


print("## Extracting molecules ##")
torsions_dict = generate_selected_torsions("optimization_inputs.json.gz")
print("## Writing selected_torsions.json ##")
import gzip
with gzip.open('selected_torsions.json.gz', 'w') as jsonfile:
    jsonfile.write(json.dumps(torsions_dict, indent=2, sort_keys=True).encode('utf-8'))

print(f'{len(torsions_dict)} torsions selected')
