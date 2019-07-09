#!/usr/bin/env python

import json
from collections import Counter, defaultdict
import tarfile

import cmiles
from openforcefield.typing.engines.smirnoff import ForceField
from openforcefield.topology import Molecule as Off_Molecule
from openforcefield.topology import Topology as Off_Topology

from qcelemental.models import Molecule



def read_aggregate_molecules(input_json):
    """ Extract the molecules and the index of them from the input json file
    aggregate molecules with the same index into a list

    Parameters
    ----------
    input_json: str,
        JSON file name to the output json of generate.py
        The data in the json file should be a list of {'initial_molecules': [..], 'cmiles_identifiers':{}}.

    Returns
    -------
    molecules_list_dict: dict
        The dictionary maps the index of a molecule to a Molecule object. e.g.
        {
            index1: [Molecule_json1a, Molecule_json1b, ..],
            index2: [Molecule_json2a, Molecule_json2b, ..],
        }

    molecule_attributes: dict
        The dicitonary maps the index of a molecule to the attributes of the molecule, e.g.
        {
            index1: {'canonical_explicit_hydrogen_smiles': .., 'canonical_isomeric_smiles': .., ..}
        }

    Note
    ----
    1. The mdata['cmiles_identifiers']['canonical_isomeric_smiles'] is selected as the index.
    2. For molecules have the same "canonical_isomeric_smiles", we use index-1, index-2 to distinguish them.
    """
    molecules_list_dict = defaultdict(list)
    molecule_attributes = {}
    # open json file
    if input_json.endswith(".tar") or input_json.endswith(".tar.gz"):
        extract_file = input_json.replace(".gz", "").replace(".tar", ".json")
        with tarfile.open(input_json, 'r') as infile:
            molecule_data_list = json.load(infile.extractfile(extract_file))
    else:
        with open(input_json) as infile:
            molecule_data_list = json.load(infile)
    # put molecules and attributes into molecules_list_dict
    molecule_hash = defaultdict(set) # use a dictionary to remove duplicates
    for mdata in molecule_data_list:
        initial_molecules = mdata['initial_molecules']
        cmiles_ids = mdata['cmiles_identifiers']
        index = cmiles_ids['canonical_isomeric_smiles']
        molecule_attributes[index] = cmiles_ids
        for m_json in initial_molecules:
            m_hash = Molecule.from_data(m_json).get_hash()
            # find duplicated molecules using their hash and skip them
            if m_hash not in molecule_hash[index]:
                molecule_hash[index].add(m_hash)
                molecules_list_dict[index].append(m_json)
    return molecules_list_dict, molecule_attributes


def smirnoff_analyze_torsions(forcefield, off_mol):
    """
    Compute the coverage of all torsions in this molecule

    Parameters
    ----------
    forcefield: openforcefield.typing.engines.smirnoff.ForceField
        The forcefield object for computing coverage
    off_mol: openforcefield.topology.Molecule
        The molecule object for computing torsions coverage

    Returns
    -------
    torsions_coverage: dict
        Key is smirks for the torsion, value is a list of torsion indices
        {SMIRKs: [(0,1,2,3), (2,4,6,7), ..] }
    """
    torsions_coverage = defaultdict(list)
    off_top = Off_Topology.from_molecules(off_mol)
    for torsion_indices, torsion_param in forcefield.label_molecules(off_top)[0]['ProperTorsions'].items():
        torsions_coverage[torsion_param.smirks].append(torsion_indices)
    return torsions_coverage


def select_torsions(molecules_list_dict, molecule_attributes, forcefield, target_coverage=3):
    """
    select index from molecules_list_dict that covers different SMIRKs

    Parameters
    ----------
    molecules_list_dict: dict
        result of read_aggregate_molecules()
    molecule_attributes: dict
        result of read_aggregate_molecules()
    forcefield: openforcefield.typing.engines.smirnoff.ForceField
        The forcefield object for computing coverage
    target_coverage: int
        target number of coverage of each torsion. After reaching this number, later torsions for the same SMIRKs will be ignored

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
    torsions_dict = {}
    smirks_torsions_counter = Counter()
    i_mol = 0
    for mol_index, mol_attr in molecule_attributes.items():
        print(f'{i_mol:<7d}: {mol_index}')
        i_mol += 1
        mapped_smiles = mol_attr['canonical_isomeric_explicit_hydrogen_mapped_smiles']
        # round trip from QCFractal molecule to OpenEye molecule then to Off Molecule
        # this is needed for now to ensure atom indices are consistent
        qcjson_mol = molecules_list_dict[mol_index][0]
        oemol = cmiles.utils.load_molecule(qcjson_mol)
        off_mol = Off_Molecule.from_openeye(oemol)
        torsions_coverage = smirnoff_analyze_torsions(forcefield, off_mol)
        for smirks, torsion_idx_list in torsions_coverage.items():
            for atom_indices in torsion_idx_list:
                if smirks_torsions_counter[smirks] < target_coverage:
                    smirks_torsions_counter[smirks] += 1
                    canonical_torsion_index = cmiles.utils.to_canonical_label(mapped_smiles, atom_indices)
                    torsions_dict[canonical_torsion_index] = {
                        'initial_molecules': molecules_list_dict[mol_index],
                        'atom_indices': [ atom_indices ],
                        'attributes': mol_attr,
                    }
                    print(f"  - torsion {atom_indices} added for smirks {smirks}")
                else:
                    print(f"  - torsion {atom_indices} skipped because {smirks} have {smirks_torsions_counter[smirks]} already")
    print("\n## Selected Torsion Coverage ##\n" + '-'*90)
    ff_torsion_param_list = forcefield.get_parameter_handler('ProperTorsions').parameters
    n_covered = 0
    for param in ff_torsion_param_list:
        count = smirks_torsions_counter[param.smirks]
        print(f"{param.smirks:80s} : {count:7d}")
        if count > 0:
            n_covered += 1
    print('-'*90)
    print(f'{n_covered} / {len(ff_torsion_param_list)} torsion SMIRKs covered')
    return torsions_dict



print("## Extracting molecules ##")
molecules_list_dict, molecule_attributes = read_aggregate_molecules("optimization_inputs.tar.gz")
print(f'{len(molecules_list_dict)} unique molecules found')

print("## Selecting torsions ##")

forcefield = ForceField('smirnoff99Frosst_experimental.offxml')

torsions_dict = select_torsions(molecules_list_dict, molecule_attributes, forcefield, target_coverage=5)

print("## Writing selected_torsions.json ##")
with open('selected_torsions.json', 'w') as jsonfile:
    json.dump(torsions_dict, jsonfile, indent=2)


print(f'{len(torsions_dict)} torsions selected')