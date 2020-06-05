"""Generates the torsion drive inputs from the molecules selected by DANCE."""
import gzip
import json
import logging
import tempfile
from collections import defaultdict
from typing import Dict, List

import cmiles
from openeye import oechem
from openforcefield.topology import Molecule, Topology
from openforcefield.typing.engines.smirnoff import ForceField

import fragmenter

INPUT_SMILES = "t142_random.smi"

# Configure logger.
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name

# smirnoff99Frosst force field
FORCE_FIELD = ForceField('test_forcefields/smirnoff99Frosst.offxml')


def calculate_mol_params(mol: oechem.OEMol) -> Dict[str, List[List[int]]]:
    """Calculates parameters of the given molecule.

    Returns a dict where the keys are parameter ids and the values are lists
    of indices where the parameter occurs (each entry in the list is itself a
    list because the parameter involves multiple atoms).
    """
    oechem.OEAddExplicitHydrogens(mol)
    off_mol = Molecule.from_openeye(mol, allow_undefined_stereo=True)
    topology = Topology.from_molecules(off_mol)
    molecule_force_list = FORCE_FIELD.label_molecules(topology)

    params = defaultdict(list)
    for _, force_dict in molecule_force_list[0].items():
        for (atom_indices, parameter) in force_dict.items():
            params[parameter.id].append(atom_indices)

    return params


def save_smiles(smiles, filename):
    """Write smiles str to smi file"""
    with open(filename, 'w') as f:
        for smi in smiles:
            f.write(smi + '\n')


def make_json(smiles):
    """
    Takes in a list of smiles strings and expands the tautomeric and isomeric
    state of the molecules and generates a .json file from these molecules. The
    functional also generates .smi files that record processed canonical smiles,
    duplicates, omega failures, cmiles failures, and skipped ions.

    Also takes in torsion indices and writes them in the `atom_indices` field.

    Copied from Jessica Maat's script:
    https://github.com/openforcefield/qca-dataset-submission/blob/master/2020-03-20-OpenFF-Gen-2-Optimization-Set-3-Pfizer-Discrepancy/01_generateOptDS.py#L350

    Input:
        smiles: List of smiles strings
    Return:
        optSmiles: List of smiles that are used as optimization inputs in the
            .json file.
    """
    with tempfile.NamedTemporaryFile('w+', suffix='.smi') as tmp:
        for line in smiles:
            tmp.writelines(line + '\n')
        tmp.seek(0)
        temp_name = tmp.name
        print(tmp.name)
        oemols = fragmenter.chemi.file_to_oemols(temp_name)

        optimization_input = []
        processed_canonical_smiles = []
        skipped = []
        duplicates = []  # duplicate states
        omega_failures = []
        cmiles_failures = []

    # SDF file for writing all conformations.
    ofs = oechem.oemolostream('optimization_inputs.sdf')

    optimization_count = 0
    for mol in oemols:
        # Filter out single atom molecules
        if mol.GetMaxAtomIdx() == 1:
            skipped.append(cmiles.utils.mol_to_smiles(mol, mapped=False))
            continue

        # Expand protonation states and stereoisomers
        states = fragmenter.states.enumerate_states(mol,
                                                    stereoisomers=False,
                                                    tautomers=False)
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
                logging.info('Found duplicate canonical SMILES {}'.format(
                    canonical_smiles))
                duplicates.append(canonical_smiles)
                continue
            else:
                processed_canonical_smiles.append(canonical_smiles)

            # Calculate indices of the parameter. We have to recalculate because
            # indices change when we use different SMILES.
            mol_from_cmiles = oechem.OEMol()
            oechem.OESmilesToMol(
                mol_from_cmiles,
                cmiles_ids['canonical_isomeric_explicit_hydrogen_mapped_smiles']
            )
            torsion_indices = calculate_mol_params(mol_from_cmiles)['t142'][0]

            # Generate molecule using mapped SMILES
            mapped_smiles = cmiles_ids[
                'canonical_isomeric_explicit_hydrogen_mapped_smiles']
            m = cmiles.utils.load_molecule(s)

            try:
                # Omega fails for some molecules.
                conformers = fragmenter.chemi.generate_conformers(m)
            except RuntimeError:
                logging.info(
                    'Omega failed to generate conformers for {}'.format(
                        cmiles_ids['canonical_isomeric_smiles']))
                # Omega failed
                omega_failures.append(cmiles_ids['canonical_isomeric_smiles'])
                continue

            qcschema_molecules = [
                cmiles.utils.mol_to_map_ordered_qcschema(conf, mapped_smiles)
                for conf in conformers.GetConfs()
            ]
            optimization_input.append({
                'cmiles_identifiers': cmiles_ids,
                'atom_indices': [torsion_indices],
                'initial_molecules': qcschema_molecules,
            })
            optimization_count += len(qcschema_molecules)
            # Write to SDF
            oechem.OEWriteMolecule(ofs, conformers)

    with gzip.open('optimization_inputs.json.gz', 'w') as f:
        f.write(json.dumps(optimization_input, indent=2).encode('utf-8'))

    ofs.close()

    save_smiles(processed_canonical_smiles, 'optimization_inputs.smi')
    save_smiles(duplicates, 'duplicates.smi')
    save_smiles(omega_failures, 'omega_failures.smi')
    save_smiles(cmiles_failures, 'cmiles_failures.smi')
    save_smiles(skipped, 'skipped_ions.smi')
    print("Number of unique molecules optimized:" + str(len(oemols)))
    print("Final optimization count is:" + str(optimization_count))

    file1 = open("finalCounts.txt", "w")  #write mode
    file1.write("Number of molecules optimized:" + str(len(oemols)) + '\n')
    file1.write("Final optimization count with expanded states is:" +
                str(optimization_count) + '\n')
    file1.close()

    opt_smiles = []
    for mol in oemols:
        opt_smiles.append(oechem.OEMolToSmiles(mol))

    return opt_smiles


def main():
    """Reads in the molecules, saves them to the JSON file."""

    # The make_json function takes in a list of SMILES strings, so we read them
    # in from the file here.
    smiles = []
    smiles_file = open(INPUT_SMILES, 'r')
    for line in smiles_file:
        smiles.append(line.split()[0])
    smiles_file.close()

    make_json(smiles)


if __name__ == "__main__":
    main()
