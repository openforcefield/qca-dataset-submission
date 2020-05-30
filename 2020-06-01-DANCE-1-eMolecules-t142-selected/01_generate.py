"""Generates the torsion drive inputs from the molecules selected by DANCE."""
import json

import cmiles
from openeye import oechem

import fragmenter

INPUT_OEB = "t142_selected.oeb"
INPUT_SMILES = "t142_selected.smi"
OUTPUT_FILENAME = "t142_selected_torsions.json"


def main():
    """Reads in the molecules, saves them to the JSON file."""
    input_stream = oechem.oemolistream(INPUT_OEB)
    smiles_file = open(INPUT_SMILES, 'r')
    torsion_dict = {}

    for mol, smiles_line in zip(input_stream.GetOEMols(), smiles_file):
        smiles = smiles_line.split()[0]
        print(smiles)
        oechem.OEAddExplicitHydrogens(mol)

        # Retrieve params stored in molecule tag.
        params = json.loads(mol.GetStringData("SMIRNOFF_PARAMS"))
        t142_indices = params["t142"][0]  # Use first occurrence of t142.

        oemol = cmiles.utils.load_molecule(smiles)
        conformers = fragmenter.chemi.generate_conformers(oemol,
                                                          max_confs=1,
                                                          strict_stereo=False,
                                                          strict_types=False)
        try:
            cmiles_ids = cmiles.get_molecule_ids(smiles, strict=False)
        except ValueError:  # Stereochemistry errors
            continue
        mapped_smiles = cmiles_ids[
            'canonical_isomeric_explicit_hydrogen_mapped_smiles']
        print(mapped_smiles)
        qcarchive_mols = [
            cmiles.utils.mol_to_map_ordered_qcschema(conf, mapped_smiles)
            for conf in conformers.GetConfs()
        ]  # Only 1 mol
        job_idx = cmiles.utils.to_canonical_label(mapped_smiles, t142_indices)

        # Map torsion indices to canonical ordered mapped SMILES
        #  mol_copy = oechem.OEMol(oemol)
        #  oechem.OEAddExplicitHydrogens(mol_copy)
        #  cmiles._cmiles_oe.canonical_order_atoms(mol_copy)
        #  dih = []
        #  for m_idx in t142_indices:
        #      print(m_idx)
        #      atom = mol_copy.GetAtom(oechem.OEHasMapIdx(m_idx + 1))
        #      dih.append(atom.GetIdx())

        torsion_dict[job_idx] = {
            'initial_molecules': qcarchive_mols,
            'dihedral': [t142_indices],
            'grid': [15],
            'cmiles_identifiers': cmiles_ids
        }

    print(f"{len(torsion_dict)} molecules saved")

    # Save output.
    with open(OUTPUT_FILENAME, 'w') as json_file:
        json.dump(torsion_dict, json_file, indent=2)

    # Clean up.
    smiles_file.close()


if __name__ == "__main__":
    main()
