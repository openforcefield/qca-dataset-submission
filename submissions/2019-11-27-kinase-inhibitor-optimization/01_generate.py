import json
import argparse
import gzip

from openeye import oechem
import cmiles
import fragmenter
from openmoltools import openeye


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Run psi4 to calculate bond orders')
    parser.add_argument('-i', '--infile', type=str,
                        help='Input molecule file')
    args = parser.parse_args()
    infile = args.infile

    oemols = fragmenter.chemi.file_to_oemols(infile)
    optimization_inputs = []

    optimization_jobs = 0
    for mol in oemols:

        name = mol.GetTitle()
        oechem.OEAddExplicitHydrogens(mol)
        smiles = cmiles.utils.mol_to_smiles(mol, mapped=False)
        cmiles_ids = cmiles.get_molecule_ids(smiles)
        mapped_smiles = cmiles_ids['canonical_isomeric_explicit_hydrogen_mapped_smiles']

        conformers = fragmenter.chemi.generate_conformers(mol)
        # Generate mol2 file for visualization and later AM1 WBO calculations
        openeye.molecule_to_mol2(conformers, tripos_mol2_filename='data/{}.mol2'.format(name), conformer=None)

        # Generate JSON geometry opt input
        qcschema_molecules = [cmiles.utils.mol_to_map_ordered_qcschema(conf, mapped_smiles) for conf in conformers.GetConfs()]
        optimization_jobs += len(qcschema_molecules)
        optimization_inputs.append({'initial_molecules': qcschema_molecules,
                                  'cmiles_identifiers': cmiles_ids})
    with gzip.open('optimization_inputs.json.gz', 'w') as f:
        f.write(json.dumps(optimization_inputs, indent=2, sort_keys=True).encode('utf-8'))

    print('Generated {} optimization jobs'.format(optimization_jobs))
