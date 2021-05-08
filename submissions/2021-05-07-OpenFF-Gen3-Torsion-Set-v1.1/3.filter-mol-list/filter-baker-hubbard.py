from multiprocessing import Pool
from typing import Tuple

import mdtraj
import numpy
from simtk import unit

from tqdm import tqdm

from openff.toolkit.topology import Molecule


def smiles_to_molecule(smiles: str) -> Molecule:
    """Attempts to parse a smiles pattern into a molecule object, guessing the
    stereochemistry if it is undefined.

    Parameters
    ----------
    smiles
        The smiles pattern to parse.

    Returns
    -------
    The parsed molecule.
    """
    from openff.toolkit.topology import Molecule
    from openff.toolkit.utils import UndefinedStereochemistryError

    try:
        molecule = Molecule.from_smiles(smiles)
    except UndefinedStereochemistryError:

        molecule: Molecule = Molecule.from_smiles(smiles, allow_undefined_stereo=True)

        stereoisomers = molecule.enumerate_stereoisomers(
            undefined_only=True, max_isomers=1
        )

        if len(stereoisomers) > 0:
            # We would ideally raise an exception here if the number of stereoisomers
            # is zero, however due to the way that the OFF toolkit perceives pyramidal
            # nitrogen stereocenters these would show up as undefined stereochemistry
            # but have no enumerated stereoisomers.
            molecule = stereoisomers[0]

    return molecule


def apply_h_bond_filter(smiles: str) -> Tuple[str, bool, bool]:

    try:

        molecule = smiles_to_molecule(smiles)
        molecule.generate_conformers(n_conformers=800, rms_cutoff=0.1 * unit.angstrom)

        if len(molecule.conformers) == 0:
            return smiles, False, True

        conformers = numpy.array(
            [
                conformer.value_in_unit(unit.nanometers).tolist()
                for conformer in molecule.conformers
            ]
        )

        topology = molecule.to_topology()._to_mdtraj()
        trajectory = mdtraj.Trajectory(conformers * unit.nanometers, topology)

        h_bonds = mdtraj.baker_hubbard(trajectory, freq=0.0, periodic=False)

        return smiles, len(h_bonds) == 0, False

    except:

        return smiles, False, True


def main():
    import argparse
    parser = argparse.ArgumentParser('Filter internal H bonding molecules')
    parser.add_argument('-i', '--input_smi', default='generated_mols_out.smi')
    parser.add_argument('-n', '--nprocesses', type=int, default=4)
    args = parser.parse_args()

    import sys
    print('# '+' '.join(sys.argv))
    
    import time
    print("Starting H-bonding molecule search...")
    start = time.time()

    
    with open(args.input_smi) as file:

        smiles_to_consider = {
            line for line in file.read().split("\n") if len(line) > 0
        }

    with Pool(processes=args.nprocesses) as pool:

        smiles_keep_error = list(
            tqdm(
                pool.imap(apply_h_bond_filter, smiles_to_consider),
                total=len(smiles_to_consider)
            )
        )

    smiles_to_keep = {
        smiles
        for smiles, keep, error in smiles_keep_error
        if keep and not error
    }
    print(f'# smiles to keep: {len(smiles_to_keep)}')
    smiles_to_remove = {
        smiles
        for smiles, keep, error in smiles_keep_error
        if not keep and not error
    }
    smiles_with_errors = {
        smiles
        for smiles, keep, error in smiles_keep_error
        if error
    }

    with open("smiles-to-keep.smi", "w") as file:
        file.write("\n".join(smiles_to_keep))
    with open("smiles-to-remove.smi", "w") as file:
        file.write("\n".join(smiles_to_remove))
    with open("smiles-with-errors.smi", "w") as file:
        file.write("\n".join(smiles_with_errors))

    end = time.time()
    print("Elapsed time %.2f seconds" % (end-start))


if __name__ == '__main__':
    main()
