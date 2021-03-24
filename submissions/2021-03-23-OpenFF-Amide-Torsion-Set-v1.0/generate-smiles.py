"""This script requires the `constructure` package to run.
"""

from openff.toolkit.topology import Molecule
from openff.toolkit.utils import UndefinedStereochemistryError

from constructure.constructors import OpenEyeConstructor
from constructure.scaffolds import SCAFFOLDS
from constructure.substituents import SUBSTITUENTS
from constructure.utilities.openeye import smiles_to_image_grid


def main():

    amide_substituents = {
        # R group on the carboxyl carbon.
        1: [
            SUBSTITUENTS["hydrogen"],
            SUBSTITUENTS["ethyl"],
            "[R]OC",
            SUBSTITUENTS["vinyl"],
            "[R]c1ccc([O-])cc1",
            "[R]c1ccc(N(=O)(=O))cc1",
            SUBSTITUENTS["nitrile"],
            SUBSTITUENTS["nitro"],
        ],
        # R groups on the nitrogen.
        2: [
            SUBSTITUENTS["hydrogen"],
            SUBSTITUENTS["methyl"],
            "[R]c1ccc([O-])cc1",
            "[R]c1ccc(N(=O)(=O))cc1",
        ],
        3: [
            SUBSTITUENTS["hydrogen"],
            SUBSTITUENTS["methyl"],
        ],
    }

    full_smiles_set = set()

    # Enumerate different amide derivatives.
    for scaffold_name in [
        "carboxylic acid amide",
        "thiocarboxylic acid amide",
        "carboxylic acid amidine",
    ]:

        scaffold = SCAFFOLDS[scaffold_name]

        substituents = {**amide_substituents}

        if len(scaffold.r_groups) == 4:
            substituents[4] = [SUBSTITUENTS["methyl"]]

        derivatives = OpenEyeConstructor.enumerate_combinations(
            scaffold, substituents, validate=False
        )

        for derivative in derivatives:

            try:
                Molecule.from_smiles(derivative)
            except UndefinedStereochemistryError:

                derivative = (
                    Molecule.from_smiles(derivative, allow_undefined_stereo=True)
                    .enumerate_stereoisomers(undefined_only=True, max_isomers=1)[0]
                    .to_smiles()
                )

            full_smiles_set.add(derivative)

    # Enumerate different urea derivatives.
    for scaffold_name in ["urea", "thiourea"]:

        scaffold = SCAFFOLDS[scaffold_name]

        derivatives = OpenEyeConstructor.enumerate_combinations(
            scaffold,
            substituents={
                1: [
                    SUBSTITUENTS["hydrogen"],
                    SUBSTITUENTS["ethyl"],
                    SUBSTITUENTS["nitrile"],
                    "[R]S(=O)(=O)O",
                ],
                2: [
                    SUBSTITUENTS["hydrogen"],
                    SUBSTITUENTS["methyl"],
                ],
                3: [
                    SUBSTITUENTS["hydrogen"],
                    SUBSTITUENTS["methyl"],
                ],
                4: [
                    SUBSTITUENTS["hydrogen"],
                    SUBSTITUENTS["methyl"],
                ]
            },
            validate=False
        )

        full_smiles_set.update(derivatives)

    with open("dataset.smi", "w") as file:
        file.write("\n".join(sorted(full_smiles_set)))

    smiles_to_image_grid([*full_smiles_set], "dataset.pdf", cols=8)


if __name__ == "__main__":
    main()
