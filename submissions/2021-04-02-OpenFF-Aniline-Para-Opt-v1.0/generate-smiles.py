from constructure.constructors import OpenEyeConstructor
from constructure.scaffolds import Scaffold
from constructure.substituents import SUBSTITUENTS
from constructure.utilities.openeye import smiles_to_image_grid


def main():

    substituents = {
        # R group on the carboxyl carbon.
        1: [
            "[R][S-]",
            "[R]N(C)(C)",
            "[R][O-]",
            "[R]N(C)",
            "[R]N",
            "[R]Nc1ccccc1",
            "[R]O",
            "[R]NO",
            "[R]OC",
            "[R]NC(=O)N",
            "[R]C(C)(C)(C)",
            "[R]C",
            "[R]C=C",
            "[R]Oc1ccccc1",
            "[R]c1ccccc1",
            "[R]SC",
            "[R]NC(=O)C",
            "[R]F",
            "[R]Sc1ccccc1",
            "[R]N=[N+]=[N-]",
            "[R]S",
            "[R]N=C=O",
            "[R]Cl",
            "[R]Br",
            "[R]OC(=O)C",
            "[R]OC(F)(F)(F)",
            "[R]OS(=O)(=O)C",
            "[R]N=C=S",
            "[R]OC(=O)C(F)(F)(F)",
            "[R]OC#N",
            "[R][N+]",
            "[R][N+](C)(C)(C)",
            "[R]C#C",
            "[R]C(Br)(Br)(Br)",
            "[R]C(=O)N",
            "[R]CS(=O)(=O)[O-]",
            "[R]C(=O)N(C)",
            "[R]C=O",
            "[R]C(=O)O",
            "[R]C(=O)OC",
            "[R]C(=O)C",
            "[R]C(Cl)(Cl)(Cl)",
            "[R]C(F)(F)(F)",
            "[R]N=C",
            "[R]S-C#N",
            "[R]C(=O)Cl",
            "[R]S(=O)(=O)N",
            "[R]C#N",
            "[R]N(=O)(=O)",
            "[R]S(=O)(=O)C(F)(F)(F)",
            "[R][N+]#N",
        ],
        # R groups on the nitrogen.
        2: [
            SUBSTITUENTS["hydrogen"],
            # SUBSTITUENTS["methyl"],
        ],
        3: [
            SUBSTITUENTS["hydrogen"],
            # SUBSTITUENTS["methyl"],
        ],
    }

    scaffold = Scaffold(
        smiles="N([R2])([R3])c1ccc([R1])cc1",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["hydrogen", "alkyl", "aryl"]
        }
    )

    derivatives = OpenEyeConstructor.enumerate_combinations(
        scaffold, substituents, validate=False
    )

    smiles_to_image_grid([*derivatives], "dataset.pdf", cols=8)

    with open("dataset.smi", "w") as file:
        file.write("\n".join(sorted(derivatives)))


if __name__ == '__main__':
    main()
