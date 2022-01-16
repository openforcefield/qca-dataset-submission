"""A runnable script containing the minimal code from the README example."""
import itertools
import json
from collections import defaultdict
from typing import List, Dict

from openff.toolkit.topology import Molecule

from constructure.constructors import OpenEyeConstructor
from constructure.scaffolds import SCAFFOLDS, Scaffold
from constructure.substituents import SUBSTITUENTS
from constructure.utilities.openeye import smiles_to_image_grid


def enumerate_substituents(
    scaffold: Scaffold, max_substituents: Dict[str, int]
) -> List[Dict[int, List[str]]]:

    substituents_by_type = defaultdict(set)

    for substituent in max_substituents:

        substituent_type = OpenEyeConstructor.classify_substituent(substituent)
        substituents_by_type[substituent_type].add(substituent)

    all_r_group_substituents = {
        i: [
            substituent
            for scaffold_type in scaffold.r_groups[i]
            for substituent in substituents_by_type[scaffold_type]
        ]
        for i in scaffold.r_groups
    }

    combinations = itertools.product(
        *(all_r_group_substituents[i + 1] for i in range(len(all_r_group_substituents)))
    )

    group_combinations = [
        {i + 1: substituent for i, substituent in enumerate(combination)}
        for combination in combinations
    ]

    # Retain only combinations which have less than the max amount of each substituent
    # (except for hydrogen)
    enumerated_substituents = []

    for combination in group_combinations:

        substituent_counts = defaultdict(int)

        for substituent in combination.values():
            substituent_counts[substituent] += 1

        substituent_counts["[R][H]"] = 0

        if any(
            count > max_substituents[substituent]
            for substituent, count in substituent_counts.items()
        ):
            continue

        enumerated_substituents.append(
            {i: [smiles] for i, smiles in combination.items()}
        )

    return enumerated_substituents


def main():

    scaffolds = {
        "alkyl": [
            "alkyne",
            "alkene",
        ],
        "hetero": [
            "aldehyde",
            "ketone",
            "ketene",
            "alcohol",
            "enol",
            "carboxylic acid",
            "carboxylic acid ester",
            "carbonic acid monoester",
            "ether",
        ],
        "aryl": [
            "benzyl",
            "p-phenylene",
        ]
    }

    max_substituents = {
        "alkyl": {
            "hydrogen": 10,
            "methyl": 1,
            "ethyl": 1,
            "isopropyl": 1,
        },
        "hetero": {
            "hydrogen": 10,
            "methyl": 10,
            "ethyl": 10,
            "isopropyl": 1,
            "phenyl": 1,
            "benzyl": 1,
        },
        "aryl": {
            "hydrogen": 10,
            "methyl": 2,
            "aldehyde": 2,
            "hydroxyl": 2,
            "acid": 2,
            "propyne": 2
        }
    }

    # smiles = {
    #     smiles
    #     for group_type in scaffolds
    #     for scaffold in scaffolds[group_type]
    #     for smiles in OpenEyeConstructor.enumerate_combinations(
    #         SCAFFOLDS[scaffold],
    #         substituents={
    #             i: [
    #                 SUBSTITUENTS[substituent]
    #                 for substituent in substituents[group_type]
    #                 if OpenEyeConstructor.classify_substituent(SUBSTITUENTS[substituent])
    #                 in SCAFFOLDS[scaffold].r_groups[i]
    #             ]
    #             for i in SCAFFOLDS[scaffold].r_groups
    #         },
    #     )
    # }

    enumerated_smiles = {
        smiles
        for group_type in scaffolds
        for scaffold in scaffolds[group_type]
        for substituent in enumerate_substituents(
            SCAFFOLDS[scaffold],
            {
                SUBSTITUENTS[substituent]: max_count
                for substituent, max_count in max_substituents[group_type].items()
            }
        )
        for smiles in OpenEyeConstructor.enumerate_combinations(
            SCAFFOLDS[scaffold],
            substituents=substituent
        )
    }

    # Select a random stereoisomer.
    stereo_smiles = set()

    for smiles in enumerated_smiles:

        molecule = Molecule.from_smiles(smiles, allow_undefined_stereo=True)

        # Limit to a single ring per molecule
        if molecule.n_rings > 1:
            continue

        stereoisomers = [
            molecule, *molecule.enumerate_stereoisomers(True, max_isomers=1)
        ]

        stereo_smiles.add(stereoisomers[-1].to_smiles())

    print(len(stereo_smiles))

    smiles_to_image_grid([*stereo_smiles], "bcc-c-o-h-v2.png", cols=8)

    with open("molecules.smi", "w") as file:
        file.write("\n".join(sorted(stereo_smiles)))


if __name__ == "__main__":
    main()
