import re
import sys
from collections import namedtuple
from pathlib import Path

import numpy as np
import qcportal  # noqa avoid zstd disaster
from openff.qcsubmit import workflow_components
from openff.qcsubmit.common_structures import QCSpec
from openff.qcsubmit.factories import OptimizationDatasetFactory
from openff.toolkit import Molecule
from qcaide import Submission


def load_ccd(filename: str) -> list[str]:
    "Load OpenEye SMILES from the CCD CIF"

    blocks = 0
    in_block = False
    smiles = list()
    with open(filename) as inp:
        for line in inp:
            if re.match(r"^data_", line):
                blocks += 1
                in_block = True
            sp = line.split()
            # previously taking first smiles arbitrarily, now looking for
            # OpenEye smiles because the cactvs ones are weird
            if (
                in_block
                and len(sp) > 5
                and "SMILES" == sp[1]
                and "OpenEye" in sp[2]
            ):
                in_block = False
                # SMILES in the sixth position, some are wrapped in quotes
                smiles.append(sp[5].removeprefix('"').removesuffix('"'))

    print(f"found {len(smiles)}/{blocks} smiles")
    return smiles


def try_from_smiles(s: str) -> Molecule | None:
    "Try to build a Molecule from SMILES, return None on failure"
    if "." in s:
        print(f"period in smiles: `{s}`", file=sys.stderr)
        return None
    try:
        mol = Molecule.from_smiles(s, allow_undefined_stereo=True)
    except Exception as e:
        print(e, file=sys.stderr)
        return None

    return mol


def elements(mol: Molecule) -> set[int]:
    return {atom.atomic_number for atom in mol.atoms}


Mol = namedtuple("M", ["smiles", "natoms", "charge", "inchi", "elements"])

with open("inchis.dat") as inp:
    inchis = {line.strip() for line in inp}

# retrieved and unzipped from
# https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz
smiles = load_ccd("components.cif")

mols = list()
failed = 0
for s in smiles:
    if (mol := try_from_smiles(s)) is None:
        failed += 1
        continue
    mols.append(
        Mol(
            charge=mol.total_charge,
            natoms=mol.n_atoms,
            smiles=mol.to_smiles(),
            inchi=mol.to_inchikey(),
            elements=elements(mol),
        )
    )


print(f"Converting {failed}/{len(smiles)} molecules failed")

weird_smiles = [
    "B1234B567B189B212B33%10B454B656B787B911C232B4%105C7612",
]

# from the strategy document
target_metals = {
    # at minimum: Pd, Fe, Zn, Mg, Cu, Li
    46,
    26,
    30,
    12,
    29,
    3,
    # also desirable: Rh, Ir, Pt, Ni, Cr, Ag
    45,
    77,
    78,
    28,
    24,
    47,
}

# Also include organic elements and halides: C, H, P, S, O, N, F, Cl, Br, but
# skipping explicit search for these for now to focus on the metals. most if
# not all of these are already covered by our existing datasets

dsmols = list()
for m in sorted(
    filter(
        lambda m: m.natoms > 10
        and abs(m.charge.magnitude) < 4.0
        and "." not in m.smiles
        and m.smiles not in weird_smiles
        and m.inchi not in inchis
        and len(m.elements & target_metals) > 0,
        mols,
    ),
    key=lambda m: m.natoms,
)[:100]:
    dsmols.append(Molecule.from_smiles(m.smiles, allow_undefined_stereo=True))

print(len(dsmols))

# Generate dataset
config = Submission.from_toml("opt.toml")
# we want to store more properties and info from the wavefunction, but that
# will come later in a separate singlepoint dataset
dataset_factory = OptimizationDatasetFactory(
    qc_specifications={
        "BP86/def2-TZVP": QCSpec(
            method="BP86",
            basis="def2-TZVP",
            spec_name="BP86/def2-TZVP",
        )
    }
)
dataset_factory.add_workflow_components(
    workflow_components.StandardConformerGenerator(
        max_conformers=10, rms_cutoff=0.5
    )
)
dataset = dataset_factory.create_dataset(
    dataset_name=config.name,
    tagline=config.name,
    description=config.description,
    molecules=dsmols,
    processors=6,
)
print(dataset.n_molecules, dataset.n_records)
dataset.metadata.submitter = config.submitter
dataset.metadata.long_description_url = (
    "https://github.com/openforcefield/qca-dataset-submission/tree/master/"
    "submissions/" + str(Path.cwd().name)
)

# verify that no new molecules appear in the output; can't use symmetric
# difference like I usually do because some molecules (~32) do get filtered
# out, I think from conformer generation failure
old_smiles = {m.to_smiles(isomeric=False) for m in dsmols}
new_smiles = {m.to_smiles(isomeric=False) for m in dataset.molecules}
assert not (new_smiles - old_smiles)


# summarize dataset for readme
confs = np.array([len(mol.conformers) for mol in dataset.molecules])

print("* Number of unique molecules:", dataset.n_molecules)
print("* Number of filtered molecules:", dataset.n_filtered)
print("* Number of conformers:", sum(confs))
print(
    "* Number of conformers per molecule (min, mean, max): "
    f"{confs.min()}, {confs.mean():.2f}, {confs.max()}"
)

masses = [
    [
        sum([atom.mass.m for atom in molecule.atoms])
        for molecule in dataset.molecules
    ]
]
print(f"* Mean molecular weight: {np.mean(np.array(masses)):.2f}")
print(f"* Max molecular weight: {np.max(np.array(masses)):.2f}")
print("* Charges:", sorted(set(m.total_charge.m for m in dataset.molecules)))

print("## Metadata")
print(f"* Elements: {{{', '.join(dataset.metadata.dict()['elements'])}}}")


def print_field(od, field):
    print(f"\t* {field}: {od[field]}")


fields = [
    "basis",
    "implicit_solvent",
    "keywords",
    "maxiter",
    "method",
    "program",
]
for spec, obj in dataset.qc_specifications.items():
    od = obj.dict()
    print("* Spec:", spec)
    for field in fields:
        print_field(od, field)
    print("\t* SCF properties:")
    for field in od["scf_properties"]:
        print(f"\t\t* {field}")

dataset.export_dataset("dataset.json.bz2")
dataset.molecules_to_file("output.smi", "smi")
dataset.visualize("dataset.pdf", columns=8)
