import tomllib
from pathlib import Path

import numpy as np
import qcportal  # noqa avoid zstd disaster
from openff.qcsubmit import workflow_components
from openff.qcsubmit.factories import OptimizationDatasetFactory
from openff.toolkit import ForceField, Molecule
from tqdm import tqdm

# Load config file and force field
ff = ForceField("openff-2.1.0.offxml")
with open("opt.toml", "rb") as f:
    config = tomllib.load(f)

# Load molecules
molecules = list()
with open("train.smi") as inp:
    for line in tqdm(inp, desc="Loading molecules"):
        if line.startswith("#"):
            continue
        pid, cmiles, *_ = line.split()
        mol = Molecule.from_smiles(cmiles, allow_undefined_stereo=True)
        labels = ff.label_molecules(mol.to_topology())[0]["ProperTorsions"]
        molecules.append(mol)

# initialize dataset factory
dataset_factory = OptimizationDatasetFactory()
dataset_factory.add_workflow_components(
    workflow_components.StandardConformerGenerator(
        max_conformers=10, rms_cutoff=0.5
    )
)

# populate dataset
dataset = dataset_factory.create_dataset(
    dataset_name=config["name"],
    tagline=config["name"],
    description=config["short_description"],
    molecules=molecules,
)
dataset.metadata.submitter = config["submitter"]
dataset.metadata.long_description_url = (
    "https://github.com/openforcefield/qca-dataset-submission/tree/master/"
    "submissions/" + str(Path.cwd().name)
)

# verify that the molecules are the same in input and output
old_smiles = {m.to_smiles(isomeric=False) for m in molecules}
new_smiles = {m.to_smiles(isomeric=False) for m in dataset.molecules}
assert not old_smiles.symmetric_difference(new_smiles)


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


# write output files
dataset.export_dataset("dataset.json.bz2")
dataset.molecules_to_file("output.smi", "smi")
dataset.visualize("dataset.pdf", columns=8)
