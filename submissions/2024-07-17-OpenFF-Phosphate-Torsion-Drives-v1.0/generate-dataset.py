import re
import tomllib
from collections import defaultdict
from pathlib import Path

import numpy as np
import qcportal  # noqa avoid zstd disaster
from openff.qcsubmit import workflow_components
from openff.qcsubmit.factories import TorsiondriveDatasetFactory
from openff.qcsubmit.utils import get_symmetry_classes, get_symmetry_group
from openff.qcsubmit.workflow_components import TorsionIndexer
from openff.toolkit import ForceField, Molecule
from tqdm import tqdm

# use Sage 2.1.0 for labeling
ff = ForceField("openff-2.1.0.offxml")

# load config
with open("td.toml", "rb") as f:
    config = tomllib.load(f)

junk = re.compile("[(,)]")

# proper torsion parameter IDs to drive
target_params = [
    "t2",
    "t3",
    "t4",
    "t45",
    "t46",
    "t50",
    "t95",
    "t123a",
    "t159",
    "t160",
]

# drive each torsion matching one of the target parameters for each molecule
molecules = list()
with open("input.smi") as inp:
    for line in tqdm(inp, desc="Tagging torsions"):
        smiles = line.strip()
        mol = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
        labels = ff.label_molecules(mol.to_topology())[0]["ProperTorsions"]
        ti = TorsionIndexer()
        sym_classes = get_symmetry_classes(mol)
        for tors, p in labels.items():
            sym = get_symmetry_group(tors[1:3], sym_classes)
            if p.id in target_params:
                ti.add_torsion(tors, sym, (-165, 180))
        mol.properties["dihedrals"] = ti
        molecules.append(mol)


# build the dataset
dataset_factory = TorsiondriveDatasetFactory()
dataset_factory.add_workflow_components(
    workflow_components.StandardConformerGenerator(max_conformers=10)
)

dataset = dataset_factory.create_dataset(
    dataset_name=config["name"],
    tagline=config["name"],
    description=config["short_description"],
    molecules=molecules,
)

dataset.metadata.submitter = "ntBre"
dataset.metadata.long_description_url = (
    "https://github.com/openforcefield/qca-dataset-submission/tree/master/"
    "submissions/" + str(Path.cwd().name)
)


# check that the molecules in the dataset match the initial molecules
old_smiles = {m.to_smiles(isomeric=False) for m in molecules}
new_smiles = {m.to_smiles(isomeric=False) for m in dataset.molecules}
# usually I check the symmetric difference, but in this case, three SMILES had
# no torsions matching the target parameters and were filtered out
assert not new_smiles - old_smiles
assert len(old_smiles) - 3 == len(new_smiles)


# summarize dataset for readme
confs = np.array([len(mol.conformers) for mol in dataset.molecules])

print("* Number of unique molecules:", dataset.n_molecules)
# With multiple torsions per unique molecule, n_molecules * confs.mean() no
# longer equals the number of conformers. instead, the number of dihedrals *
# confs.mean() should equal the number of conformers. The dataset contains one
# record per driven torsion (rather than combining multiple dihedrals into the
# same record), so n_records is the same as manually adding up len(dihedrals)
# for each record.
print("* Number of driven torsions:", dataset.n_records)
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


# export the dataset
dataset.export_dataset("dataset.json.bz2")
dataset.molecules_to_file("output.smi", "smi")
dataset.visualize("dataset.pdf", columns=8)
