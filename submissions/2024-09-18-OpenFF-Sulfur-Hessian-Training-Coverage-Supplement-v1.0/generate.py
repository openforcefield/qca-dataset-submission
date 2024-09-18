from pathlib import Path

import numpy as np
import qcportal  # noqa avoid zstd disaster
from openff.qcsubmit.factories import BasicDatasetFactory
from openff.qcsubmit.results import OptimizationResultCollection
from openff.qcsubmit.results.filters import (
    RecordStatusEnum,
    RecordStatusFilter,
)
from openff.qcsubmit.utils import _CachedPortalClient, portal_client_manager
from openff.toolkit import ForceField
from qcaide import Submission
from qcportal.singlepoint import SinglepointDriver

# Load config file and force field
ff = ForceField("openff-2.1.0.offxml")
config = Submission.from_toml("opt.toml")

client = _CachedPortalClient("https://api.qcarchive.molssi.org", ".")
opt = OptimizationResultCollection.from_server(
    client,
    datasets=["OpenFF Sulfur Optimization Training Coverage Supplement v1.0"],
)

print(f"Retrieved {opt.n_results} results")

with portal_client_manager(lambda _: client):
    recs = opt.filter(
        RecordStatusFilter(status=RecordStatusEnum.complete)
    ).to_records()

print(f"Filtered to {len(recs)} completed records")

# these are the final_molecules from the optimization records
molecules = [mol for _rec, mol in recs]

# initialize dataset factory
dataset_factory = BasicDatasetFactory()

# populate dataset
dataset = dataset_factory.create_dataset(
    dataset_name=config.name,
    molecules=molecules,
    description=config.description,
    tagline=config.name,
)
dataset.driver = SinglepointDriver.hessian
dataset.metadata.submitter = config.submitter
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
