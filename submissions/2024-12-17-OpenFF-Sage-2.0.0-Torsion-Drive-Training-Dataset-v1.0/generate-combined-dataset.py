import numpy as np
import json
import yaml
import requests
import logging
from typing import Type
from collections import Counter, defaultdict

logging.getLogger("openff").setLevel(logging.ERROR)

from openff.qcsubmit.results import TorsionDriveResultCollection
from openff.qcsubmit.datasets import TorsiondriveDataset

from openff.qcsubmit.factories import TorsiondriveDatasetFactory
from openff.qcsubmit.common_structures import MoleculeAttributes
from openff.qcsubmit.utils import _CachedPortalClient, portal_client_manager

from openff.toolkit.utils import OpenEyeToolkitWrapper, ToolkitRegistry
from openff.units import unit


def pull_record_id_cmiles(Opt: TorsionDriveResultCollection):
    """Pull CMILES strings associated with each molecule record

    Parameters
    ----------
    Opt
        A *ResultCollection object providing a reference to, and
        allows the retrieval of, data from a single result record
        stored in a QCFractal instance.
    
    """
    rec_ids_cmiles = {}
    for _, results in Opt.entries.items():
        tmp_rec_ids_cmiles = {result.record_id: result.cmiles for result in results}
        lx, ly = len(tmp_rec_ids_cmiles), len(rec_ids_cmiles)
        rec_ids_cmiles.update(tmp_rec_ids_cmiles)
        if len(rec_ids_cmiles) != lx + ly:
            raise ValueError("Multiple servers share record ids")

    return rec_ids_cmiles



# _________ Pull Dataset into TorsiondriveDataset object ____________
print("# Creating dataset")

file = requests.get(
    "https://raw.githubusercontent.com/openforcefield/openff-sage/37a36e7eeaf6cdca795847089a288bdff168c08a/data-set-curation/quantum-chemical/data-sets/1-2-0-td-set.json"
)
try:
    filtered_and_combined = TorsionDriveResultCollection.parse_raw(file.content)
except Exception:
    raise ImportError(f"Failed to import: {file.url}")

# Use cached portal client
client = _CachedPortalClient(
    "https://api.qcarchive.molssi.org:443", ".cache"
)
with portal_client_manager(lambda _: client):
    rec_and_mol = filtered_and_combined.to_records()
rec_and_cmiles = pull_record_id_cmiles(filtered_and_combined)

print("Number of results: ", filtered_and_combined.n_results, flush=True)
print("Finished converting to records", flush=True)

dataset_factory1 = TorsiondriveDatasetFactory()
provenance1 = dataset_factory1.provenance(ToolkitRegistry([OpenEyeToolkitWrapper]))

with open("ds_info.json") as f:
    dataset_information = json.load(f)

dataset1 = TorsiondriveDataset(
    dataset_name=dataset_information["dataset_name"],
    dataset_tagline=dataset_information["dataset_tagline"],
    description=dataset_information["description"],
    provenance=provenance1,
)
dataset1.metadata.submitter = dataset_information["metadata.submitter"]
dataset1.metadata.long_description_url = dataset_information["metadata.long_description_url"]

cmiles_count = defaultdict(Counter)
molecules = []
for record, molecule  in rec_and_mol:
    cmiles = rec_and_cmiles[record.id]

    if cmiles not in cmiles_count:
        molecules.append(molecule)
    hash = record.initial_molecules[0].get_hash()
    cmiles_count[cmiles][hash] += 1

    conformer_dihedral = "-{}-{}".format(
        list(cmiles_count[cmiles].keys()).index(hash),
        cmiles_count[cmiles][hash]
    )

    dataset1.add_molecule(
        index=rec_and_cmiles[record.id] + conformer_dihedral,
        molecule=None,
        extras=record.extras,
        keywords=record.specification.keywords,
        attributes=MoleculeAttributes.from_openff_molecule(molecule),
        initial_molecules=record.initial_molecules, 
        dihedrals=record.specification.keywords.dihedrals,
    )

# Check that the molecules are identical
opt_hashes = {
    mol.get_hash() 
    for rec, _mol in rec_and_mol
    for mol in rec.initial_molecules
}

new_hashes = {
    qcemol.identifiers.molecule_hash
    for moldata in dataset1.dataset.values()
    for qcemol in moldata.initial_molecules
}

print("Dataset creation preserved Molecules? ", opt_hashes == new_hashes)

# _________ Pull Statistics from Dataset ____________

lx = len(cmiles_count)
n_confs1, n_heavy_atoms1, masses1, unique_charges1 = np.zeros(lx), [], np.zeros(lx), np.zeros(lx)
for i, (cmiles, hashes) in enumerate(cmiles_count.items()):
    n_confs1[i] = len(hashes)
    n_heavy_atoms1.append(molecules[i].to_rdkit().GetNumHeavyAtoms())
    masses1[i] = sum([atom.mass.m for atom in molecules[i].atoms])
    unique_charges1[i] = molecules[i].total_charge.m_as(unit.elementary_charge)
unique_charges1 = sorted(set(unique_charges1))

elements1 = sorted(dataset1.metadata.dict()['elements'])


# _________ Write Output ____________

print("\n# Heavy Atom Counts")
counts1 = Counter(n_heavy_atoms1)
for n_heavy in sorted(counts1):
    print(f"{str(n_heavy):>3}: {counts1[n_heavy]}")

print("\n\n# Output for README\n")
print("* Description: {}".format(dataset1.description))
print("* Purpose: {}".format(dataset1.dataset_tagline))
print("* Name: {}".format(dataset1.dataset_name))
print("* Number of unique molecules: {}".format(len(cmiles_count)))
print("* Number of filtered molecules:", dataset1.n_filtered)
print("* Number of driven torsions: {}".format(dataset1.n_records))
print("* Number of conformers:", int(sum(n_confs1)))
print(
    "* Number of conformers (min, mean, max): {:.2f}, {:.2f}, {:.2f}".format(
        min(n_confs1), np.mean(n_confs1), max(n_confs1)
    )
)
print(
    "* Molecular weight (min, mean, max): {:.2f}, {:.2f}, {:.2f}".format(
        min(masses1), np.mean(masses1), max(masses1)
    )
)
print("* Charges: {}".format(", ".join([str(x) for x in unique_charges1])))
print("* Submitter: {}".format(dataset1.metadata.submitter))

print("\n## Metadata")
print(f"* Elements: {{{', '.join(elements1)}}}")

fields = ["basis", "implicit_solvent", "keywords", "maxiter", "method", "program"]
for spec, obj in dataset1.qc_specifications.items():
    od = obj.dict()
    print("* QC Specifications:", spec)
    for field in fields:
        print(f"  * {field}: {od[field]}")
    print("  * SCF Properties:")
    for field in od["scf_properties"]:
        print(f"    * {field}")

dataset1.export_dataset("dataset.json.bz2")
dataset1.molecules_to_file("dataset.smi", "smi")
dataset1.visualize("dataset.pdf", columns=8)
