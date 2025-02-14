import numpy as np
import json
import requests
import date
from collections import Counter, defaultdict

from deepdiff import DeepDiff
import periodictable

from qcportal import PortalClient
from qcportal.optimization import OptimizationDatasetEntry
from qcportal.torsiondrive import TorsiondriveDatasetEntry
DatasetEntry = {"optimization": OptimizationDatasetEntry, "torsiondrive": TorsiondriveDatasetEntry}

from qcfractal.snowflake import FractalSnowflake
snowflake = FractalSnowflake()
client = snowflake.client()
#ADDRESS = "https://api.qcarchive.molssi.org:443/"
#client = PortalClient(ADDRESS, cache_dir=".")


# _________ Pull Record IDs of Relevant Datasets ____________
print("Getting record ids")

file = requests.get(
    "https://raw.githubusercontent.com/openforcefield/openff-sage/37a36e7eeaf6cdca795847089a288bdff168c08a/data-set-curation/quantum-chemical/data-sets/1-2-0-opt-set-v3.json"
)
data = json.loads(file.content)
provenance = data["provenance"]
# list with: {type, record_id, cmiles, inchi_key}
entry_dicts = data["entries"][ADDRESS]
dataset_type = entry_dicts[0]["type"]


# _________ Get Record IDs / Make Entries ____________
print("Getting records")
records = client.get_records([int(x["record_id"]) for x in entry_dicts], missing_ok=False)


# ________ Get / Check Specification ______________
print("Retrieving and comparing specifications")

spec = records[0].specification
diff_dict = Counter()
for record in records:
    diff_dict[str(DeepDiff(spec, record.specification))] += 1
    
if "{}" not in diff_dict and len(diff_dict) > 1:
    raise ValueError("Record specifications are not all the same!")


# ________ Get Molecule Entries ______________

cmiles_by_record_id = {
    int(x["record_id"]): {"cmiles": x["cmiles"], "mol": None} 
    for x in entry_dicts
}
entries = []
for record in records:
    entries.append(DatasetEntry[dataset_type](
        name=cmiles_by_record_id[record.id]["cmiles"],
        initial_molecule=record.initial_molecule,
    ))
    cmiles_by_record_id[record.id]["mol"] = record.initial_molecule


# _________ Initialize New Dataset ____________
print("Initializing new dataset")

with open("ds_info.json") as f:
    dataset_information = json.load(f)

dataset = client.add_dataset(
    dataset_type,
    dataset_name=dataset_information["dataset_name"],
    tagline=dataset_information["dataset_tagline"],
    description=dataset_information["description"],
    provenance=provenance,
    default_tag="openff",
    owner_user="openffbot",
    metadata={
        "submitter": dataset_information["metadata.submitter"],
        "creation_data": date.today(),
        'collection_type': 'OptimizationDataset',
        'long_description_url': dataset_information["metadata.long_description_url"],
        "short description": dataset_information["dataset_tagline"],
        "dataset_name": dataset_information["dataset_name"],
        "elements": provenance['applied-filters']['ElementFilter-3']['allowed_elements'],
    },
)


# ________ Add Specification to Dataset ______________
print("Adding specification to dataset")

dataset.add_specification(
    name='default', 
    specification=spec
)


# ________ Add Entries to Dataset ______________
print("Adding entries to dataset")

dataset.add_entries(entries)


# _________ Pull Statistics from Dataset ____________
print("Generating Molecular Statistics")

cmiles_count = defaultdict(Counter)
molecules = []
for recid, x in cmiles_by_record_id.items():
    cmiles = x["cmiles"]

    if cmiles not in cmiles_count:
        molecules.append(x["mol"])
    hash = x["mol"].get_hash()
    cmiles_count[cmiles][hash] += 1

#    conformer_dihedral = "-{}-{}".format(
#        list(cmiles_count[cmiles].keys()).index(hash),
#        cmiles_count[cmiles][hash]
#    )

lx = len(cmiles_count)
n_confs, n_heavy_atoms, masses, unique_charges = np.zeros(lx), [], np.zeros(lx), np.zeros(lx)
elements = []
for i, (cmiles, hashes) in enumerate(cmiles_count.items()):
    n_confs[i] = len(hashes)
    n_heavy_atoms.append(len([x for x in molecules[0].symbols if x != "H"]))
    elements.extend(list(set([x for x in molecules[0].symbols])))
    masses[i] = sum([getattr(periodictable, x).mass for x in molecules[0].symbols])
    unique_charges[i] = molecules[i].molecular_charge
    
unique_charges = sorted(set(unique_charges))

elements = sorted(list(set(elements)))


# _________ Write Output ____________

print("\n# Heavy Atom Counts")
counts1 = Counter(n_heavy_atoms)
for n_heavy in sorted(counts1):
    print(f"{str(n_heavy):>3}: {counts1[n_heavy]}")

print("\n\n# Output for README\n")
print("* Description: {}".format(dataset.description))
print("* Purpose: {}".format(dataset.dataset_tagline))
print("* Name: {}".format(dataset.dataset_name))
print("* Number of unique molecules: {}".format(len(cmiles_count)))
print("* Number of filtered molecules:", dataset.n_filtered)
print("* Number of driven torsions: {}".format(dataset.n_records))
print("* Number of conformers:", int(sum(n_confs)))
print(
    "* Number of conformers (min, mean, max): {:.2f}, {:.2f}, {:.2f}".format(
        min(n_confs), np.mean(n_confs), max(n_confs)
    )
)
print(
    "* Molecular weight (min, mean, max): {:.2f}, {:.2f}, {:.2f}".format(
        min(masses), np.mean(masses), max(masses)
    )
)
print("* Charges: {}".format(", ".join([str(x) for x in unique_charges])))
print("* Submitter: {}".format(dataset.metadata.submitter))

print("\n## Metadata")
print(f"* Elements: {{{', '.join(elements)}}}")

fields = ["basis", "implicit_solvent", "keywords", "maxiter", "method", "program"]
for spec, obj in dataset.qc_specifications.items():
    od = obj.dict()
    print("* QC Specifications:", spec)
    for field in fields:
        print(f"  * {field}: {od[field]}")
    print("  * SCF Properties:")
    for field in od["scf_properties"]:
        print(f"    * {field}")

# Do I need to output these?
dataset.export_dataset("dataset.json.bz2")
dataset.molecules_to_file("dataset.smi", "smi")
dataset.visualize("dataset.pdf", columns=8)
