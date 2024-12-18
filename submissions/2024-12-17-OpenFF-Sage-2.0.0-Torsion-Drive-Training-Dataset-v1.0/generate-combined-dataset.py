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


def pull_record_id_cmiles(Opt: Type[TorsionDriveResultCollection]):
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
        # TODO: Check if updating dic would change the number of records
        rec_ids_cmiles.update(tmp_rec_ids_cmiles)

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

#with open('dataset_information.json') as f: dataset_information = json.load(f)
dataset_information = {
    "dataset_name": "OpenFF Sage 2.0.0 Torsion Drive Training Dataset v1.0",
    "dataset_tagline": "B3LYP-D3BJ/DZVP conformers applicable to drug-like molecules for OpenFF 2.0.0 Sage",
    "description": "A quantum chemical (QC) dataset curated to train the OpenFF 2.0.0 Sage torsion potentials. This QC dataset with the OpenFF default level of theory, B3LYP-D3BJ/DZVP, is used to benchmark Sage geometries and energetics. These optimized conformer geometries where used to train one dimensional torsional profiles. This Generation 2 dataset increases chemical diversity when compared to Generation 1, which are of value to our industry partners. Large molecules (>20 heavy atoms) were also included, including more flexible molecules and a greater degree of conformational variation which provide intramolecular interactions. This is the complete optimization dataset used for training OpenFF 2.0.0 Sage, consisting of the following datasets: 'OpenFF Gen 2 Torsion Set 1 Roche', 'OpenFF Gen 2 Torsion Set 2 Coverage', 'OpenFF Gen 2 Torsion Set 3 Pfizer Discrepancy', 'OpenFF Gen 2 Torsion Set 4 eMolecules  - Discrepancy', 'OpenFF Gen 2 Torsion Set 5 Bayer' and 'OpenFF Gen 2 Torsion Set 6 supplemental 2'. The `HydrogenBondFilter(method='baker-hubbard')` filter was applied, and the following record IDs were dropped due to issues with ForceBalance: 6098580, 2703504, 2703505, 18045478. Further information can be found in the curation scripts for the linked repositories.",
    "metadata.submitter": "Jennifer A Clark",
    "metadata.long_description_url": "https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/OpenFF-Sage-2.0.0-Torsion-Drive-Training-Dataset-v1.0"
}

dataset1 = TorsiondriveDataset(
    dataset_name=dataset_information["dataset_name"],
    dataset_tagline=dataset_information["dataset_tagline"],
    description=dataset_information["description"],
    provenance=provenance1,
)
dataset1.metadata.submitter = dataset_information["metadata.submitter"]
dataset1.metadata.long_description_url = dataset_information["metadata.long_description_url"]


# Have to add records this way to avoid a round trip through the toolkit.
# rec_and_mole is a linear representation of all the molecules
# records_by_cmiles groups together conformers

cmiles_count = {}
molecules = []
for record, molecule  in rec_and_mol:
    cmiles = rec_and_cmiles[record.id]

    if cmiles not in cmiles_count:
        molecules.append(molecules)
        cmiles_count[cmiles] = defaultdict(Counter)
    cmiles_count[cmiles][1][record.initial_molecules[0].get_hash()] += 1

    dataset1.add_molecule(
        index=rec_and_cmiles[record.id],
        molecule=None,
        extras=record.extras,
        keywords=record.specification.keywords,
        attributes=MoleculeAttributes.from_openff_molecule(molecule),
        initial_molecules=record.initial_molecules, 
        dihedrals=record.specification.keywords.dihedrals,
    )
    
# Note this: base_record.record_type

# Check that the molecules are identical
#opt_hashes = {rec.initial_molecules.get_hash() for rec, _mol in rec_and_mol}
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


# cmiles_count = Counter()

# _________ Pull Statistics from Dataset ____________

lx = len(cmiles_count)
n_confs1, n_heavy_atoms1, masses1, unique_charges1 = np.zeros(lx), [], np.zeros(lx), np.zeros(lx)
for i,cmiles, (mol, hashes) in enumerate(cmiles_count):
    n_confs1[i] = len(hashes)
    n_heavy_atoms1.append(mol.to_rdkit().GetNumHeavyAtoms())
    masses1[i] = sum([atom.mass.m for atom in mol.atoms])
    unique_charges1[i] = mol.total_charge.m_as(unit.elementary_charge)
unique_charges1 = sorted(set(unique_charges1))

#n_confs1 = [len(x) for (x, _) in cmiles_count] # the number of conformers and torsion drives for each
#n_heavy_atoms1 = np.array(
#    [mol.to_rdkit().GetNumHeavyAtoms() for mol in dataset1.molecules] # not weighted by conformer
#)
#masses1 = np.array(
#    [sum([atom.mass.m for atom in mol.atoms]) for mol in dataset1.molecules]
#)
#unique_charges1 = sorted(set(m.total_charge.m_as(unit.elementary_charge) for m in dataset1.molecules))

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
# With multiple torsions per unique molecule, n_molecules * confs.mean() no
# longer equals the number of conformers. instead, the number of dihedrals *
# confs.mean() should equal the number of conformers. The dataset contains one
# record per driven torsion (rather than combining multiple dihedrals into the
# same record), so n_records is the same as manually adding up len(dihedrals)
# for each record.
print("* Number of driven torsions: {}".format(dataset1.n_records))
print("* Number of conformers:", sum(n_confs1))
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
print("* Charges: {}".format(", ".join(unique_charges1)))
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
