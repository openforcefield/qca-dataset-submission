import numpy as np
import requests
import logging
from typing import Type
from collections import Counter

logging.getLogger("openff").setLevel(logging.ERROR)

from openff.qcsubmit.results import OptimizationResultCollection
from openff.qcsubmit.datasets import OptimizationDataset
from openff.qcsubmit.factories import OptimizationDatasetFactory
from openff.qcsubmit.common_structures import MoleculeAttributes

from openff.toolkit.utils import OpenEyeToolkitWrapper, ToolkitRegistry
from openff.units import unit


def pull_record_id_cmiles(Opt):
    """Pull CMILES strings associated with each molecule record

    Parameters
    ----------
    Opt
        An OptimizationResultCollection object providing a reference to, and
        allows the retrieval of, data from a single optimization result record
        stored in a QCFractal instance.
    
    """
    rec_ids_cmiles = {}
    for _, results in Opt.entries.items():
        tmp_rec_ids_cmiles = {result.record_id: result.cmiles for result in results}
        # TODO: Check if updating dic would change the number of records
        rec_ids_cmiles.update(tmp_rec_ids_cmiles)

    return rec_ids_cmiles



# _________ Pull Dataset into OptimizationDataset object ____________
print("# Creating dataset")

file = requests.get(
    "https://raw.githubusercontent.com/openforcefield/openff-sage/37a36e7eeaf6cdca795847089a288bdff168c08a/data-set-curation/quantum-chemical/data-sets/1-2-0-opt-set-v3.json"
)
filtered_and_combined = OptimizationResultCollection.parse_raw(file.content)
rec_and_mol = filtered_and_combined.to_records()
rec_and_cmiles = pull_record_id_cmiles(filtered_and_combined)

print("Number of results: ", filtered_and_combined.n_results, flush=True)
print("Finished converting to records", flush=True)

dataset_factory1 = OptimizationDatasetFactory()
provenance1 = dataset_factory1.provenance(ToolkitRegistry([OpenEyeToolkitWrapper]))

dataset1 = OptimizationDataset(
    dataset_name="OpenFF Sage 2.0.0 Training Optimization Dataset v1.0",
    dataset_tagline="B3LYP-D3BJ/DZVP conformers applicable to drug-like molecules for OpenFF 2.0.0 Sage",
    description=(
        "A quantum chemical (QC) dataset curated to train OpenFF 2.0.0 Sage, "
        "with reparametrized Lennard-Jones (LJ) and valence parameters, the latter "
        "relevant to this dataset. This QC dataset with the OpenFF default level of "
        "theory, B3LYP-D3BJ/DZVP, is used to benchmark Sage geometries and energetics. "
        "These optimized conformer geometries where used in conjunction with the QC "
        "dataset used to train one dimensional torsional profiles. This Generation 2 "
        "dataset increases chemical diversity when compared to Generation 1, which are "
        "of value to our industry partners. Large molecules (>20 heavy atoms) were also "
        "included, including more flexible molecules and a greater degree of conformational "
        "variation which provide intramolecular interactions.\nThis is the complete optimization "
        "dataset used for training OpenFF 2.0.0 Sage, consisting of the following datasets: "
        "'OpenFF Gen 2 Opt Set 1 Roche', 'OpenFF Gen 2 Opt Set 2 Coverage', 'OpenFF Gen 2 Opt Set 3 Pfizer Discrepancy', "
        "'OpenFF Gen 2 Opt Set 4 eMolecules  - Discrepancy', and 'OpenFF Gen 2 Opt Set 5 Bayer'.\nThe "
        "following filters were applied: RecordStatusFilter(status=RecordStatusEnum.complete), "
        "ConnectivityFilter(tolerance=1.2), UndefinedStereoFilter(), ConformerRMSDFilter(max_conformers=10), "
        "and ElementFilter(allowed_elements=['H', 'C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I'])"
        "Further information can be found in the curation scripts for the linked repositories."
    ),
    provenance=provenance1,
)
dataset1.metadata.submitter = "Jennifer A Clark"
dataset1.metadata.long_description_url = (
    "https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/"
    "2024-12-12-OpenFF-Sage-2.0.0-Training-Optimization-Dataset-v1.0"
)


# Have to add records this way to avoid a round trip through the toolkit.
# rec_and_mole is a linear representation of all the molecules
# records_by_cmiles groups together conformers
records_by_cmiles = {}
for record, molecule in rec_and_mol:
    cmiles = rec_and_cmiles[record.id]
    if cmiles in records_by_cmiles:
        records_by_cmiles[cmiles].append((record, molecule))
    else:
        records_by_cmiles[cmiles] = [(record, molecule)]

for cmiles, records in records_by_cmiles.items():
    base_record, base_molecule = records[0]
    base_molecule._conformers = [m.conformers[0] for _, m in records]

    # Need to input initial molecules so that QCA will recognize the inputs as the same.
    initial_mols = [rec.initial_molecule for rec, _ in records]

    dataset1.add_molecule(
        index=cmiles,
        molecule=None,
        initial_molecules=initial_mols, 
        attributes=MoleculeAttributes.from_openff_molecule(base_molecule),
        extras=base_record.extras,
        keywords=base_record.specification.keywords,
    )

# Check that the molecules are identical
opt_hashes = {rec.initial_molecule.get_hash() for rec, _mol in rec_and_mol}

new_hashes = {
    qcemol.identifiers.molecule_hash
    for moldata in dataset1.dataset.values()
    for qcemol in moldata.initial_molecules
}

print("Dataset creation preserved Molecules? ", opt_hashes == new_hashes)



# _________ Pull Statistics from Dataset ____________

n_confs1 = np.array([mol.n_conformers for mol in dataset1.molecules])

n_heavy_atoms1 = np.array(
    [mol.to_rdkit().GetNumHeavyAtoms() for mol in dataset1.molecules]
)

masses1 = np.array(
    [sum([atom.mass.m for atom in mol.atoms]) for mol in dataset1.molecules]
)

unique_charges1 = [
    str(charge)
    for charge in sorted(
        set(
            [
                mol.total_charge.m_as(unit.elementary_charge)
                for mol in dataset1.molecules
            ]
        )
    )
]

# _________ Write Output ____________

print("\n# Heavy Atom Counts")
counts1 = Counter(n_heavy_atoms1)
for n_heavy in sorted(counts1):
    print(f"{str(n_heavy):>3}: {counts1[n_heavy]}")

print("\n\n# Output for README\n")
print("* Description: {}".format(dataset1.description))
print("* Purpose: {}".format(dataset1.dataset_tagline))
print("* Name: {}".format(dataset1.dataset_name))
print("* Number of unique molecules: {}".format(dataset1.n_molecules))
print("* Number of filtered molecules:", dataset1.n_filtered)
print("* Number of conformers: {}".format(dataset1.n_records))
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
print(f"* Elements: {{{', '.join(sorted(dataset1.metadata.dict()['elements']))}}}")

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
