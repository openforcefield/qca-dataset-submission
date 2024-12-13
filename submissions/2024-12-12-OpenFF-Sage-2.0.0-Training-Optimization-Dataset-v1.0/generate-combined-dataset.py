
from openff.qcsubmit.results import OptimizationResultCollection
from openff.qcsubmit.datasets import OptimizationDataset
from openff.qcsubmit.factories import OptimizationDatasetFactory
from openff.qcsubmit.common_structures import MoleculeAttributes
import numpy as np
from openff.toolkit.utils import OpenEyeToolkitWrapper, ToolkitRegistry
from openff.units import unit
import logging
logging.getLogger("openff").setLevel(logging.ERROR)

filtered_and_combined = OptimizationResultCollection.parse_file(
    "/Users/jennifer.clark/bin/openff-sage/data-set-curation/quantum-chemical/data-sets/1-2-0-opt-set-v3.json",
#    "/Users/jennifer.clark/bin/openff-sage/data-set-curation/quantum-chemical/data-sets/1-2-0-td-set.json",
)
rec_and_mol = filtered_and_combined.to_records()

print('Number of results: ',filtered_and_combined.n_results,flush=True)
print('Finished converting to records',flush = True)

# Need to input initial molecules so that QCA will recognize the inputs as the same.
initial_mols = [rec[0].initial_molecule for rec in rec_and_mol]

dataset_factory1 = OptimizationDatasetFactory()
provenance1 = dataset_factory1.provenance(ToolkitRegistry([OpenEyeToolkitWrapper]))

dataset1 = OptimizationDataset(
    dataset_name="OpenFF Sage 2.0.0 Training Optimization v1.0",
    dataset_tagline="B3LYP-D3BJ/DZVP conformers applicable to drug-like molecules for OpenFF 2.0.0 Sage",
    description=("A quantum chemical (QC) dataset curated to train OpenFF 2.0.0 Sage, "
    "with reparametrized Lennard-Jones (LJ) and valence parameters, the latter "
    "relevent to this dataset. This QC dataset with the OpenFF default level of "
    "theory, B3LYP-D3BJ/DZVP, is used to benchmark Sage geometries and energetics. "
    "These optimized conformer geometries where used in conjunction with the QC "
    "dataset used to train one dimensional torsional profiles. This Generation 2 "
    "dataset increases chemical diversity when compared to Generation 1, which are "
    "of value to our industry partners. Large molecules (>20 heavy atoms) were also "
    "included, including more flexible molecules and a greater degree of conformational "
    "variation which provide intramolecular interactions."),
    provenance=provenance1,
)
dataset1.metadata.submitter = "Jennifer A Clark"
dataset1.metadata.long_description_url = (
        "https://github.com/openforcefield/qca-dataset-submission/tree/master/"
        "submissions/"
        "2024-12-12-OpenFF-Sage-2.0.0-Training-Optimization-Dataset-v1.0"
    )


# Have to add records this way to avoid a round trip through the toolkit.
records_by_cmiles= {}
for record, molecule in rec_and_mol:
    cmiles = molecule.to_smiles(isomeric=True, explicit_hydrogens=True, mapped=True)
    if cmiles in records_by_cmiles.keys():
        records_by_cmiles[cmiles].append((record, molecule))
    else:
        records_by_cmiles[cmiles]=[(record, molecule)]

for records in records_by_cmiles.values():
    base_record, base_molecule = records[0]
    base_molecule._conformers = [m.conformers[0] for _, m in records]

    dataset1.add_molecule(
        index=base_molecule.to_smiles(
            isomeric=True, explicit_hydrogens=False, mapped=False
        ),
        molecule=None,
        initial_molecules=[rec.initial_molecule for rec, _ in records],
        attributes=MoleculeAttributes.from_openff_molecule(base_molecule),
        extras=base_record.extras,
        keywords=base_record.specification.keywords,
    )

# Check that the molecules are identical
opt_hashes = {
        rec.initial_molecule.get_hash() for rec, _mol in rec_and_mol
    }

new_hashes = {
    qcemol.identifiers.molecule_hash
    for moldata in dataset1.dataset.values()
    for qcemol in moldata.initial_molecules
}

print('Molecules are the same? ',opt_hashes==new_hashes)

n_confs1 = np.array(
    [mol.n_conformers for mol in dataset1.molecules]
)

n_heavy_atoms1 = np.array(
    [mol.to_rdkit().GetNumHeavyAtoms() for mol in dataset1.molecules]
)

masses1 = np.array([
    sum([atom.mass.m for atom in mol.atoms])
    for mol in dataset1.molecules
])

elements1 = set(
    atom.symbol
    for mol in dataset1.molecules
    for atom in mol.atoms
)

unique_charges1 = [str(charge) for charge in sorted(set([
    mol.total_charge.m_as(unit.elementary_charge)
    for mol in dataset1.molecules
]))]

from collections import Counter

print("# heavy atoms")
counts1 = Counter(n_heavy_atoms1)
for n_heavy in sorted(counts1):
    print(f"{str(n_heavy):>3}: {counts1[n_heavy]}")

###### Write Output
print('* Description: {}'.format(dataset1.description))
print('* Purpose: {}'.format(dataset1.dataset_tagline))
print('* Name: {}'.format(dataset1.dataset_name))
print('* Number of unique molecules: {}'.format(dataset1.n_molecules))
print("* Number of filtered molecules:", dataset1.n_filtered)
print('* Number of conformers: {}'.format(dataset1.n_records))
print('* Number of conformers (min, mean, max): {:.2f}, {:.2f}, {:.2f}'.format(min(n_confs1),np.mean(n_confs1),max(n_confs1)))
print('* Molecular weight (min, mean, max): {:.2f}, {:.2f}, {:.2f}'.format(min(masses1),np.mean(masses1),max(masses1)))
print('* Charges: {}'.format(' '.join(unique_charges1)))
print('* Submitter: {}'.format(dataset1.metadata.submitter))

print("## Metadata")
print(f"* Elements: {{{', '.join(dataset1.metadata.dict()['elements'])}}}")

def print_field(od, field): print(f"  * {field}: {od[field]}")

fields = ["basis", "implicit_solvent", "keywords", "maxiter", "method", "program"]
for spec, obj in dataset1.qc_specifications.items():
    od = obj.dict()
    print("* QC Specifications:", spec)
    for field in fields:
        print_field(od, field)
    print("  * SCF Properties:")
    for field in od["scf_properties"]:
        print(f"    * {field}")

dataset1.export_dataset("dataset.json.bz2")
dataset1.molecules_to_file('dataset.smi', 'smi')
dataset1.visualize("dataset.pdf", columns=8)
