from qcportal import PortalClient
from qcelemental.models.results import WavefunctionProtocolEnum
from openff.qcsubmit.results import OptimizationResultCollection,BasicResultCollection
from openff.qcsubmit.datasets import OptimizationDataset
from openff.qcsubmit.factories import OptimizationDatasetFactory
from openff.qcsubmit.common_structures import Metadata, QCSpec, MoleculeAttributes
from collections import defaultdict
import numpy as np
from openff.toolkit import Molecule
import tqdm
from openff.toolkit.utils import OpenEyeToolkitWrapper, ToolkitRegistry
from openff.units import unit
import itertools
import multiprocess
from openeye import oechem
import logging
# # suppress stereochemistry warnings
logging.getLogger("openff").setLevel(logging.ERROR)

from openeye import oechem

filtered_and_combined = OptimizationResultCollection.parse_file('filtered_and_combined_nagl2_opt.json')
print('Number of results: ',filtered_and_combined.n_results,flush=True)

rec_and_mol = filtered_and_combined.to_records()
print('Finished converting to records',flush = True)

# Need to input initial molecules so that QCA will recognize the inputs as the same.
initial_mols = [rec[0].initial_molecule for rec in rec_and_mol]

dataset_factory1 = OptimizationDatasetFactory()
provenance1 = dataset_factory1.provenance(ToolkitRegistry([OpenEyeToolkitWrapper]))

dataset1 = OptimizationDataset(
    dataset_name="OpenFF NAGL2 Training Optimization Dataset v4.0",
    dataset_tagline="B3LYP-D3BJ/DZVP conformers of diverse fragment molecules.",
    description=('Description will be updated later'), # will add at the end
    provenance=provenance1
)
dataset1.metadata.submitter = "amcisaac"
dataset1.metadata.long_description_url = (
        "https://github.com/openforcefield/qca-dataset-submission/tree/master/"
        "submissions/"
        "2024-11-27-OpenFF-NAGL2-Training-Optimization-Dataset-v4.0"
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

print('* Name: {}'.format(dataset1.dataset_name))
print('* Number of unique molecules: {}'.format(dataset1.n_molecules))
print('* Number of conformers: {}'.format(dataset1.n_records))
print('* Number of conformers (min, mean, max): {:.2f}, {:.2f}, {:.2f}'.format(min(n_confs1),np.mean(n_confs1),max(n_confs1)))
print('* Molecular weight (min, mean, max): {:.2f}, {:.2f}, {:.2f}'.format(min(masses1),np.mean(masses1),max(masses1)))
print('* Charges: {}'.format(' '.join(unique_charges1)))

print("## Metadata")
print(f"* Elements: {{{', '.join(dataset1.metadata.dict()['elements'])}}}")

def print_field(od, field): print(f"  * {field}: {od[field]}")

fields = ["basis", "implicit_solvent", "keywords", "maxiter", "method", "program"]
for spec, obj in dataset1.qc_specifications.items():
    od = obj.dict()
    print("* Spec:", spec)
    for field in fields:
        print_field(od, field)
    print("  * SCF properties:")
    for field in od["scf_properties"]:
        print(f"    * {field}")

dataset1.metadata.long_description=(("A dataset containing molecules from the "
        "`MLPepper RECAP Optimized Fragments v1.0` and `MLPepper RECAP Optimized Fragments v1.0 Add Iodines` "
        "with additional conformers and optimized at the OpenFF default level of theory (B3LYP-D3BJ/DZVP). "
        "The dataset is intended to be used for calculating single point energies and properties, "
        "which will then be used to train our second-generation graph neural network charge model (NAGL2)."
        "This is the final dataset, with part 1 and part 2 combined, and with errored records, geometric rearrangements, and stereochemistry issues filtered out.\n\n"
        "For each molecule, a set of up to 5 conformers were generated by:\n"
        "  * generating a set of up to 1000 conformers with a RMS cutoff of 0.1 Å "
        "using the OpenEye backend of the OpenFF toolkit\n"
        "  * applying ELF conformer selection (max 5 conformers) using OpenEye\n\n"
        "Dataset information:\n"
        "* Number of unique molecules: {}\n"
        "* Number of conformers: {}\n"
        "* Number of conformers (min, mean, max): {:.2f}, {:.2f}, {:.2f}\n"
        "* Molecular weight (min, mean, max): {:.2f}, {:.2f}, {:.2f}\n"
        "* Charges: {}".format(dataset1.n_molecules,dataset1.n_records,min(n_confs1),np.mean(n_confs1),max(n_confs1),min(masses1),np.mean(masses1),max(masses1),' '.join(unique_charges1))
        ))


dataset1.export_dataset("dataset.json.bz2")
dataset1.molecules_to_file('dataset.smi', 'smi')
dataset1.visualize("dataset.pdf", columns=8)
