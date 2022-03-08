from collections import defaultdict

from openff.qcsubmit.common_structures import Metadata, QCSpec, MoleculeAttributes
from openff.qcsubmit.datasets import BasicDataset
from openff.qcsubmit.results import OptimizationResultCollection
from openff.qcsubmit.results.filters import ConnectivityFilter
from qcelemental.models import DriverEnum
from qcelemental.models.results import WavefunctionProtocolEnum
from qcportal import FractalClient


def main():

    optimization_results = OptimizationResultCollection.from_server(
        client=FractalClient(),
        datasets=["OpenFF ESP Industry Benchmark Set v1.0"],
        spec_name="HF/6-31G*"
    )
    filtered_records = optimization_results.filter(
        ConnectivityFilter()
    ).to_records()

    dataset = BasicDataset(
        dataset_name="OpenFF ESP Industry Benchmark Set v1.1",
        dataset_tagline="HF/6-31G* conformers of public industry benchmark molecules.",
        description="A dataset containing molecules from the `OpenFF Industry Benchmark "
        "Season 1 v1.1` with conformers taken from the "
        "`OpenFF ESP Industry Benchmark Set v1.0` optimization (HF/6-31G*) set and with "
        "wavefunctions stored and any conformers that underwent a connectivity change "
        "filtered out.",
        driver=DriverEnum.energy,
        metadata=Metadata(
            submitter="simonboothroyd",
            long_description_url=(
                "https://github.com/openforcefield/qca-dataset-submission/tree/master/"
                "submissions/"
                "2022-02-02-OpenFF-ESP-Industry-Benchmark-Set-v1.1-single-point"
            )
        ),
        qc_specs={
            "HF/6-31G*": QCSpec(
                program="psi4",
                method="hf",
                basis="6-31G*",
                spec_name="HF/6-31G*",
                spec_description=(
                    "The standard HF/6-31G* basis used to derive RESP style charges."
                ),
                store_wavefunction=WavefunctionProtocolEnum.orbitals_and_eigenvalues
            )
        }
    )

    records_by_cmiles = defaultdict(list)

    for record, molecule in filtered_records:

        records_by_cmiles[
            molecule.to_smiles(isomeric=True, explicit_hydrogens=True, mapped=True)
        ].append((record, molecule))

    for records in records_by_cmiles.values():

        base_record, base_molecule = records[0]
        base_molecule._conformers = [m.conformers[0] for _, m in records]

        dataset.add_molecule(
            index=base_molecule.to_smiles(
                isomeric=True, explicit_hydrogens=False, mapped=False
            ),
            molecule=base_molecule,
            attributes=MoleculeAttributes.from_openff_molecule(base_molecule),
            extras=base_record.extras,
            keywords=record.keywords,
        )

    # Export the data set.
    dataset.export_dataset("dataset.json.xz")
    dataset.molecules_to_file('dataset.smi', 'smi')

    dataset.visualize("dataset.pdf", columns=8)

    print(dataset.qc_specifications)

    print("N MOLECULES", dataset.n_molecules)
    print("N RECORDS", dataset.n_records)


if __name__ == '__main__':
    main()
