import os
from tempfile import NamedTemporaryFile

import openeye

from openeye import oechem
from openff.recharge.esp.storage import MoleculeESPStore, MoleculeESPRecord
from openff.recharge.utilities.geometry import reorder_conformer
from openforcefield.topology import Molecule
from qcelemental.models.results import WavefunctionProtocolEnum

from qcsubmit.common_structures import QCSpec, PCMSettings
from qcsubmit.factories import BasicDatasetFactory


def molecule_from_record(record: MoleculeESPRecord) -> Molecule:
    """Converts an ``openff-recharge`` ESP record to to an Open Force Field
    molecule."""

    oe_molecule = oechem.OEMol()
    oechem.OESmilesToMol(oe_molecule, record.tagged_smiles)
    ordered_conformer = reorder_conformer(oe_molecule, record.conformer)

    # Clear the records index map.
    for atom in oe_molecule.GetAtoms():
        atom.SetMapIdx(0)

    oe_molecule.DeleteConfs()
    oe_molecule.NewConf(oechem.OEFloatArray(ordered_conformer.flatten()))

    with NamedTemporaryFile(suffix=".mol2") as file:

        # Workaround for stereochemistry being incorrectly perceived.
        molecule = Molecule.from_openeye(oe_molecule, allow_undefined_stereo=True)

        molecule.to_file(file.name, "mol2")
        molecule = molecule.from_file(file.name)

    return molecule


def main():

    # Load in the set of molecules to add to the data set.
    esp_store = MoleculeESPStore()

    molecules = [
        molecule_from_record(esp_record)
        for smiles in esp_store.list()
        for esp_record in esp_store.retrieve(smiles)
    ]

    # Store the conformers as SDF files.
    os.makedirs("conformers", exist_ok=True)

    for i, molecule in enumerate(molecules):
        molecule.to_file(os.path.join("conformers", f"{i}.sdf"), "SDF")

    # Generate the data set to submit.
    factory = BasicDatasetFactory(
        qc_specifications={
            "resp-2-vacuum": QCSpec(
                method="pw6b95",
                basis="aug-cc-pV(D+d)Z",
                spec_name="resp-2-vacuum",
                spec_description=(
                    "The quantum chemistry specification used in the RESP2 publication "
                    "for the vacuum (i.e. no PCM) calculations."
                ),
                store_wavefunction=WavefunctionProtocolEnum.orbitals_and_eigenvalues
            ),
            "resp-2-water": QCSpec(
                method="pw6b95",
                basis="aug-cc-pV(D+d)Z",
                spec_name="resp-2-water",
                spec_description=(
                    "The quantum chemistry specification used in the RESP2 publication "
                    "for the aqueous (i.e. with PCM) calculations."
                ),
                store_wavefunction=WavefunctionProtocolEnum.orbitals_and_eigenvalues,
                implicit_solvent=PCMSettings(
                    units="angstrom",
                    cavity_Type="GePol",
                    cavity_Area=0.3,
                    cavity_Scaling=True,
                    cavity_RadiiSet="Bondi",
                    cavity_Mode="Implicit",
                    medium_SolverType="CPCM",
                    medium_Solvent="Water",
                )
            ),
        }
    )

    data_set = factory.create_dataset(
        dataset_name="OpenFF BCC Refit Study COH v1.0",
        molecules=molecules,
        description="A data set curated for the initial stage of the on-going OpenFF "
        "study which aims to co-optimize the AM1BCC bond charge correction (BCC) "
        "parameters against an experimental training set of density and enthalpy of "
        "mixing data points and a QM training set of electric field data."
        "\n\n"
        "The initial data set is limited to only molecules composed of C, O, H. This "
        "limited scope significantly reduces the number of BCC parameters which must "
        "be retrained, thus allowing for easier convergence of the initial "
        "optimizations."
        "\n\n"
        "The included molecules are those included in the experimental data set as "
        "well as an additional set chosen to ensure that each BCC parameter to train "
        "has been sufficiently (at least five instances) represented and exercised."
        "\n\n"
        "The conformers included in the set where generated using version 0.0.1a4 of "
        "the openff-recharge package. The exact conformer generation settings are "
        "attached as provenance.",
        tagline="C,H,O single point training data for BCC refits.",
    )

    # Validate that the data set matches expectations.
    assert data_set.n_molecules == 94
    assert data_set.n_records == 215

    # Attach the conformer generation provenance
    data_set.provenance["openff-recharge"] = "0.0.1a4"
    data_set.provenance["conformer-generation"] = (
        '{"method": "omega-elf10", "sampling_mode": "dense", "max_conformers": 5}'
    )
    data_set.provenance["openeye"] = openeye.__version__

    # Correct the dataset metadata.
    data_set.metadata.submitter = "simonboothroyd"
    data_set.metadata.long_description_url = (
        "https://github.com/openforcefield/qca-dataset-submission/tree/master/"
        "submissions/"
        "2020-10-30-OpenFF-BCC-Refit-Study-COH"
    )

    # Export the data set.
    data_set.export_dataset("dataset.json.xz")
    data_set.molecules_to_file('molecules.smi', 'smi')

    data_set.visualize("dataset.pdf", columns=8)


if __name__ == '__main__':
    main()
