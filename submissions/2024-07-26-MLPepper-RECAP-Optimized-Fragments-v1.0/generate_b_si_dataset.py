from openff.qcsubmit.factories import OptimizationDatasetFactory
from openff.qcsubmit import workflow_components
import h5py
from openff.toolkit import Molecule

# read the molecules from the hdf5 file and remove all stereo info as we will expand the states anyway
def main():
    f = h5py.File("pubchem-boron-silicon.hdf5", "r")
    print("Extracting molecules")
    with open("pubchem-boron-silicon.smi", "w") as output:
        for name, entry in f.items():
            mapped_smiles = entry["smiles"][0].decode("utf-8")
            molecule = Molecule.from_mapped_smiles(mapped_smiles=mapped_smiles, allow_undefined_stereo=True)
            output.write(f"{molecule.to_smiles(isomeric=False, mapped=False, explicit_hydrogens=False)} {name}\n")

    f.close()
            


    factory = OptimizationDatasetFactory()
    # setup the aimnet2 specs and tags
    factory.clear_qcspecs()
    factory.add_qc_spec(
        method="wb97m-d3",
        basis=None,
        program="aimnet2",
        spec_name="aimnet2-wb97m-d3",
        spec_description="aiment2 wb97m-d3 opt spec"
    )
    factory.compute_tag = "aimnet2"
    # build the workflow

    fragment = workflow_components.RECAPFragmenter()
    factory.add_workflow_components(fragment)

    stereo_expander = workflow_components.EnumerateStereoisomers()
    factory.add_workflow_components(stereo_expander)
    
    conf_gen = workflow_components.StandardConformerGenerator(max_conformers=5, toolkit="rdkit")
    factory.add_workflow_components(conf_gen)

    # should we remove molecules which don't have Boron and Silicon do we need the other fragments?

    dataset = factory.create_dataset(
        dataset_name="PubChem_B-Si_fragments",
        description="A dataset of RECAP fragmented molecules from the SPICE2 Boron and Silicon molecules extracted from PubChem.",
        tagline="PubChem Boron and Silicon molecule fragments from SPICE2.",
        molecules="pubchem-boron-silicon.smi"
    )
    # loop through and remove molecules without B or Si
    new_dataset = {}
    for name, entry in dataset.dataset.items():
        off_mol = entry.get_off_molecule(include_conformers=False)
        atomic_numbers = set([atom.atomic_number for atom in off_mol.atoms])
        if 5 in atomic_numbers or 14 in atomic_numbers:
            new_dataset[name] = entry
        else: 
            continue

    dataset.dataset = new_dataset
    factory.export_workflow("opt_factory.json")
    dataset.visualize("dataset.pdf")
    dataset.export_dataset("opt_dataset_B_Si.json.bz2")
    dataset.molecules_to_file("dataset.smi", "smi")

if __name__ == "__main__":
    main()
