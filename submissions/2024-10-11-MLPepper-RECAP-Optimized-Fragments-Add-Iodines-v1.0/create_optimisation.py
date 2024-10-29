from openff.qcsubmit.factories import OptimizationDatasetFactory
from openff.qcsubmit import workflow_components
from openff.toolkit import Molecule
from qcportal import PortalClient

def main():
    
    #create iodines from mlpepper iodines. 
    iodine_smiles = []
    smiles = []
    with open("./dataset_mlpepper.smi", "r") as input:
        for mols in input:
            smiles.append(mols.strip('\n'))

    for molecule in smiles:
        if 'Cl' in molecule:
            iodine_smiles.append(molecule.replace('Cl','I'))
        if 'Br' in molecule:
            iodine_smiles.append(molecule.replace('Br','I'))
            
    with open('./iodine_dataset.smi','w+') as file:
        for iodine in iodine_smiles:
            file.write(iodine + '\n')

    client = PortalClient(address="", username="", password="")
    
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

    stereo_expander = workflow_components.EnumerateStereoisomers()
    factory.add_workflow_components(stereo_expander)
    
    conf_gen = workflow_components.StandardConformerGenerator(max_conformers=5, toolkit="rdkit")
    factory.add_workflow_components(conf_gen)
        
    dataset = factory.create_dataset(
        dataset_name="PubChem_I_fragments",
        description="A dataset of RECAP fragmented molecules from the SPICE2 Iodines replaced with Chlorines molecules extracted from PubChem.",
        tagline="PubChem Iodines molecule fragments from SPICE2.",
        molecules="iodine_smiles.smi"
    )
    
    # loop through and remove molecules without I
    new_dataset = {}
    for name, entry in dataset.dataset.items():
        off_mol = entry.get_off_molecule(include_conformers=False)
        atomic_numbers = set([atom.atomic_number for atom in off_mol.atoms])
        if 53 in atomic_numbers:
            new_dataset[name] = entry
        else: 
            continue
    
    dataset.dataset = new_dataset
    factory.export_workflow("opt_factory.json")
    dataset.visualize("opt_dataset.pdf")
    dataset.export_dataset("opt_dataset_I.json.bz2")
    dataset.molecules_to_file("iodine_prepared.smi", "smi")
    
    dataset.submit(client)
    
if __name__ == "__main__":
    main()