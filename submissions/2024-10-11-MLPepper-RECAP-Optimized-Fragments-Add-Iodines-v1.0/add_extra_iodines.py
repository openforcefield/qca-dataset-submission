from openff.qcsubmit.factories import OptimizationDatasetFactory
from openff.qcsubmit.datasets import load_dataset
from openff.qcsubmit.results import OptimizationResultCollection

from openff.qcsubmit import workflow_components
from openff.toolkit import Molecule
from qcportal import PortalClient

def main():
    
    ADDRESS = ""
    USERNAME = ""
    PASSWORD = ""

    client = PortalClient(address=ADDRESS, username=USERNAME, password=PASSWORD)
    
    existing_dataset = load_dataset('dataset.json.bz2')
    
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

    # fragment = workflow_components.RECAPFragmenter()
    # factory.add_workflow_components(fragment)

    stereo_expander = workflow_components.EnumerateStereoisomers()
    factory.add_workflow_components(stereo_expander)
    
    conf_gen = workflow_components.StandardConformerGenerator(max_conformers=5, toolkit="rdkit")
    factory.add_workflow_components(conf_gen)

    # should we remove molecules which don't have Boron and Silicon do we need the other fragments?
    
    new_iodines_dataset = factory.create_dataset(
        dataset_name="PubChem_I_fragments",
        description="A dataset of RECAP fragmented molecules from the SPICE2 Iodines replaced with Chlorines molecules extracted from PubChem.",
        tagline="PubChem Iodines molecule fragments from SPICE2.",
        molecules="extra_iodines.smi"
    )
    
    total_dataset = new_iodines_dataset + existing_dataset

    factory.export_workflow("opt_factory.json")
    total_dataset.visualize("dataset.pdf")
    total_dataset.export_dataset("opt_dataset_I_add_new.json.bz2")
    total_dataset.molecules_to_file("iodine_prepared_add_new.smi", "smi")
    
    #will this resubmit the complete ones?
    total_dataset.submit(client)

if __name__ == "__main__":
    main()