from openff.qcsubmit.results import OptimizationResultCollection
from openff.qcsubmit.results.filters import ConnectivityFilter, RecordStatusFilter, RecordStatusEnum
from openff.qcsubmit.common_structures import QCSpec, SCFProperties
# from openff.qcsubmit.common_structures import DDXSettings, SCFProperties, Metadata
from qcportal import PortalClient

DATASET_NAME  = ""
ADDRESS = ""
USERNAME = ""
PASSWORD = ""
client = PortalClient(address=ADDRESS, username=USERNAME, password=PASSWORD)

result_dataset = OptimizationResultCollection.from_server(client=client, datasets=DATASET_NAME, spec_name="aimnet2-wb97m-d3")
print(result_dataset)
print(f"Dataset has {result_dataset.n_results} raw results.")

status_filter = RecordStatusFilter(status=RecordStatusEnum.complete)
result_dataset = result_dataset.filter(status_filter)

print(f"Dataset has {result_dataset.n_results} complete records.")

# remove connectivity issues
connect_filter = ConnectivityFilter()
result_dataset = result_dataset.filter(connect_filter)

print(f"Dataset has {result_dataset.n_results} after removing connection issues")
# mock some data for the dataset will be removed later

print(f"Building single point dataset")
basic_dataset = result_dataset.create_basic_dataset(
    dataset_name="ESP 50k opt Iodines",
    description="A combined dataset of the recap 50k molecules filtered for Cls and Brs which are then converted to I",
    tagline="ESP DATASET",
    driver="energy",
    qc_specifications=[QCSpec(
        method="wb97x-d",
        basis="def2-tzvpp",
        program="psi4",
        spec_name="wb97x-d/def2-tzvpp",
        spec_description="wb97x-d/def2-tzvpp gas",
    )]
)

basic_dataset.export_dataset("iodine_filtered.json.bz2")
basic_dataset.molecules_to_file("iodine_filtered.smi", "smi")
