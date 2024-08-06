from openff.qcsubmit.results import OptimizationResultCollection
from openff.qcsubmit.results.filters import ConnectivityFilter, RecordStatusFilter, RecordStatusEnum
from openff.qcsubmit.common_structures import QCSpec, SCFProperties
from qcportal import PortalClient

DATASET_NAME  = ""
ADDRESS = ""
USERNAME = ""
PASSWORD = ""
client = PortalClient(address=ADDRESS, username=USERNAME, password=PASSWORD)

result_dataset = OptimizationResultCollection.from_server(client=client, datasets=DATASET_NAME, spec_name="aimnet2")

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
    dataset_name="ESP 50k opt",
    description="A combined dataset of the recap 50k molecules and Lily's Br subset each  molecule has 5 conformers which are optimised with AIMNET2 (wb97m-d3)",
    tagline="ESP DATASET",
    driver="energy",
    qc_specifications=[QCSpec(
        method="HF",
        basis="6-31G*",
        program="psi4",
        spec_name="HF/6-31G*-gas",
        spec_description="Standard HF gas",
    )]
)
basic_dataset.export_dataset("esp_50k_Br_singlepoint_dataset.json.bz2")
