from openff.qcsubmit.datasets import OptimizationDataset
from openff.qcsubmit.serializers import deserialize
from qcelemental.models.procedures import OptimizationProtocols

ds = OptimizationDataset.parse_obj(deserialize('dataset.json.bz2'))

ds.protocols = OptimizationProtocols(trajectory='final')
ds.export_dataset('dataset.json.bz2')
