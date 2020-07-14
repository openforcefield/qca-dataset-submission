"""
Script to run validation using qcsubmit on dataset.json files
"""

from qcsubmit.datasets import BasicDataset, OptimizationDataset, TorsiondriveDataset
from qcsubmit.serializers import deserialize



def create_dataset(file_name):
    """
    Create a dataset from the input json.
    """
    data = deserialize(file_name)



# def validate