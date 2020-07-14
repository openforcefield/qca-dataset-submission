"""
Script to run validation using qcsubmit on dataset.json files
"""

from qcsubmit.datasets import BasicDataset, OptimizationDataset, TorsiondriveDataset
from qcsubmit.serializers import deserialize
import json
import sys

datasets = {"BasicDataset": BasicDataset, "OptimizationDataset": OptimizationDataset, "TorsiondriveDataset": TorsiondriveDataset}


def validate_dataset(file_name):
    """
    Create a dataset from the input json.
    """
    print(file_name)
    data = deserialize(file_name)
    dataset_type = data["dataset_type"]
    # now try and make the correct dataset
    dataset_class = datasets.get(dataset_type, None)
    if dataset_class is not None:
        # validate the dataset
        dataset = dataset_class.parse_file(file_name)
        dataset._get_missing_basis_coverage(raise_errors=True)

    else:
        raise RuntimeError(f"The dataset type {dataset_type} is not supported.")


if __name__ == "__main__":
    file_names = json.loads(sys.argv[0])
    dataset_paths = []
    for file in file_names:
        if "dataset.json" in file:
            dataset_paths.append(file)

    # validate all datasets
    for dataset in dataset_paths:
        validate_dataset(dataset)