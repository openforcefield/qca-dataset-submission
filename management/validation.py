"""
Script to run validation using qcsubmit on dataset.json files
"""

import copy
import json
import os
from argparse import ArgumentParser


import pandas as pd
from github import Github
from qcsubmit.datasets import (BasicDataset, OptimizationDataset,
                               TorsiondriveDataset)
from qcsubmit.exceptions import (DatasetInputError, DihedralConnectionError,
                                 LinearTorsionError, MissingBasisCoverageError,
                                 MolecularComplexError)
from qcsubmit.serializers import deserialize

datasets = {
    "BasicDataset": BasicDataset,
    "OptimizationDataset": OptimizationDataset,
    "TorsiondriveDataset": TorsiondriveDataset}

check_mark = ":heavy_check_mark:"
cross = ":x:"
missing = ":heavy_exclamation_mark:"
REPO_NAME = 'openforcefield/qca-dataset-submission'

def get_data(file_name):
    """
    Return the deserialized dataset file.
    """
    return deserialize(file_name)


def create_dataset(dataset_data):
    dataset_type = dataset_data["dataset_type"]
    dataset_class = datasets.get(dataset_type, None)
    if dataset_class is not None:
        return dataset_class.parse_obj(dataset_data)
    else:
        raise RuntimeError(f"The dataset type {dataset_type} is not supported.")


def validate_dataset(dataset_data):
    """
    Create a dataset from the data and run normal validation on each molecule.
    Catch each of the error types and report them.
    """
    errors = {
        "cmiles": [],
        "dihedrals": [],
        "linear": [],
        "complex": []
    }
    data_copy = copy.deepcopy(dataset_data)
    # remove the entries so they can be checked one by one
    entries = data_copy.pop("dataset")
    # remove the scf props and meta data as this will be checked in a different step
    del data_copy["scf_properties"]
    del data_copy["metadata"]
    dataset = create_dataset(data_copy)

    # now check each entry
    for entry in entries.values():
        try:
            dataset.add_molecule(**entry, molecule=None)
        except DatasetInputError:
            # this mean the cmiles is not valid
            errors["cmiles"].append(entry.index)
        except DihedralConnectionError:
            # the torsion is not connected
            errors["dihedrals"].append(entry.index)
        except LinearTorsionError:
            errors["linear"].append(entry.index)
        except MolecularComplexError:
            errors["complex"].append(entry.index)

    return errors


def check_metadata(dataset_data):
    # remove the dataset and scf props
    data_copy = copy.deepcopy(dataset_data)
    del data_copy["scf_properties"]
    del data_copy["dataset"]
    dataset = create_dataset(data_copy)
    try:
        dataset.metadata.validate_metadata(raise_errors=True)
        return check_mark
    except DatasetInputError:
        return cross


def check_scf_props(dataset_data):
    # remove the metadata and dataset
    data_copy = copy.deepcopy(dataset_data)
    del data_copy["dataset"]
    del data_copy["metadata"]
    try:
        _ = create_dataset(data_copy)
        return check_mark
    except DatasetInputError:
        return cross


def check_basis_coverage(dataset_data):
    "Note this will only work if the full dataset can be loaded."

    try:
        dataset = create_dataset(dataset_data)
    except Exception:
        return missing
    try:
        dataset._get_missing_basis_coverage(raise_errors=True)
        return check_mark
    except MissingBasisCoverageError:
        return cross


def get_meta_info(dataset_data):
    return [dataset_data.get("dataset_name", missing),
            dataset_data.get("dataset_type", missing),
            dataset_data.get("method", missing),
            dataset_data.get("basis", missing)]


def main_validation(dataset_names):
    """
    Generate a report dataframe for each dataset found.
    """
    dataset_dataframe = {}
    # these are the attributes checked
    index = ["**Dataset Name**", "**Dataset Type**", "**Method**", "**Basis**", "**Valid Cmiles**",
               "**Connected Dihedrals**", "**No Linear Torsions**", "**No Molecular Complexes**", "**Complete Metatdata**",
               "**Valid SCF Properties**", "**Full Basis Coverage**"]
    _error_order = ["cmiles", "dihedrals", "linear", "complex"]

    for dataset_name in dataset_names:
        dataset_validators = []
        # get the data from the dataset
        data = get_data(dataset_name)
        # get the metadata
        dataset_validators.extend(get_meta_info(data))
        # check the first set of entry errors
        errors = validate_dataset(data)
        for error in _error_order:
            if errors[error]:
                dataset_validators.append(cross)
            else:
                dataset_validators.append(check_mark)
        # now check the metadata
        dataset_validators.append(check_metadata(data))
        # now check the scf
        dataset_validators.append(check_scf_props(data))
        # now check the basis
        dataset_validators.append(check_basis_coverage(data))

        dataset_dataframe[dataset_name] = dataset_validators
    # now make the dataframe
    df = pd.DataFrame(data=dataset_dataframe, index=index)
    comment = f"""
    ## QCSubmit Validation Report
    
    {df.to_markdown()}
    """

    # postprocess due to raw spacing above
    comment = "\n".join([substr.strip() for substr in comment.split('\n')])

    return comment


def main():

    parser = ArgumentParser(description="Validation methods for QCSubmit datasets.")
    parser.add_argument("dataset_files", help="This is the dataset file that should be validated.")
    parser.add_argument("pull_number", type=int, help="This is the PR number that danger bot will report on.")

    args = parser.parse_args()


    # now work out what is to be validated
    file_names = json.loads(args.dataset_files)
    dataset_paths = []
    for file in file_names:
        if "dataset.json" in file:
            dataset_paths.append(file)

    comment = main_validation(dataset_paths)

    # now we need the pr and to add the comment.
    g = Github(os.environ['GH_TOKEN'])
    repo = g.get_repo(REPO_NAME)
    pr = repo.get_pull(args.pull_number)
    pr.create_issue_comment(comment)

    # now we need to work out if the workflow should fail
    if missing in comment or cross in comment:
        raise DatasetInputError("The datasets have errors please see report for details.")


if __name__ == "__main__":
    main()


