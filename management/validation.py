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
    Convert them into the display output.
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
            errors["cmiles"].append(entry["index"])
            # remove the index error after the qcsubmit patch
        except (DihedralConnectionError, IndexError):
            # the torsion is not connected
            errors["dihedrals"].append(entry["index"])
        except LinearTorsionError:
            errors["linear"].append(entry["index"])
        except MolecularComplexError:
            errors["complex"].append(entry["index"])

    report = {
        "**Valid Cmiles**": cross if errors["cmiles"] else check_mark,
        "**Connected Dihedrals**": cross if errors["dihedrals"] else check_mark,
        "**No Linear Torsions**": cross if errors["linear"] else check_mark,
        "**No Molecular Complexes**": cross if errors["complex"] else check_mark
    }
    return report


def check_metadata(dataset_data):
    # remove the dataset and scf props
    data_copy = copy.deepcopy(dataset_data)
    del data_copy["scf_properties"]
    del data_copy["dataset"]
    dataset = create_dataset(data_copy)
    try:
        dataset.metadata.validate_metadata(raise_errors=True)
        report = check_mark
    except DatasetInputError:
        report = cross

    return {"**Complete Metatdata**": report}


def check_scf_props(dataset_data):
    # remove the metadata and dataset
    data_copy = copy.deepcopy(dataset_data)
    del data_copy["dataset"]
    del data_copy["metadata"]
    try:
        _ = create_dataset(data_copy)
        report = check_mark
    except DatasetInputError:
        report = cross

    return {"**Valid SCF Properties**": report}


def check_basis_coverage(dataset_data):
    "Note this will only work if the full dataset can be loaded."

    try:
        dataset = create_dataset(dataset_data)
    except Exception:
        return {"**Full Basis Coverage**": missing}
    try:
        dataset._get_missing_basis_coverage(raise_errors=True)
        report = check_mark
    except MissingBasisCoverageError:
        report = cross

    return {"**Full Basis Coverage**": report}


def get_meta_info(dataset_data):
    return {"**Dataset Name**": dataset_data.get("dataset_name", missing),
            "**Dataset Type**": dataset_data.get("dataset_type", missing),
            "**Method**": dataset_data.get("method", missing),
            "**Basis**": dataset_data.get("basis", missing)}


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
        dataset_validators = {}
        # get the data from the dataset
        data = get_data(dataset_name)
        # get the metadata
        dataset_validators.update(get_meta_info(data))
        # check the first set of entry errors
        dataset_validators.update(validate_dataset(data))
        # now check the metadata
        dataset_validators.update(check_metadata(data))
        # now check the scf
        dataset_validators.update(check_scf_props(data))
        # now check the basis
        dataset_validators.update(check_basis_coverage(data))

        dataset_dataframe[dataset_name] = pd.Series(dataset_validators)
    # now make the dataframe
    df = pd.DataFrame(data=dataset_dataframe)

    # now construct the comment with the qcsubmit version info
    version = get_version_info()
    comment = f"""
    ## QCSubmit Validation Report
    
    {df.to_markdown()}

    <details>
    <summary><b>QCSubmit</b> version information(<i>click to expand</i>)</summary>
    <!-- have to be followed by an empty line! -->

    {version.to_markdown()}
    </details>   
    """

    # postprocess due to raw spacing above
    comment = "\n".join([substr.strip() for substr in comment.split('\n')])

    return comment


def get_version_info():
    """
    Get the version info for the packages used to validate the submission.
    """
    import importlib
    report = {}
    # list the core packages here
    packages = ["qcsubmit", "openforcefield", "basis_set_exchange", "qcelemental"]
    for package in packages:
        module = importlib.import_module(package)
        report[package] = pd.Series({"version": module.__version__})

    # now try openeye else use rdkit
    try:
        import openeye
        report["openeye"] = pd.Series({"version": openeye.__version__})
    except ImportError:
        import rdkit
        report["rdkit"] = pd.Series({"version": rdkit.__version__})

    return pd.DataFrame(report).transpose()


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


