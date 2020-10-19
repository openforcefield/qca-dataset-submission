"""
Script to run validation using qcsubmit on dataset_fragments.json files
"""

import copy
import json
import os
import glob
from argparse import ArgumentParser


import pandas as pd
from github import Github
from qcsubmit.datasets import (BasicDataset, OptimizationDataset,
                               TorsiondriveDataset)
from qcsubmit.exceptions import (DatasetInputError, DihedralConnectionError,
                                 LinearTorsionError,
                                 MolecularComplexError, QCSpecificationError, ConstraintError)
from qcsubmit.serializers import deserialize
from qcsubmit.utils import update_specification_and_metadata
import qcportal as ptl

datasets = {
    "dataset": BasicDataset,
    "optimizationdataset": OptimizationDataset,
    "torsiondrivedataset": TorsiondriveDataset}

check_mark = ":heavy_check_mark:"
cross = ":x:"
missing = ":heavy_exclamation_mark:"
REPO_NAME = 'openforcefield/qca-dataset-submission'


def get_data(file_name):
    """
    Return the deserialized dataset file.
    """
    return deserialize(file_name=file_name)


def create_dataset(dataset_data):
    dataset_type = dataset_data["dataset_type"]
    dataset_class = datasets.get(dataset_type.lower(), None)
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
        "complex": [],
        "constraints": [],
    }
    data_copy = copy.deepcopy(dataset_data)
    # remove the entries so they can be checked one by one
    entries = data_copy.pop("dataset")
    # remove the scf props and meta data as this will be checked in a different step
    del data_copy["scf_properties"]
    del data_copy["metadata"]
    del data_copy["qc_specifications"]
    dataset = create_dataset(data_copy)

    # now check each entry
    for entry in entries.values():
        try:
            dataset.add_molecule(**entry, molecule=None)
        except DatasetInputError:
            # this mean the cmiles is not valid
            errors["cmiles"].append(entry["index"])
            # remove the index error after the qcsubmit patch
        except DihedralConnectionError:
            # the torsion is not connected
            errors["dihedrals"].append(entry["index"])
        except LinearTorsionError:
            errors["linear"].append(entry["index"])
        except MolecularComplexError:
            errors["complex"].append(entry["index"])
        except ConstraintError:
            errors["constraints"].append(entry["index"])

    report = {
        "**Valid Cmiles**": cross if errors["cmiles"] else check_mark,
        "**Connected Dihedrals**": cross if errors["dihedrals"] else check_mark,
        "**No Linear Torsions**": cross if errors["linear"] else check_mark,
        "**No Molecular Complexes**": cross if errors["complex"] else check_mark,
        "**Valid Constraints**": cross if errors["constraints"] else check_mark
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


def check_qcspec_coverage(dataset_data):
    """
    For each qcspec try and load it into the dataset, catch an error if not valid.
    Also load all elements into the dataset and check coverage.
    """
    data_copy = copy.deepcopy(dataset_data)
    # remove any data that could cause an error
    del data_copy["dataset"]
    del data_copy["scf_properties"]
    metadata = data_copy.pop("metadata")
    qc_specs = data_copy.pop("qc_specifications")
    # make the empty dataset and add the elements back
    dataset = create_dataset(data_copy)
    dataset.metadata.elements = metadata["elements"]
    dataset.clear_qcspecs()
    # now try and add each spec
    spec_report = {}
    for spec in qc_specs.values():
        try:
            dataset.add_qc_spec(**spec)
            validated = check_mark
        except QCSpecificationError:
            validated = cross

        spec_report[spec["spec_name"]] = {"**Specification Name**": spec["spec_name"], "**Method**": spec["method"],
                                          "**Basis**": spec["basis"], "**Wavefunction Protocol**": spec["store_wavefunction"], "**Validated**": validated}

    # now get the basis coverage
    all_coverage = dataset._get_missing_basis_coverage(raise_errors=False)
    # now we need to update each report
    for key, report in spec_report.items():
        coverage = all_coverage.get(key, missing)
        if coverage == missing:
            spec_report[key]["**Full Basis Coverage**"] = coverage
        elif coverage:
            spec_report[key]["**Full Basis Coverage**"] = cross
        else:
            spec_report[key]["**Full Basis Coverage**"] = check_mark

    return spec_report


def get_meta_info(dataset_data):
    elements = dataset_data.get("metadata", {}).get("elements", missing)
    if elements != missing:
        elm_str = " ,".join(elements)
    else:
        elm_str = elements

    return {"**Dataset Name**": dataset_data.get("dataset_name", missing),
            "**Dataset Type**": dataset_data.get("dataset_type", missing),
            "**Elements**": elm_str,
            }


def check_compute_request(dataset_data):
    """
    Check the compute request this will access the archive and check the element coverage and any specs already ran.
    """
    qc_specs = dataset_data.pop("qc_specifications")
    dataset = create_dataset(dataset_data)
    client = ptl.FractalClient()
    # now update the dataset with client elements and specs
    updated_dataset = update_specification_and_metadata(dataset=dataset, client=client)
    if updated_dataset.metadata.validate_metadata(raise_errors=False):
        # fill in missing fields
        updated_dataset.metadata.long_description_url = "https://test.org"
    # now we need to try and add each spec this will raise errors if the spec has already been stored
    spec_report = {}
    for spec in qc_specs.values():
        try:
            updated_dataset.add_qc_spec(**spec)
            validated = check_mark
        except QCSpecificationError:
            validated = cross

        spec_report[spec["spec_name"]] = {"**Specification Name**": spec["spec_name"], "**Method**": spec["method"],
                                          "**Basis**": spec["basis"],
                                          "**Wavefunction Protocol**": spec["store_wavefunction"],
                                          "**Validated**": validated}

    # now get the basis coverage
    all_coverage = dataset._get_missing_basis_coverage(raise_errors=False)
    # now we need to update each report
    for key, report in spec_report.items():
        coverage = all_coverage.get(key, missing)
        if coverage == missing:
            spec_report[key]["**Full Basis Coverage**"] = coverage
        elif coverage:
            spec_report[key]["**Full Basis Coverage**"] = cross
        else:
            spec_report[key]["**Full Basis Coverage**"] = check_mark

    return updated_dataset.dict(), spec_report


def main_validation(dataset_names):
    """
    Generate a report dataframe for each dataset found.
    """
    dataset_dataframe = {}
    qcspec_dataframe = {}

    for dataset_name in dataset_names:
        dataset_validators = {}
        qc_coverage = None
        # get the data from the dataset
        data = get_data(dataset_name)
        # check if there is a dataset else this might be a compute request
        if not data["dataset"] or "compute.json" in dataset_name:
            data, qc_coverage = check_compute_request(data)

        # get the metadata
        dataset_validators.update(get_meta_info(data))
        # check the first set of entry errors
        dataset_validators.update(validate_dataset(data))
        # now check the metadata
        dataset_validators.update(check_metadata(data))
        # now check the scf
        dataset_validators.update(check_scf_props(data))
        # now check the qcspec if not already done
        if qc_coverage is None:
            qc_coverage = check_qcspec_coverage(data)
        for key, coverage in qc_coverage.items():
            name = f"{dataset_name}/{key}"
            qcspec_dataframe[name] = pd.Series(coverage)

        dataset_dataframe[dataset_name] = pd.Series(dataset_validators)

    # now make the dataframe
    metatadata = pd.DataFrame(data=dataset_dataframe)
    qcspec = pd.DataFrame(data=qcspec_dataframe)

    # now construct the comment with the qcsubmit version info
    version = get_version_info()
    comment = f"""
    ## QCSubmit Validation Report
    
    {metatadata.to_markdown()}

    ## QC Specification Report

    {qcspec.to_markdown()}

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
        # this covers files that are deleted and picked up by the file change check
        if os.path.exists(file):
            if "dataset.json" in file:
                dataset_paths.append(file)
            elif glob.fnmatch.fnmatch(os.path.basename(file), "compute*.json"):
                dataset_paths.append(file)
        else:
            continue

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


