"""
Script to run validation using qcsubmit on dataset_fragments.json files
"""

import copy
import json
import os
import glob
from argparse import ArgumentParser
import importlib.util
from pprint import pprint

from github import Github
import pandas as pd
import qcportal as ptl
from qcportal.torsiondrive.record_models import TorsiondriveSpecification
from qcportal.optimization.record_models import OptimizationSpecification
from qcportal.singlepoint.record_models import QCSpecification

from openff.qcsubmit.datasets import (BasicDataset, OptimizationDataset,
                                      TorsiondriveDataset,
                                      update_specification_and_metadata)
from openff.qcsubmit.common_structures import SCFProperties, Metadata
from openff.qcsubmit.exceptions import (DatasetInputError, DihedralConnectionError,
                                 LinearTorsionError,
                                 MolecularComplexError, QCSpecificationError, ConstraintError, PCMSettingError)
from openff.qcsubmit.serializers import deserialize

scaffold_validation_path = os.path.join(os.path.dirname(__file__), "scaffold_validation.py")
spec = importlib.util.spec_from_file_location("scaffold_validation", scaffold_validation_path)
val_scfld = importlib.util.module_from_spec(spec)
spec.loader.exec_module(val_scfld)

datasets = { # For dataset*.json
    "dataset": BasicDataset,
    "optimizationdataset": OptimizationDataset,
    "torsiondrivedataset": TorsiondriveDataset
}
ds_specifications = { # For scaffold*.json
    "singlepoint": QCSpecification,
    "optimization": OptimizationSpecification,
    "torsiondrive": TorsiondriveSpecification,
}

check_mark = ":fire:"
cross = ":x:"
uncertain = ":grey_question:"
missing = ":heavy_exclamation_mark:"

QCFRACTAL_URL = "https://api.qcarchive.molssi.org:443/"

REPO_NAME = 'openforcefield/qca-dataset-submission'
DATASET_GLOB = "dataset*.json*"
COMPUTE_GLOB = "compute*.json*"
SCAFFOLD_GLOB = "scaffold*.json*"


def get_data(file_name):
    """Return the deserialized dataset file.
    
    Note that if ``scf_properties`` is at the top level of the output dictionary
    the associated dictionary value is added to each of the ``qc_specification``
    items.
    
    Parameters
    ----------
    file_name : str
        File name and path leading to a .json, or .json.bz2 file.
        
    Returns
    -------
    dict
        Dictionary of dataset contents.
    """
    data = deserialize(file_name=file_name)
    # fix for scf properties moving from dataset to spec
    if "scf_properties" in data and "qc_specifications" in data:
        scf_properties = data.pop("scf_properties")

        # the spec also changed from dict to list at some point
        if isinstance(data["qc_specifications"], dict):
            for spec in data["qc_specifications"].values():
                spec["scf_properties"] = scf_properties
        else:
            for spec in data["qc_specifications"]:
                spec["scf_properties"] = scf_properties
    return data


def create_dataset(dataset_data):
    """Create a QCSubmit dataset from a dictionary

    Parameters
    ----------
    dataset_data : dict
        Dictionary used to create QCSubmit dataset

    Raises
    ------
    RuntimeError
        If dataset type is not supported

    Returns
    -------
    openff.qcsubmit.datasets.*Dataset
        QCSubmit dataset
    """
    if "type" in dataset_data:
        dataset_type = dataset_data["type"]
    elif "dataset_type" in dataset_data:
        dataset_type = dataset_data["dataset_type"]

    dataset_class = datasets.get(dataset_type.lower(), None)
    if dataset_class is not None:
        return dataset_class.parse_obj(dataset_data)
    else:
        raise RuntimeError(f"The dataset type {dataset_type} is not supported.")


def create_spec_report(spec, validated, extras):
    """Create specification report

    Parameters
    ----------
    spec : dict
        Specification dictionary
    validated : bool
        Whether the specification has been validated
    extras : dict
        Additional information to print

    Returns
    -------
    dict
        Report to print
    """

    solvent = spec.get("implicit_solvent", None)

    data = {
        "**Specification Name**": spec["spec_name"],
        "**Method**": spec["method"],
        "**Basis**": spec["basis"],
        "**Wavefunction Protocol**": spec["store_wavefunction"],
        "**Implicit Solvent**": json.dumps(solvent) if solvent is not None else solvent,
        "**Keywords**": json.dumps({} if not spec.get("keywords", None) else spec["keywords"]),
        "**Validated**": validated
    }
    data.update(extras)
    return data

def validate_dataset(dataset_data, flag_scaffold=False):
    """Create a dataset from the data and run normal validation on each molecule.
    
    Catch each of the error types and report them. Convert them into the display output.
    
    Parameters
    ----------
    dataset_data : dict
        Input dataset dictionary
    flag_scaffold : bool, default=False
        Flag whether a qcportal.external.scaffold structure is being used, otherwise a QCSubmit dataset is at hand.

    Returns
    -------
    dict
        Output report dictionary
    """
    errors = {
        "cmiles": [],
        "dihedrals": [],
        "linear": [],
        "complex": [],
        "bonds": [],
        "constraints": [],
        "coordinates": [],
        "total_charge": [],
        "constraint": [],
    }
    data_copy = copy.deepcopy(dataset_data)
    # remove the entries so they can be checked one by one
    entries = data_copy.pop("entries") if flag_scaffold else data_copy.pop("dataset")
    
    if flag_scaffold:
        dataset_type = dataset_data["metadata"]["dataset_type"]
        for entry in entries.values():
            val_result = val_scfld.check_entry(entry, dataset_type)
            for key, check in val_result.items():
                if check:
                    errors[key].append(entry["name"])

    else: # QCSubmit dataset
        # remove the scf props and meta data as this will be checked in a different step
        # try as this is now a per spec property
        data_copy.pop("qc_specifications", None)
        del data_copy["metadata"]
        dataset = create_dataset(data_copy)
    
        # now check each entry
        for entry in entries.values():
            try:
                dataset.add_molecule(**entry, molecule=None)
            except DatasetInputError:
                # this means the cmiles is not valid
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
            except ConstraintError:
                errors["constraint"].append(entry["index"])
        
        errors["coordinates"] = []
        errors["total_charge"] = []
        errors["bonds"] = []

    # print out the errors and the index of any entries which fall into this type
    print("Errors and entries:")
    for error, entries in errors.items():
        print(f"Error type: {error}, entries: {entries}")
    report = {
        "**Valid Cmiles**": cross if errors["cmiles"] else check_mark,
        "**Connected Dihedrals**": cross if errors["dihedrals"] else check_mark,
        "**No Linear Torsions**": cross if errors["linear"] else check_mark,
        "**No Molecular Complexes**": cross if errors["complex"] else check_mark if (not errors["bonds"] or None) else uncertain,
        "**Valid Constraints**": cross if errors["constraints"] else check_mark if (not errors["bonds"] or None) else uncertain,
        "**Total Charge**": cross if errors["total_charge"] else check_mark,
        "**Valid Coordinates**": cross if errors["coordinates"] else check_mark
    }
    return report


def check_metadata(dataset_data, flag_scaffold=False):
    """Ensure dataset metadata meets standards

    Parameters
    ----------
    dataset_data : dict
        Input dataset dictionary
    flag_scaffold : bool, default=False
        Flag whether a qcportal.external.scaffold structure is being used, otherwise a QCSubmit dataset is at hand.

    Returns
    -------
    dict
        Output report dictionary

    Raises
    ------
    DatasetInputError
        Missing metadata fields
    """
    # build and validate the metadata object
    if flag_scaffold:
        try:
            fields = ["description", "tagline", "provenance", "tags",]
            extras_fields = ["submitter", "collection_type", "creation_date", "long_description", "short_description",
                      "long_description_url", "elements", "dataset_name"]
            for field in fields:
                if field not in dataset_data["metadata"] or len(dataset_data["metadata"][field]) == 0:
                    raise ValueError(f"Dataset, {dataset_data['metadata']['name']}, metadata is missing: {field}")
            extras = dataset_data["metadata"]["extras"] if "extras" in dataset_data["metadata"] else dataset_data["metadata"]["metadata"]
            for field in extras_fields:
                if field not in extras or len(extras[field]) == 0:
                    raise ValueError(f"Dataset, {dataset_data['metadata']['name']}, metadata.extras is missing: {field}")
            report = check_mark
        except Exception as e:
            print(str(e))
            report = cross
    else:
        metadata = Metadata.parse_obj(dataset_data["metadata"])
        try:
            metadata.validate_metadata(raise_errors=True)
    
            # QCSubmit no longer validates that `long_description_url` is `None` so we
            # need to check it here.
            if metadata.long_description_url is None:
                raise DatasetInputError(
                    "The metadata has the following incomplete fields "
                    "['long_description_url']"
                )
    
            report = check_mark
        except DatasetInputError:
            report = cross

    return {"**Complete Metatdata**": report}


def check_scf_props(spec):
    # check each scf_prop in the spec
    for scf_property in spec["scf_properties"]:
        try:
            SCFProperties(scf_property)
        except QCSpecificationError:
            return {"**Valid SCF Properties**": cross}
    return {"**Valid SCF Properties**": check_mark}


def check_qcspec_coverage(dataset_data, flag_scaffold=False):
    """For each qcspec try and load it into the dataset, catch an error if not valid.
    
    Also load all elements into the dataset and check coverage.
    
    Parameters
    ----------
    dataset_data : dict
        Input dataset dictionary
    flag_scaffold : bool, default=False
        Flag whether a qcportal.external.scaffold structure is being used, otherwise a QCSubmit dataset is at hand.

    Returns
    -------
    dict
        Output report dictionary
        
    """
    data_copy = copy.deepcopy(dataset_data)
    if flag_scaffold:
        del data_copy["entries"]
        elements = data_copy["metadata"]["extras"]["elements"] if "extras" in data_copy["metadata"] else data_copy["metadata"]["metadata"]["elements"]
        qc_specs = data_copy.pop("specifications")
        spec_report = {}
        for name, spec_dict in qc_specs.items():
            spec = spec_dict["specification"]
            if "optimization_specification" in spec_dict["specification"]: # TD
                keywords = spec_dict["specification"]["optimization_specification"]["qc_specification"]["keywords"]
                program = spec_dict["specification"]["optimization_specification"]["qc_specification"]["program"]
            elif "qc_specification" in spec_dict["specification"]: # Opt
                keywords = spec_dict["specification"]["qc_specification"]["keywords"]
                program = spec_dict["specification"]["qc_specification"]["program"]
            else: # SP
                keywords = spec_dict["specification"]["keywords"]
                program = spec_dict["specification"]["program"]

            valid_scf_props = val_scfld.check_scf_keywords(keywords)
            valid_scf_props = {"**Valid SCF Properties**": check_mark if valid_scf_props else cross}

            try:
                dataset_type = data_copy["metadata"]["dataset_type"]
                if ("ddx" in keywords or "pcm" in keywords) and program.lower() != "psi4":
                    raise ValueError("Implicit solvent is only supported in Psi4")
                _ = ds_specifications[dataset_type](**copy.deepcopy(spec))
                validated = check_mark
            except Exception as e:
                print(str(e))
                validated = cross

            spec_report[name] = val_scfld.create_spec_report(spec_dict, validated, valid_scf_props)
            
        # now get the basis coverage
        all_coverage = val_scfld.check_basis_coverage(qc_specs, elements)
    else:
        # remove any data that could cause an error
        del data_copy["dataset"]
        metadata = data_copy.pop("metadata")
        qc_specs = data_copy.pop("qc_specifications")
        # make the empty dataset and add the elements back
        dataset = create_dataset(data_copy)
        dataset.metadata.elements = metadata["elements"]
        dataset.clear_qcspecs()
        # now try and add each spec
        spec_report = {}
        for spec in qc_specs.values():
            valid_scf_props = check_scf_props(spec)
            try:
                # remove the scf props so they dont cause issues
                # we already validated them explicitly above in the call to `check_scf_props`
                # not necessary to validate them twice, and we want to distinguish them from other errors
                del spec["scf_properties"]
                dataset.add_qc_spec(**spec)
                validated = check_mark
            except (QCSpecificationError, PCMSettingError):
                validated = cross
    
            spec_report[spec["spec_name"]] = create_spec_report(spec, validated, valid_scf_props)
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


def get_meta_info(dataset_data, flag_scaffold=False):
    """Add metadata information to the report dictionary

    Parameters
    ----------
    dataset_data : dict
        Input dataset dictionary
    flag_scaffold : bool, default=False
        Flag whether a qcportal.external.scaffold structure is being used, otherwise a QCSubmit dataset is at hand.

    Returns
    -------
    dict
        Output report dictionary
    """
    if flag_scaffold:
        elements = dataset_data.get("metadata", {}).get("extras", {}).get("elements", missing)
    else:
        elements = dataset_data.get("metadata", {}).get("elements", missing)
    elm_str = " ,".join(elements) if elements != missing else elements

    if flag_scaffold:
        dataset_type = dataset_data["metadata"]["dataset_type"]
        dataset_name = dataset_data["metadata"].get("name", missing)
    else:
        dataset_type = dataset_data["type"] if "type" in dataset_data else dataset_data["dataset_type"]
        dataset_name = dataset_data.get("dataset_name", missing)

    return {"**Dataset Name**": dataset_name,
            "**Dataset Type**": dataset_type,
            "**Elements**": elm_str,
            }


def check_compute_request(dataset_data):
    """Check the compute request.
    
    This will access the archive and check the element coverage and any specs already ran.
    
    Parameters
    ----------
    dataset_data : dict
        Dictionary derived from json file
        
    Returns
    -------
    updated_dataset : dict
        Updated version of input
    spec_report : dict
        Dictionary representing report to output
    """

    qc_specs = dataset_data.pop("qc_specifications")
    dataset = create_dataset(dataset_data)
    client = ptl.PortalClient(QCFRACTAL_URL)

    # now update the dataset with client elements and specs
    updated_dataset = update_specification_and_metadata(dataset=dataset,
                                                        client=client)

    # now we need to try and add each spec; this will raise errors if the spec
    # has already been stored
    spec_report = {}
    for spec in qc_specs.values():
        valid_scf_props = check_scf_props(spec)
        try:
            updated_dataset.add_qc_spec(**spec)
            validated = check_mark
        except QCSpecificationError:
            validated = cross

        spec_report[spec["spec_name"]] = create_spec_report(spec, validated, valid_scf_props)

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
    """Generate a report dataframe for each dataset found.
    
    Parameters
    ----------
    dataset_names : list[str]
        List of dataset names and paths for "dataset*.json*", "compute*.json*", 
        and "scaffold*.json*" to validate.
    
    Returns
    -------
    str
        Output comment compiles for all datasets in the PR.
    
    """
    dataset_dataframe = {}
    qcspec_dataframe = {}

    for dataset_name in dataset_names:
        flag_scaffold = glob.fnmatch.fnmatch(os.path.basename(dataset_name), SCAFFOLD_GLOB)
        dataset_validators = {}
        qc_coverage = None
        # get the data from the dataset
        data = get_data(dataset_name)
        # check if there is a dataset else this might be a compute request
        if ( # will not enter if scaffold
            "dataset" in data and not data["dataset"] and not flag_scaffold or 
            glob.fnmatch.fnmatch(os.path.basename(dataset_name), COMPUTE_GLOB)
        ):
            data, qc_coverage = check_compute_request(data)

        # get the metadata
        dataset_validators.update(get_meta_info(data, flag_scaffold=flag_scaffold))
        # check the first set of entry errors
        dataset_validators.update(validate_dataset(data, flag_scaffold=flag_scaffold))
        # now check the metadata
        dataset_validators.update(check_metadata(data, flag_scaffold=flag_scaffold))
        # now check the qcspec if not already done
        if qc_coverage is None:
            qc_coverage = check_qcspec_coverage(data, flag_scaffold=flag_scaffold)
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
    packages = ["openff.qcsubmit", "openff.toolkit", "basis_set_exchange", "qcelemental"]
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
            if glob.fnmatch.fnmatch(os.path.basename(file), DATASET_GLOB):
                dataset_paths.append(file)
            elif glob.fnmatch.fnmatch(os.path.basename(file), COMPUTE_GLOB):
                dataset_paths.append(file)
            elif glob.fnmatch.fnmatch(os.path.basename(file), SCAFFOLD_GLOB):
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


