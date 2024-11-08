
This file outlines the standards and data requirements needed for submitting a dataset to QCArchive.
This ensures that we have a consistent data model for downstream processes.

STANDARDS version: 4 (adopted 2024.11.XX)

We distinguish between standards for the datasets (i.e. the actual data), and the standards for training/fitting and benchmarking/testing force fields.

# Dataset standards

Each new dataset, or version bump, must have a corresponding directory in `qca-dataset-submission/submissions`.
As a rule, we never modify existing datasets, so any changes must be made using a version bump (e.g. removing records with an underlying problem).

## Each molecule must have the following information:
- Canonical isomeric explicit hydrogen mapped SMILES
- Coordinates
- Total charge

## Each dataset's metadata on QCA must have the following: 
- Name (which matches the submission directory on `qca-dataset-submission`
- Version 
- A short description
- A long description, including the following:
    - Intended usage of the dataset
    - Level of theory
    - Elements covered
    - Charges covered
    - Min, mean, and max molecular weight of the molecules in the dataset
    - Zenodo DOI (for public release datasets only, see below on FF releases)
    - Description should be kept up-to-date if e.g. new compute specs are added, using the [`set_description`](https://molssi.github.io/QCFractal/user_guide/records/base.html#qcportal.dataset_models.BaseDataset.set_description) function in QCPortal.
- A link/URL/reference pointing to the provenance to reproduce dataset (e.g. the GH submission folder)
- A github submitter username
- Each entry has canonical isomeric explicit hydrogen mapped SMILES

## Each dataset subdirectory in QCA Dataset submission must have the following information:
- Full environment used to generate dataset
- Dataset generation script(s) including the following:
    - How conformers were generated
    - Origin of SMILES strings
    - Enumerated stereo flag (True/False)
    - Enumerated tautomers flag (True/False)
    - Enumeration keywords/settings
    - Code used to set up the dataset

 - A descriptive README file, as described below

## Each dataset README must contain the following information:

- All information specified in the dataset description/metadata
- A changelog (if there are changes, see the revision section below)
- The name of the person who selected/sourced the molecules 
- A description of the meaning of the entry/molecule keys/names, if not generated automatically by QCSubmit
- A set of elements that the dataset contains
- A set of charges that the dataset contains
- The min, mean and max molecular weight of molecules the dataset contains 

## Each revision must use the following procedure:

A revision means creating a new dataset based on an existing one, with the intent of fixing/improving it.

- The changelog must be copied from the current version (if it exists), and a new entry added for the new version. The changelog should describe any changes between the new version and previous version, including QCArchive record ids of any records that were removed or modified, with an explanation as to why they were removed.
- The dataset version information is updated, and the description updated if necessary
- The README is updated
- A notebook or script used to manipulate the dataset and generate the revision
    - Each file of record should be named with the version that it first addresses, e.g. `submit-v3.0`
        - If `submit-v3.0` also includes a `v3.1 and v3.2` update, then a new file for `v3.2.1` would be `submit-v3.2.1.py`
    - Python notebooks should have version changes in order, and can be run incrementally.
- Update the `index` of datasets on the GH repository

## QCArchive-specific dataset standards

- Compute tag is `openff`
- The computations used for the FITTING standards use the specification named `default`

## Dataset naming and versioning

Each dataset shall be versioned according to `vX.Y.Z`. X refers to the major version, Y to the minor version, and Z to a micro version.
- The naming of a dataset should have the following structure:

    `"OpenFF <descriptive and uniquely-identifying name> v<version number>"`

- The first submission of a dataset will have a version `"v3.0"`

* The major version shall indicate the STANDARDS that the dataset conforms to. Datasets which are not intended to conform to any STANDARDS should start with 0, e.g. `"v0.1"`. 

* Datasets with versions starting with `"v1.x"` and `"v2.x"` do not follow any official STANDARDS, and thus should be considered `"v0.x"`. 

- A minor version change (e.g. `"v3.1"`) represents a minor addition and/or fixes problems:
	- Adding molecules or removing/replacing records with an underlying problem
	- Adding compute specifications
	- Errors/bugs in the molecule specification
	- Changes necessary to adhere to the STANDARDS (i.e. changes necessary to placate the NONE compliance status)

- A micro version change (e.g. `"v3.1.1"`) represents a cosmetic change, or a change that is based on dynamic information that does not change the underlying data:
    - Cosmetic changes
	- Updating the dataset status

This allows the ability to record the version update in the changelog, but not the actual dataset name as this would require a new dataset in QCArchive. 

A best-effort is made to ensure that a dataset follows its underlying STANDARDS. One must assume that the newest version of a dataset best conforms to these STANDARDS, and the same promise may not hold for earlier versions. The changelog should address any changes made to improve compliance.

# Fitting standards

- Reference level of theory: `B3LYP-D3BJ/DZVP`
- Geometry optimization: `geomeTRIC` using the TRIC coordinate system
- QM program: `Psi4`

For unconstrained geometries, all molecules must have:

- Wiberg Bond Orders (parameter interpolation)
- Hessian (vibrational/force constant fitting)

Pre-submission filtering:

- Unless explicitly specified in the submission descriptions, torsion drives must be on four connected atoms
- Torsions driving a ring will give a warning, and torsions in a a ring of  3, 4, 5, or 6 atoms is considered an error
- Warnings will be given if an atom does not have a complete valence set

Post-submission filtering:

- OpenFF toolkit ingestion with strict stereo checking in RDKit
- Hydrogen bonding (for torsion drives only)
- Torsion drives on rings or other high barrier issues
- CMILES (topology) change

## Force field releases

Upon fitting for a new force field release, for the purpose of paper publication, public reference, etc, all molecules should be placed in a single dataset (per type). This gives a single reference for these data instead of many references. All filtering must be done prior, such that all molecules in the release datasets pass all post-submission filters and correspond exactly to the data used for the force field fit.

The format of these dataset names must use the following format:

    `"OpenFF SMIRNOFF <friendly name> <ff version>"`

for example, all datasets (optimizations, torsion drives, and Hessians) with the name `"OpenFF SMIRNOFF Sage 2.0.0"` would refer to all data used to train `openff-2.0.0.offxml` in the `openforcefields` package.

Besides the regular information from the other datasets, these fitting datasets must have:

- `DOI`

### Force field benchmarking

All force field benchmarking will be done using the [YAMMBS](https://github.com/openforcefield/yammbs/tree/main/yammbs) package. Standards will be introduced as our process for these benchmarks is finalized.

