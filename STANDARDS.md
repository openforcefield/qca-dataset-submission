
This file outlines the standards and data requirements needed for submitting a dataset to QCArchive.
This ensures that we have a consistent data model for downstream processes.

STANDARDS version: 3 (adopted 2020.12.11)

We distinguish between standards for the datasets (i.e. the actual data), and the standards for training/fitting and benchmarking/testing force fields.

# Dataset standards

## Each molecule must have the following information:
- Canonical isomeric explicit hydrogen mapped SMILES
- Provenance information of SMILES generation (NEW)
- Coordinates
- Total charge

## Each dataset must have the following information:
- Name
- Version (NEW)
- A short description
- A long description
- A link/URL/reference pointing to the provenance to reproduce dataset (e.g. the GH submission folder)
- A changelog (NEW)
	- The changelog is a python-like dictionary of the form `{ version: entry }`
	- Each entry in the changelog has the python dictionary form `{ person: str, date: str, description: str }`
- A github submitter username
- The name of the person who selected/sourced the molecules (NEW)
- A description of the meaning of the entry/molecule keys/names (NEW)
- Each entry has canonical isomeric explicit hydrogen mapped SMILES
- Provenance info of CMILES generation
- A set of elements that the dataset contains
- A set of charges that the dataset contains (NEW)
- The mean and max molecular weight of molecules the dataset contains (NEW)
- Enumerated stereo flag (True/False) (NEW)
- Enumerated tautomers flag (True/False) (NEW)
- Enumeration provenance info (NEW)
- Computation blacklist (known failures) (NEW)
- Dataset status flags (NEW):
	- Complete: `COMPLETE` `INCOMPLETE`
	- State: `DONE` `WORKING` `PAUSED`
	- Compliance: `V3` `NONE`
	- COMPLETE/DONE/V3; all molecules were successful
	- INCOMPLETE/DONE/V3; some molecules were not calculated successfully, and won't be retried
	- INCOMPLETE/WORKING/V3; in progress
	- INCOMPLETE/DONE/NONE; the dataset does not conform to the standards, and can't be fixed
	- INCOMPLETE/WORKING/NONE; the dataset does not conform to the standards, but is working anyways
	- INCOMPLETE/PAUSED/NONE; the dataset does not conform to the standards, and calculations have been suspended

## Each dataset README must contain the following information

All information specified in the dataset

## Each revision must use the following procedure:

A revision means creating a new dataset based on an existing one, with the intent of fixing/improving it.

- The changelog must be copied from the current version, and a new entry added in the new version
- The dataset version information is updated
- The README is updated
- A notebook or record of the lines of python used to manipulate the dataset responsible for the revision
	- Each file of record should be named with the version that it first addresses, e.g. `submit-v3.0`
		- If `submit-v3.0` also includes a `v3.1 and v3.2` update, then a new file for `v3.2.1` would be `submit-v3.2.1.py`
	- Python notebooks should have version changes in order, and can be run incrementally.
- Update the `index` of datasets on the GH repository

## QCArchive-specific dataset standards

- Compute tag is `openff`
- The computations used for the FITTING standards use the specification named `default`

## Dataset naming and versioning

Each dataset shall be versioned.
- The naming of a dataset should have the following structure:

    `"OpenFF <descriptive and uniquely-identifying name> v<version number>"`

- The first submission of a dataset will have a version `"v3.0"`

* The major version shall indicate the STANDARDS that the dataset conforms to. Datasets which are not intended to conform to any STANDARDS should start with 0, e.g. `"v0.1"`. 

* Datasets with versions starting with `"v1.x"` and `"v2.x"` do not follow any official STANDARDS, and thus should be considered `"v0.x"`.

- A minor version change (e.g. `"v3.1"`) represents a minor addition and/or fixes problems:
	- Adding molecules
	- Adding compute specifications
	- Errors/bugs in the molecule specification
	- Changes necessary to adhere to the STANDARDS (i.e. changes necessary to placate the NONE compliance status)

- A micro version change (e.g. `"v3.1.1"`) represents a cosmetic change, or a change that is based on dynamic information that does not change the underlying data:
    - Cosmetic changes
	- Updating the blacklist
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
- Hydrogen bonding
- Torsion drives on rings or other high barrier issues
- CMILES (topology) change

## Force field releases

Upon fitting for a new force field release, for the purpose of paper publication, public reference, etc, all molecules should be placed in a single dataset (per type). This gives a single reference for these data instead of many references. Filtering must be done prior, such that all molecules in the release datasets pass all post-submission filters.

The format of these dataset names must use the following format:

    `"OpenFF SMIRNOFF <friendly name> <ff version>"`

for example, all datasets (optimizations, torsion drives, and Hessians) with the name `"OpenFF SMIRNOFF Sage 2.0.0"` would refer to all data used to train `openff-2.0.0.offxml` in the `openforcefields` package.

Besides the regular information from the other datasets, these fitting datasets must have:

- `DOI`

### Force field benchmarking

These use an MM compute specification. The name of the specification must be the name of the FF file, e.g. `openff-2.0.0`

In these cases, the following pre-submission checks must be successful:

- OpenFF toolkit ingestion with strict stereo checking in RDKit
- The OpenFF toolkit can successfully create an OpenMM system

