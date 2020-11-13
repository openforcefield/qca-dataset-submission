
This file outlines the standards and data requirements needed for submitting a dataset to QCArchive.
This ensures that we have a consistent data model for downstream processes.

STANDARDS version: 3-DRAFT (2020-11-13)

We distinguish between standards for the datasets (i.e. the actual data), and the standards for training force fields.

# Dataset Standards

## Each molecule must have the following information:
- Canonical isomeric explict hydrogen mapped SMILES
- Provenence info of SMILES generation (NEW)
- Coordinates
- Total Charge
- Coordinates must be in CMILES order (discuss)

## Each dataset must have the following information:
- Name
- Version (NEW)
- A short description
- A long description
- A link/URL/reference pointing to the provenence to reproduce dataset (e.g. the GH submission folder)
- A changelog (NEW)
- A github submitter username
- The name of the person who selected/sourced the molecules (NEW)
- A description of the meaning of the entry/molecule keys/names (NEW)
- Each entry has canonical isomeric explict hydrogen mapped SMILES
- Provenence info of CMILES generation (NEW)
- A set of elements that the dataset contains
- A set of charges that the dataset contains (NEW)
- The mean and max molecular weight of molecules the dataset contains (NEW)
- Enumerated stereo flag (True/False) (NEW)
- Enumerated tautomers flag (True/False) (NEW)
- Enumeration provenence info (NEW)
- Computation blacklist (known failures) (NEW)
* Dataset status (NEW):
	- COMPLETE/DONE; all molecules were successful
	- INCOMPLETE/DONE; some molecules were not calculated successfully, and won't be retried
	- INCOMPLETE/WORKING; in progress
	- NONCOMPLIANT/DONE; the dataset does not conform to the standards, and can't be fixed
	- NONCOMPLIANT/WORKING; the dataset does not conform to the standards, is working anyways
	- NONCOMPLIANT/PAUSED; the dataset does not conform to the standards, and calculations have been suspended

## Each dataset README must contain the following information

All information specified in the dataset

## Each revision must use the following procedure:

A revision means creating a new dataset based on an existing one, with the intent of fixing/improving it.

- The changelog must be copied from the current version, and a new entry added in the new version
- The dataset version information is updated
- The README is updated
- A notebook or record of the lines of python used to manipulate the dataset responsible for the revision
	- If the revision should have the new version in the filename e.g. `update-v3.1.1.py` (this does not work for notebooks; discuss)
- Update the `index` of datasets on the GH repo

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
	- Errors/bugs in the molecule specification
	- Changes necessary to adhere to the STANDARDS (i.e. changes necessary to placate the NONCOMPLIANT status)

- A patch (e.g. `"v3.1.1"`) represents a cosmetic change, or a change that is based on dynamic information that does not change the underlying data:
    - Cosmetic changes
	- Updating the blacklist
	- Updating the dataset status

(It might be useful to make the version bump in the changelog, but not the actual dataset name; this would depend on the ability to changing a dataset metadata of an existing dataset. This would allow existing scripts to not need to update to use a new dataset name when e.g. we update the status to COMPLETE/DONE or NONCOMPLIANT/PAUSED->INCOMPLETE/WORKING).

A best-effort is made to ensure that a dataset follows its underlying STANDARDS. One must assume that the newest version of a dataset best conforms to these STANDARDS, and the same promise may not hold for earlier versions. The changelog should address any changes made to improve compliance.

Each version increment should take the information from the previous `changelog` field, and add a new entry of the form { "version": description } that explains the modifications made to the dataset. Each dataset should therefore have the complete changelog.

# Fitting standards

Reference level of theory: `B3LYP/DZVP`
Geometry optimizations: `geomeTRIC` using the TRIC coordinate system
QM program: `Psi4`

For unconstrained geometries, all molecules must have:
- Wiberg Bond Orders (parameter interpolation)
- Hessian (frequency fitting) (is this too restrictive? discuss)

Pre-submission filtering:
- Unless explicitly specified in the submission descriptions, torsion drives must be on 4 connected atoms
- Torsions driving a ring will give a warning, and torsions in a a ring of  3,4,5,6 atoms is considered an error
- Warnings will be given if an atom does not have a complete valence set

Post-submission filtering: (we need testers/validators for these)
	- OpenFF toolkit ingestion with strict stereo checking in RDKit
	- Hydrogen bonding
	- Torsion drives on rings or other high barrier issues
	- CMILES change (no way to build new CMILES since connectivity cannot be trusted, also this means we can't trust per-molecule CMILES until this filtering check passes)

## Force Field Releases

Upon fitting for a new force field release, for the purpose of paper publication, public reference, etc, all molecules should be placed in a single dataset (per type). This gives a single reference for these data instead of many references. Filtering must be done prior, such that all molecules in the release datasets pass all post-submission filters.

The format of these dataset names must use the following format:

    `"OpenFF SMIRNOFF <friendly name> <ff version>"`

for example, all datasets (optimizations, torsion drives, and Hessians) with the name `"OpenFF SMIRNOFF Sage 2.0.0"` would refer to all data used to train `openff-2.0.0.offxml` in the `openforcefields` package.

Besides the regular information from the other datasets, these fitting datasets must have:

- `DOI`
- Dataset standards version (this document)

