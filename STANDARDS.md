This file outlines the standards and requirements needed for submitting a dataset to QCArchive.
This ensures that we have a consistent data model for downstream processes.

STANDARDS version: 1

Requires:
	- QCEngine 0.13
	- Other QCA pieces needed to understand the options below

# Required fields 

A human readable description of the fields and why they are necessary

## Required fields layout

- All datasets
	- The name
	- The metadata
	* The `default` specification:
		- `PSI4/B3LYP-D3BJ/DZVP`
			- scf_properties
				- `wiberg_lowdin_indices`
				- `mayer_indices`
	* The entry:
		- `attributes`
			- `canonical_isomeric_hydrogen_mapped_cmiles`
		- The molecules
			- A valid QCSchema molecule with the follwing also required:
				- `extras`
					- `canonical_isomeric_hydrogen_mapped_cmiles`
				- `connectivity`

- Optimizations
	- The qc_spec options

	- The optimizer options
		- `program: geometric`
			- The geometric options

	
- TorsionDrives
	- The qc_spec options

	- The optimizer options
		- `program: geometric`
			- The geometric options

- Hessians
	- The qc_spec options

- GridOptimizations
	- The qc_spec options

* Training sets
	- Necessary contributed values

* Benchmarking sets
	- Necessary contributes values

## Job specifications (level of theory; settings)

OpenFF depends on a QCSpecification named "default" which corresponds to `B3LYP-D3BJ/DZVP` in Psi4. Submissions may have multiple specifications, but must include the `default`.

# Best practices

* If any calculations are (to be redone from another collection, re-use the old input (coordinates, atom ordering etc) as this will avoid running the calculation again and will just create new references in the database to the old results and should help keep the cost of the calculations down.  

# Dataset naming and versioning

Each dataset shall be versioned.
- The naming of a dataset should have the following structure:

    `"OpenFF <descriptive and uniquely-identifying name> v<version number>"`

- The first submission of a dataset will have a version `"v1.0"`

* The major version shall indicate the STANDARDS that the dataset conforms to. Datasets which are not intended to conform to these STANDARDS should start with 0, e.g. `"v0.1"`

- A dataset with the suffix `"-beta"` is not to be used for production work.

- A minor version change (e.g. `"v1.1"`) means cosmetic or minor additions/problems were addressed
    - Cosmetic changes
	- Errors/bugs in the molecule specification
	- Changes necessary to adhere to the STANDARDS

# Tags indicate status

A tag `"complete"` indicates that a dataset is completed as far as OpenFF is concerned.
This means that any errors remaining are known to be acceptable or impossible to fix.
It also means that no additional work is being done on the dataset to get it to completion.

A tag `"inflight"` indicates that a dataset is not completed as far as OpenFF is concerned.
This means that any errors remaining are being actively addressed.

All datasets should also feature a `"openff"` tag.

# Force Field Releases 

When a new force field is released, a dataset corresponding to all results used for the force field fitting should be created.
This gives a single reference for these data instead of many references.
The format of these dataset names is:

    `"OpenFF Force Field <friendly name> <ff version>"`

# Group

The dataset's group should be set to `"OpenFF"`.

# Molecule validation

* See ["Molecule submission checklist"](https://github.com/openforcefield/qcsubmit/issues/9)

* Unique keys in `ds.data.records` must not reference the same entry.
* Unless explicitly specified in the submission descriptions, torsion drives must be on 4 connected atoms
* Torsions driving a ring will give a warning, and torsions in a a ring of  3,4,5,6 atoms is considered an error
* Warnings will be given if an atom does not have a complete valence set

# Standard functions and modules for entry preparation

* QCSubmit (https://github.com/openforcefield/qcsubmit)

# Related/ongoing discussions

## Required fields

* See ["Fields that should be required for OpenFF submissions"](https://github.com/openforcefield/qcsubmit/issues/3)

