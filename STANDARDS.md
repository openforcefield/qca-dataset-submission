# Dataset Standards

This file outlines the standards and requirements needed for submitting a dataset to QCArchive.
This ensures that we have a consistent data model for downstream processes.

## Required fields 

Current list:
* Ensure all submissions have cmiles, most important are mapped hydrogen smiles
* Ensure the WBO is requested for all submissions, this should be included in the scf properties list using the flag `wiberg_lowdin_indices`

## Best practices
* If any calculations are to be redone from another collection, re-use the old input (coordinates, atom ordering etc) as this will avoid running the calculation again and will just create new references in the database to the old results and should help keep the cost of the calculations down.  

## Dataset naming and versioning

Each dataset shall be versioned.
- The naming of a dataset should have the following structure:

    `"OpenFF <descriptive and uniquely-identifying name> v<version number>"`

- The first submission of a dataset will have a version `"v1.0"`

- A dataset with the suffix `"-beta"` is not to be used for production work.

- A minor version change (e.g. `"v1.1"`) means cosmetic or minor additions/problems were addressed
    - mispelling
    - addition of a e.g. Wiberg bond orders

## Tags indicate status

A tag `"complete"` indicates that a dataset is completed as far as OpenFF is concerned.
This means that any errors remaining are known to be acceptable or impossible to fix.
It also means that no additional work is being done on the dataset to get it to completion.

A tag `"inflight"` indicates that a dataset is not completed as far as OpenFF is concerned.
This means that any errors remaining are being actively addressed.

All datasets should also feature a `"openff"` tag.

## Force Field Releases

When a new force field is released, a dataset corresponding to all results used for the force field fitting should be created.
This gives a single reference for these data instead of many references.
The format of these dataset names is:

    `"OpenFF Force Field <friendly name> <ff version>"`

## Group

The dataset's group should be set to `"OpenFF"`.

## Molecule validation

* See ["Molecule submission checklist"](https://github.com/openforcefield/qcsubmit/issues/9)

## Standard functions and modules for entry preparation

* QCSubmit (https://github.com/openforcefield/qcsubmit)

## Related/ongoing discussions

### Required fields

* See ["Fields that should be required for OpenFF submissions"](https://github.com/openforcefield/qcsubmit/issues/3)
