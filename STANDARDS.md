This file outlines the standards and requirements needed for submitting a dataset to QCArchive. This ensures that we have a consistent data model for downstream processes.

# Required fields 

Current list:
* Ensure all submissions have cmiles, most important are mapped hydrogen smiles
* Ensure the WBO is requested for all submissions, this should be included in the scf properties list using the flag wiberg_lowdin_indices

# Best practices
* If any calculations are to be redone from another collection re-use the old input (coordinates, atom ordering etc) used as this will avoid running the calculation again and will just create new references in the database to the old results and should help keep the cost of the calculations down.  

# Dataset naming and versioning

Each dataset shall be versioned.

# Molecule validation

* See "Molecule submission checklist" https://github.com/openforcefield/qcsubmit/issues/9

# Standard functions and modules for entry preparation

* QCSubmit (https://github.com/openforcefield/qcsubmit)

# Related/ongoing discussions

## Required fields

* See "Fields that should be required for OpenFF submissions" https://github.com/openforcefield/qcsubmit/issues/3
