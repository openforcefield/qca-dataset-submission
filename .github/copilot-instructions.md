# GitHub Copilot Code Review Instructions

When performing an automated review of pull requests in this repository, the review tool (such as GitHub Copilot Code Review) should prioritize compliance with the dataset standards as defined in our STANDARDS.md (a file found in the top level of this repository, which a reviewer should read before performing a review) and associated policies.  
The goal is to provide feedback that helps maintain a consistent, reproducible, and high-quality dataset submission workflow.

## Key Areas to Check

### 1. Dataset Submission Structure

Confirm that any new dataset (or version bump) is placed in its own directory under qca-dataset-submission/submissions/YYY-MM-DD-<dataset name with spaces replaced with hyphens>.  

Ensure existing dataset*.json files or scaffold*.json files are not being directly modified; instead changes should go via a version bump in the same directory.  

Check for existence of required sub-directories and files (dataset generation scripts, README.md, conda environment yaml).

### 2. Molecule and Metadata Standards in *.py and *.ipynb* File

Metadata is defined in the *.ipynb or *.py file used to create the dataset. Do not review dataset*.json.bz2 or scaffold*.json.bz2 LFS files.

For each molecule entry: verify that it has a canonical isomeric explicit hydrogen-mapped SMILES ("CMILES") string, coordinate data, and total charge.  

For each dataset metadata entry: verify that it includes name, version, short description, long description with usage, method, basis, elements & charges covered, min/mean/max molecular weight, DOI (if applicable), provenance link, GitHub submitter username.  

Ensure metadata names use spaces for dataset names, and versioning uses the vX.Y.Z scheme as specified.

The dataset url should be “https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/<YYYY-MM-DD>-<dataset name with spaces replaced with hyphens>“

The dataset url should be included in the metadata as 'long_description_url'

The repository README.md should be updated to have a new row at the bottom of one of the tables. If a singlepoint dataset type, update the table in Basic Datasets section; if an optimization dataset type, update the table in the Optimization Datasets section; if a torsiondrive dataset type, update the table in the TorsionDrive Datasets section. A new row should include:
| <dataset name> | [<YYYY-MM-DD>-<dataset name with spaces replaced with hyphens>](<dataset url>)                                 | <dataset short description (e.g. tagline)> | <elements in the dataset taken from the metadata> | |

### 3. README and Changelog

Check that README.md in the PR dataset directory contains: dataset name, tagline (e.g. short description), long description (all matching metadata), changelog (if version > initial), name of person who sourced molecules, description of molecule keys/names, list of elements, list of charges, min/mean/max molecular weight. The QCEngine process and subprocess with the keywords used in the *.ipynb or *.py file used to generate the dataset.

The title at the top of the README.md should be identical to the dataset name.

If the submission directory already exists, the PR is a revision. For revisions: ensure changelog is updated with new entry describing the changes (including record IDs removed/modified and explanations).

Ensure all files in the directory are listed and defined in the README.md.

### 4. Versioning & Naming Convention

Confirm versioning is consistent: major version X corresponds to the STANDARDS version (e.g., v4.x for STANDARDS version 4) defined in STANDARDS.md in the top level of this repo.  

Check dataset naming structure: dataset names should start with “OpenFF” (or “OpenFF SMIRNOFF” for force-field releases) if the version is >=v4.0, the dataset not should also use spaces.

For datasets not conforming to any version of STANDARDS, version should start with v0.x.  This is common in cases where the molecules do not have valid CMILES (add in other common reasons for 0.X versions)

Minor version bumps (v4.1, v4.2, etc) should reflect minor additions/fixes rather than major redesigns.


## Tone and Format of Comments

Begin the review with a warning that your review is experimental and should be superseded by the contents of the STANDARDS.md doc.

Provide a short summary at the top of the review comment (e.g., “The new dataset directory <insert dataset name> mostly follows the standards, but ...”)  

Then list specific issues or improvements as bullet points, referencing file paths and dataset names when possible.  

If the change is small and fully compliant, you may respond: “Looks good — all mandatory fields present and versioning correct. 🔥”

If certain checks cannot be automatically verified (e.g., “Were conformers generated as per enumeration keywords?”) note them as manual verification needed.

Provide all needed information but avoid being verbose.

## Scope / Focus

Focus primarily on structural and metadata compliance with the STANDARDS.md document. The information may be formatted as a nested list or structured text; exact formatting is not prescribed.

Do not (in this automated review) attempt deep scientific correctness (e.g., verifying actual theory computation results) — leave that to domain experts.  

Ensure that submission workflows are reproducible and consistent (scripts (e.g. *.ipynb or *.py), environment yaml, README.md present).  

Flag missing or mis-named files, inconsistent versioning, absent required metadata, naming issues, missing changelog entries for versions beyond the initial submission.

## Feedback & Follow-Up

For each flagged issue, include a note on how to fix it (e.g., “Add the DOI in the long description field of the metadata and update README accordingly”).  

If the pull request includes a version bump, verify the change directory exists, version increments correctly, metadata version updated, README updated, and index of datasets in repository updated.  

If no issues in adherence to the standards can be found, indicate positive verification.

## Comment Resolution and Verification

When reviewing a PR or subsequent comments:
- If a maintainer has responded to your comment with clarification or justification, mark the comment as resolved if their response adequately addresses the concern
- If a file has been updated to address your comment, verify the change and mark as resolved
- Do NOT leave comments unresolved when:
  - The maintainer explains why the current approach is compliant with STANDARDS.md
  - The STANDARDS.md explicitly allows for the approach taken
  - Your initial comment was based on a misinterpretation of the standards
- Focus on objective compliance with STANDARDS.md, not on style preferences or suggestions beyond the requirements