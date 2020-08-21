# OpenFF QCArchive Dataset Submission

All datasets submitted to [QCArchive](https://qcarchive.molssi.org/) via this repository conform to the [Dataset Standards](./STANDARDS.md).
The procedures outlined below are designed to make this as straightforward as possible for submitters.

## Making a submission

**Note**: Making submissions requires being a maintainer on this repository.
          Open an issue if you would like to contribute to this project and make submissions to QCArchive.

To make a new submission, make a PR from a branch on this repository:

1. Clone this repository to your local machine.

2. Create and switch to a new branch:

        git checkout -b <new-submission-branch-name>

3. Create a directory in `./submissions`:

        mkdir submissions/<yyyy-mm-dd>-<new-submission-name>

4. Use QCSubmit to prepare your `dataset.json` in that directory:



5. When validation checks pass and a reviewer has approved your submission, it will be merged.
   Submission to the public QCArchive will follow automatically.


## The lifecycle of a dataset

Following submission, datasets are carried to completion by a well-defined [Dataset Lifecycle](./LIFECYCLE.md).
Much of this lifecycle is automated, with status indicated on the [Dataset Tracking](https://github.com/openforcefield/qca-dataset-submission/projects/1?fullscreen=true) board.

No further action from the submitter is required, unless the submission is placed in "Requires Scientific Review".
In this case, the submitter will be contacted to decide on next steps.


## Adding molecules to a submission



## Adding compute to a submission


