# The Lifecycle of a Dataset Submission

All Open Force Field datasets submitted to QCArchive undergo well-defined *lifecycle*.

![Dataset Lifecycle](assets/lifecycle-diagram.png)

Each labeled rectangle in the lifecycle represents a *state*. 
A submission PR changes state according to the arrows.
Changes in state may be performed by automation or manually by a human when certain critera are met.

The lifecycle process is described below, with [bracketed] items indicating the agent of action, one of:
- [GHA]: Github Actions
- [Board]: Github Project Board
- [Human]: A maintainer of the `qca-dataset-submission` repository.

1. A PR is created against `qca-dataset-submission` by a submitter.
    - the template is filled out with informational sections according to the [PR template](.github/pull_request_template.md)
    - [GHA]  `validation` operates on all `dataset.json` files found in the PR; performs validation checks
        - comment made based on validation checks
        - reruns on each push

2. Add card for the PR to [Dataset Tracking](https://github.com/openforcefield/qca-dataset-submission/projects/1) board.
    - [Human]  add 'tracking' tag to PR
    - [GHA]  [`lifecycle-backlog`](.github/workflows/lifecycle-backlog.yml) will add card to ["Backlog"](https://github.com/openforcefield/qca-dataset-submission/projects/1#column-9577334) state for PR if not yet there.

3. When dataset is ready for submission to public QCArchive (validations pass, submitters and reviewers satisfied), PR is merged.
    - [Board]  PR card will move to state ["Queued for Submission"](https://github.com/openforcefield/qca-dataset-submission/projects/1#column-9577335) immediately.
    - [GHA]  [`lifecycle-backlog`](.github/workflows/lifecycle-backlog.yml) will move PR card to state ["Queued for Submission"](https://github.com/openforcefield/qca-dataset-submission/projects/1#column-9577335) if merged and in state ["Backlog"](https://github.com/openforcefield/qca-dataset-submission/projects/1#column-9577334)
    - [GHA]  [`lifecycle-submission`](.github/workflows/lifecycle-submission.yml) will attempt to submit the dataset
        - if successful, will move card to state ["Error Cycling"](https://github.com/openforcefield/qca-dataset-submission/projects/1#column-9577365); add comment to PR
        - if failed, will keep card queued; add comment to PR; attempt again next execution 

4. COMPLETE, INCOMPLETE, ERROR numbers reported for `Optimizations`, `TorsionDrives`
    - [GHA]  [`lifecycle-error-cycle`](.github/workflows/lifecycle-error-cycle.yml) will collect the above statistics for state ["Error Cycling"](https://github.com/openforcefield/qca-dataset-submission/projects/1#column-9577365) PRs
        - will restart all errored `Optimizations` and `TorsionDrives`
        - will move PR to state ["Archived/Complete"](https://github.com/openforcefield/qca-dataset-submission/projects/1#column-9577372) if no ERROR, INCOMPLETE, all COMPLETE

5. PR will remain in state ["Error Cycling"](https://github.com/openforcefield/qca-dataset-submission/projects/1#column-9577365) until moved to ["Requires Scientific Review"](https://github.com/openforcefield/qca-dataset-submission/projects/1#column-9577358) or until all tasks COMPLETE
    - [Human]  if errors appear persistent,  move to state ["Requires Scientific Review"](https://github.com/openforcefield/qca-dataset-submission/projects/1#column-9577358)
    - discussion should be had on PR for next version
    - [Human]  once decided, state moved to ["End of Life"](https://github.com/openforcefield/qca-dataset-submission/projects/1#column-9577336)

6. [GHA]  `lifecycle-end-of-life` will add tag 'end-of-life' to dataset in QCArchive for PR in ["End of Life"](https://github.com/openforcefield/qca-dataset-submission/projects/1#column-9577336)

7. [GHA]  `lifecycle-archived-complete` will add tag 'archived-complete' to dataset in QCArchive for PR in ["Archived/Complete"](https://github.com/openforcefield/qca-dataset-submission/projects/1#column-9577372)


### Issues with the approach

1. Do we loop through PRs with the `tracking` tag?
   How do we make the connection to the corresponding directory?
    - only need to make this connection for 
        - submission
        - collection name, dataset name.

    - what happens if multiple PRs use the same directory?
        - at worst, error cycling would happen in series for the same dataset
        - generally avoidable if we put one into scientific review, end of life, etc.
            - or just untag `tracking`

2. Multiple dataset submissions / directories made in a single PR
    - discouraged, but allowed
        - scientific review and error cycling performed *on all*, at same time
        - essentially ties them together

3. For legacy datasets, add `dataset.json` file with name, type.
    - use this to indicate it is a dataset
    - make new PR for each
    - move immediately to error cycling
