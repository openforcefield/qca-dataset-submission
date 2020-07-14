# The Lifecycle of a Dataset Submission


Writing a Github Actions workflow that will run a script on a schedule.
The script will:
- given a list of datasets, get status of each one from public QCArchive
- write that status out in the mapped anchor issue/PR for that dataset


On creation of PR:
- validation tests run

On merge of PR:
- submission to public QCA
- creation of tracking issue


Make it flag based.
If a repo has some kind of metadata flag, indicating that it's being tracked, track it.
Tracking issue should follow a convention.


## Dataset states

<strike>
1. States are managed with the Treant category 'state'.
    - 'submitted'
    - 'error cycling'
    - 'requires scientific review'
    - 'end of life'
    - 'archival/complete'

2. States map directly to states in qca-dataset-submission "Dataset Tracking" project.
</strike>

3. States move forward.
   They do not move backward.
    - we want this to be a soft-assumption of our model
    - can be violated, but generally shouldn't be

4. Automation will:
    1. Gather tagged 'tracking' Treants.
    2. Instantiate `DataSet` object.
    3. Check state of corresponding issue on project.
    4. Apply operations appropriate for that state.
        - these should be idempotent, since state 
          can be changed on project board.
        - security implications?
            - should be fine, because board not editable by non-members

5. State advancment will be a manual process for now.
    - more of it will be automated with time.

## Proposed operational pathway

`lifecycle` will run daily
- could separate out into multiple

1. A PR is created against `qca-dataset-submission` by a submitter.
    - the template is filled out with informational sections (TBD)
    - [GHA] `validation` operates on `spec.json`; performs validation checks
        - comment made based on validation checks
        - reruns on each push

2. Trigger for addition to `Dataset Tracking` project.
    - add 'tracking' tag to PR
    - [GHA] `lifecycle` will add card to "Backlog" state for PR if not yet there.

3. When dataset is ready for submission to public QCA, merge PR.
    - [board automation] PR card will move to state "Queued for Submission" 
    - [GHA] `lifecycle` will move PR card to state "Queued for Submission" if merged and in state "Backlog"
    - [GHA] `lifecycle` will attempt to submit the dataset
        - if successful, will move card to state "Error Cycling"; add comment to PR
        - if failed, will keep card queued; add comment to PR; attempt again next day

4. COMPLETE, INCOMPLETE, ERROR numbers reported for `Optimizations`, `TorsionDrives`
    - [GHA] `lifecycle` will collect the above statistics for state "Error Cycling" PRs
    - [GHA] `lifecycle` will restart all errored `Optimizations` and `TorsionDrives` for state "Error Cycling" PRs
    - [GHA] `lifecycle` will move PR to state "Archived/Complete" if no ERROR, INCOMPLETE, all COMPLETE

5. PR will remain in state "Error Cycling" until moved to "Requires Scientific Review" or until all tasks COMPLETE
    - if errors appear persistent, human manager should move to state "Requires Scientific Review"
    - discussion should be had on PR for next version
    - once decided, state moved to "End of Life"

6. [GHA] `lifecycle` should mark datasets in states "End of Life" and "Archived/Complete" in some meaningful way on public QCArchive

### Issues with the approach

1. Do we loop through PRs with the `tracking` tag?
   How do we make the connection to the corresponding directory?
    - only need to make this connection for 
        - submission
        - collection name, dataset name.

    Should we specify that this needs to be in a comment?
        - for a `tracking`-tagged PR without a comment, 

    Could we parse the PR data for directory impacted?
        - what happens if multiple PRs use the same directory?
            - at worst, error cycling would happen in series for the same dataset
            - generally avoidable if we put one into scientific review, end of life, etc.
                - or just untag `tracking`

2. Multiple dataset submissions / directories made in a single PR
    - discouraged, but allowed
        - scientific review and error cycling performed *on all*, at same time
        - essentially ties them together

3. For legacy datasets, add flagpost file with name.
    - use this to indicate it is a dataset
    - make new PR for each, adding flagpost
    - move immediately to error cycling
    - test on one of the recent torsiondrive series

4. For new datasets, flagpost is `spec.json`
    - name directly extracted from there
