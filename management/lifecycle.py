#!/usr/bin/env python

import os
import json
import traceback
from collections import defaultdict
from datetime import datetime

from github import Github

REPO_NAME = "openforcefield/qca-dataset-submission"
DATASET_FILENAME = "dataset.json"


class Submission:
    """A submission, corresponding to a single PR, possibly multiple datasets.

    """

    def __init__(self, pr, ghapi, repo=None):
        """Create a new Submission instance that performs operations on PR
        card state based on data in the PR itself.

        Since a submission can have multiple DataSets tied together, this
        allows for control of card state based on what's going on in the
        collection of DataSets the PR is linked to.

        Parameters
        ----------
        pr : github.PullRequest
            PullRequest corresponding to the dataset submission.
        ghapi : github.Github
            An authenticated Github Python API client object.
        repo : str
            Github repo where datasets are tracked.

        """
        self.pr = pr
        self.ghapi = ghapi

        if repo is None:
            self.repo = ghapi.get_repo(REPO_NAME)
        else:
            self.repo = repo

        self.datasets = self._gather_datasets()

    def _gather_datasets(self):
        files = self.pr.get_files()
        datasets = [
            file.filename for file in files if DATASET_FILENAME in file.filename
        ]
        return datasets

    @staticmethod
    def _get_board_card_state(board, pr):
        pr_state = None
        pr_card = None
        for state, cards in board.items():
            for card in cards:
                if card.get_content().number == pr.number:
                    pr_state = state
                    pr_card = card
                    break

        return pr_card, pr_state

    @staticmethod
    def _get_column(repo, column):
        proj = [
            proj for proj in repo.get_projects() if proj.name == "Dataset Tracking"
        ][0]

        cols = list(proj.get_columns())
        return [col for col in cols if col.name == column][0]

    def set_backlog(self):
        backlog = self._get_column(self.repo, "Backlog")
        backlog.create_card(content_id=self.pr.id, content_type="PullRequest")

    def execute_state(self, board=None, states=None):
        """Based on current state of the PR, perform appropriate actions.

        """
        if board is None:
            board = _get_full_board(self.repo)

        pr_card, pr_state = self._get_board_card_state(board, self.pr)

        # if card not on board, then it starts in the Backlog
        if pr_state is None:
            pr_state = self.set_backlog()

            # reload board, since we just added this card
            board = _get_full_board(self.repo)
            pr_card, pr_state = self._get_board_card_state(board, self.pr)

        # exit early if states specified, and this PR is not
        # in one of those
        if states is not None:
            if pr_state not in states:
                return

        if pr_state == "Backlog":
            return self.execute_backlog(pr_card, pr_state)
        elif pr_state == "Queued for Submission":
            return self.execute_queued_submit(pr_card, pr_state)
        elif pr_state == "Error Cycling":
            return self.execute_errorcycle(pr_card, pr_state)
        elif pr_state == "Requires Scientific Review":
            return self.execute_requires_scientific_review(pr_card, pr_state)
        elif pr_state == "End of Life":
            return self.execute_end_of_life(pr_card, pr_state)
        elif pr_state == "Archived/Complete":
            return self.execute_archived_complete(pr_card, pr_state)

    def resolve_new_state(self, dataset_results):
        """If new state agreed upon by dataset results, that state is returned.
        Otherwise, returns `None`.
        """
        # get unique states recommended by datasets for this PR
        # may not always be the same, say, if e.g. submission fails for one
        # of many datasets in this submission
        new_card_state = set(res["new_state"] for res in dataset_results)

        # if all datasets agree on the new card state, we change to that state
        if len(new_card_state) == 1:
            new_state = list(new_card_state)[0]
            return new_state
        else:
            return None

    def evolve_state(self, pr_card, pr_state, new_state):
        # no need to move if we are already in the new state
        if pr_state != new_state:
            state_col = self._get_column(self.repo, new_state)
            pr_card.move(position="top", column=state_col)

    def execute_backlog(self, pr_card, pr_state):
        """If PR is in the backlog and is merged, it will get moved to the
        queued for submission state.

        """
        if self.pr.is_merged():

            comment = f"""
            ## Lifecycle - Backlog

            Merged dataset moved from "Backlog" to "Queued for Submission".

            """

            # postprocess due to raw spacing above
            comment = "\n".join([substr.strip() for substr in comment.split("\n")])

            # submit comment
            self.pr.create_issue_comment(comment)

            self.evolve_state(pr_card, pr_state, "Queued for Submission")

            return {"new_state": "Queued for Submission"}
        else:
            return {"new state": "Backlog"}

    def execute_queued_submit(self, pr_card, pr_state):
        """Submit datasets, perhaps with some retry logic.

        """
        results = []
        for dataset in self.datasets:
            print(f"Processing dataset '{dataset}'")
            ds = DataSet(dataset, self, self.ghapi)
            results.append(ds.execute_queued_submit())

        new_state = self.resolve_new_state(results)
        if new_state is not None:
            self.evolve_state(pr_card, pr_state, new_state)

    def execute_errorcycle(self, pr_card, pr_state):
        """Error cycle each dataset

        """
        results = []
        for dataset in self.datasets:
            print(f"Processing dataset '{dataset}'")
            ds = DataSet(dataset, self, self.ghapi)
            results.append(ds.execute_errorcycle())

        new_state = self.resolve_new_state(results)
        if new_state is not None:
            self.evolve_state(pr_card, pr_state, new_state)

        if new_state == "Archived/Complete":
            for dataset in self.datasets:
                dataset.comment_archived_complete()


class DataSet:
    """A dataset submitted to QCArchive.
    
    A dataset has a lifecycle with well-defined states.
    This class represents the current state of a dataset,
    and provides the machinery for execution of lifecycle processes based on that state.

    All lifecycle state is stored on Github in the original PR for the submission,
    mapped onto states in the "Datset Tracking" project board.
    
    """

    def __init__(self, dataset, submission, ghapi, repo=None):
        """Create new DataSet instance linking a submission dataset to its PR.

        Parameters
        ----------
        dataset : path-like
            Path to dataset submission file.
        submission : Submission
            Submission instance corresponding to the dataset submission.
        ghapi : github.Github
            An authenticated Github Python API client object.
        repo : str
            Github repo where datasets are tracked.

        """
        self.dataset = dataset
        self.submission = submission
        self.pr = submission.pr
        self.ghapi = ghapi

        if repo is None:
            self.repo = ghapi.get_repo(REPO_NAME)
        else:
            self.repo = repo

    def _parse_spec(self):
        with open(self.dataset, "r") as f:
            spec = json.load(f)

        dataset_name = spec["dataset_name"]
        dataset_type = spec["dataset_type"]

        return dataset_name, dataset_type

    def _get_qca_client(self):
        import qcportal as ptl

        client = ptl.FractalClient(
            username=os.environ["QCA_USER"], password=os.environ["QCA_KEY"]
        )

        return client

    def _get_meta(self):
        import pandas as pd

        datehr = datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC")
        dataset_name, dataset_type = self._parse_spec()

        meta = {
            "**Dataset Name**": dataset_name,
            "**Dataset Type**": dataset_type,
            "**UTC Datetime**": datehr,
        }

        return pd.DataFrame(pd.Series(meta, name=""))

    def _version_info_report(self):
        version = get_version_info()

        comment = f"""
        <details>
        <summary><b>QCSubmit</b> version information(<i>click to expand</i>)</summary>
        <!-- have to be followed by an empty line! -->

        {version.to_markdown()}
        </details>
        """

        return comment

    def execute_queued_submit(self, max_retries=3):
        """Submit the dataset, perhaps with some retry logic.

        """
        from qcsubmit.serializers import deserialize

        client = self._get_qca_client()

        # load dataset into QCSubmit class
        ds = deserialize(self.dataset)
        dataset_qcs = create_dataset(ds)

        try:
            # submit to QCArchive
            output = dataset_qcs.submit(client)
            self._queued_submit_report(output, success=True)
        except:
            self._queued_submit_report(traceback.format_exc(), success=False)
            return {"new_state": "Queued for Submission"}
        else:
            return {"new_state": "Error Cycling"}

    def _queued_submit_report(self, output, success):
        success_text = "**SUCCESS**" if success else "**FAILED**"

        comment = f"""
        ## Lifecycle - QCSubmit Submission Report : {success_text}

        {self._get_meta().to_markdown()}
        
        Response from public QCArchive:

        ```
        {output}
        ```

        ----------
        {self._version_info_report()}

        """

        # postprocess due to raw spacing above
        comment = "\n".join([substr.strip() for substr in comment.split("\n")])

        # submit comment
        self.pr.create_issue_comment(comment)

    def execute_errorcycle(self, restart=False):
        """Obtain complete, incomplete, error stats for dataset and report.
        
        For suspected random errors, we perform restarts.

        If dataset complete, move state to "Archived/Complete".

        """
        client = self._get_qca_client()

        dataset_name, dataset_type = self._parse_spec()
        ds = client.get_collection(dataset_type, dataset_name)

        if dataset_type.lower() == "TorsionDriveDataset".lower():
            complete = self._errorcycle_torsiondrive(ds, client)

        elif dataset_type.lower() == "OptimizationDataset".lower():
            complete = self._errorcycle_optimization(ds, client)

        elif dataset_type.lower() == "GridOptimizationDataset".lower():
            complete = self._errorcycle_gridopt(ds, client)

        elif dataset_type.lower() == "Dataset".lower():
            complete = self._errorcycle_dataset(ds, client)

        if complete:
            return {"new_state": "Archived/Complete"}
        else:
            return {"new_state": "Error Cycling"}

    def comment_archived_complete(self):
        comment = f"""
        ## Lifecycle - Archived/Complete

        {self._get_meta().to_markdown()}

        **Dataset Complete!**

        """

        # postprocess due to raw spacing above
        comment = "\n".join([substr.strip() for substr in comment.split("\n")])

        # submit comment
        self.pr.create_issue_comment(comment)

    def _errorcycle_torsiondrive(self, ds, client):
        import management as mgt

        tdrs, df_tdr = self._errorcycle_torsiondrive_get_tdr_errors(ds, client)
        opts, df_tdr_opt = self._errorcycle_torsiondrive_get_tdr_opt_errors(ds, client)

        opt_error_counts = mgt.count_unique_optimization_error_messages(
            opts, full=True, pretty_print=True, tolerate_missing=True
        )

        self._errorcycle_torsiondrive_report(df_tdr, df_tdr_opt, opt_error_counts)

        if (df_tdr[["RUNNING", "INCOMPLETE", "ERROR"]].sum().sum() == 0) and (
            df_tdr_opt[["INCOMPLETE", "ERROR"]].sum().sum() == 0
        ):
            complete = True
        else:
            # restart errored torsiondrives and optimizations
            self._errorcycle_restart_optimizations(opts, client)
            self._errorcycle_restart_torsiondrives(tdrs, client)
            complete = False

        return complete

    def _errorcycle_torsiondrive_get_tdr_errors(self, ds, client):
        import pandas as pd
        import management as mgt

        # gather torsiondrive results
        results = defaultdict(dict)
        for spec in ds.list_specifications().index.tolist():
            tdrs = mgt.get_torsiondrives(ds, spec, client)

            for status in ["COMPLETE", "RUNNING", "INCOMPLETE", "ERROR"]:
                results[spec][status] = len(
                    [tdr for tdr in tdrs if tdr.status == status]
                )

        df = pd.DataFrame(results).transpose()
        df.index.name = "specification"
        return tdrs, df

    def _errorcycle_torsiondrive_get_tdr_opt_errors(self, ds, client):
        import pandas as pd
        import management as mgt

        # gather torsiondrive optimization results
        results = defaultdict(dict)
        for spec in ds.list_specifications().index.tolist():
            opts = mgt.merge(mgt.get_torsiondrive_optimizations(ds, spec, client))

            for status in ["COMPLETE", "INCOMPLETE", "ERROR"]:
                results[spec][status] = len(
                    [opt for opt in opts if opt.status == status]
                )

        df = pd.DataFrame(results).transpose()
        df.index.name = "specification"
        return opts, df

    def _errorcycle_torsiondrive_report(self, df_tdr, df_tdr_opt, opt_error_counts):
        comment = f"""
        ## Lifecycle - Error Cycling Report

        {self._get_meta().to_markdown()}

        All errored tasks and services will be restarted.
        Errored states prior to restart reported below.

        ### `TorsionDriveRecord` current status

        {df_tdr.to_markdown()}

        ### `OptimizationRecord` current status

        {df_tdr_opt.to_markdown()}

        #### `OptimizationRecord` Error Tracebacks:

        ```
        {opt_error_counts}
        ```

        ----------
        {self._version_info_report()}

        """

        # postprocess due to raw spacing above
        comment = "\n".join([substr.strip() for substr in comment.split("\n")])

        # submit comment
        self.pr.create_issue_comment(comment)

    def _errorcycle_restart_optimizations(self, opts, client):
        import management as mgt

        # TODO: add some nuance for the types of optimzations
        # we will *not* restart, such as SCF convergence issues
        mgt.restart_optimizations(opts, client)

    def _errorcycle_restart_torsiondrives(self, tdrs, client):
        import management as mgt

        mgt.restart_torsiondrives(tdrs, client)

    def _errorcycle_optimization(self, ds, client):
        import management as mgt

        opts, df_opt = self._errorcycle_optimization_get_opt_errors(ds, client)

        opt_error_counts = mgt.count_unique_optimization_error_messages(
            opts, full=True, pretty_print=True, tolerate_missing=True
        )

        self._errorcycle_optimization_report(df_opt, opt_error_counts)

        if df_opt[["INCOMPLETE", "ERROR"]].sum().sum() == 0:
            complete = True
        else:
            # restart errored optimizations
            self._errorcycle_restart_optimizations(opts, client)
            complete = False

        return complete

    def _errorcycle_optimization_get_opt_errors(self, ds, client):
        import pandas as pd
        import management as mgt

        # gather optimization results
        results = defaultdict(dict)
        for spec in ds.list_specifications().index.tolist():
            opts = mgt.get_optimizations(ds, spec, client)

            for status in ["COMPLETE", "INCOMPLETE", "ERROR"]:
                results[spec][status] = len(
                    [opt for opt in opts if opt.status == status]
                )

        df = pd.DataFrame(results).transpose()
        df.index.name = "specification"
        return opts, df

    def _errorcycle_optimization_report(self, df_opt, opt_error_counts):
        comment = f"""
        ## Lifecycle - Error Cycling Report

        {self._get_meta().to_markdown()}

        All errored tasks will be restarted.
        Errored states prior to restart reported below.

        ### `OptimizationRecord` current status

        {df_opt.to_markdown()}

        #### `OptimizationRecord` Error Tracebacks:

        ```
        {opt_error_counts}
        ```

        ----------
        {self._version_info_report()}

        """

        # postprocess due to raw spacing above
        comment = "\n".join([substr.strip() for substr in comment.split("\n")])

        # submit comment
        self.pr.create_issue_comment(comment)

    def execute_requires_scientific_review(self):
        pass

    def execute_end_of_life(self):
        pass

    def execute_archived_complete(self):
        pass


def create_dataset(dataset_data):
    from qcsubmit.datasets import BasicDataset, OptimizationDataset, TorsiondriveDataset

    datasets = {
        "BasicDataset": BasicDataset,
        "OptimizationDataset": OptimizationDataset,
        "TorsiondriveDataset": TorsiondriveDataset,
    }

    dataset_type = dataset_data["dataset_type"]
    dataset_class = datasets.get(dataset_type, None)
    if dataset_class is not None:
        return dataset_class.parse_obj(dataset_data)
    else:
        raise RuntimeError(f"The dataset type {dataset_type} is not supported.")


def _get_full_board(repo):
    proj = [proj for proj in repo.get_projects() if proj.name == "Dataset Tracking"][0]
    board = {col.name: [card for card in col.get_cards()] for col in proj.get_columns()}
    return board


def _get_tracking_prs(repo):
    prs = [
        pr
        for pr in repo.get_pulls(state="all")
        if "tracking" in list(map(lambda x: x.name, pr.labels))
    ]
    return prs


def get_version_info():
    """
    Get the version info for the packages used to validate the submission.
    """
    import importlib
    import pandas as pd

    report = {}
    # list the core packages here
    packages = ["qcsubmit", "openforcefield", "basis_set_exchange", "qcelemental"]
    for package in packages:
        module = importlib.import_module(package)
        report[package] = pd.Series({"version": module.__version__})

    # now try openeye else use rdkit
    try:
        import openeye

        report["openeye"] = pd.Series({"version": openeye.__version__})
    except ImportError:
        import rdkit

        report["rdkit"] = pd.Series({"version": rdkit.__version__})

    return pd.DataFrame(report).transpose()


def main():
    """Map PRs tagged with 'tracking' into corresponding datasets.

    """
    import argparse

    parser = argparse.ArgumentParser(
        description="Process PRs according to dataset lifecycle"
    )
    parser.add_argument(
        "--states",
        type=str,
        nargs="*",
        help="states to limit processing to; if not provided, use all states",
    )
    parser.add_argument(
        "--prs",
        type=int,
        nargs="*",
        help="PR numbers to limit processing to; if not provided, all labeled PRs processed",
    )

    args = parser.parse_args()

    if args.states:
        states = args.states
    else:
        states = None

    if args.prs:
        prnums = args.prs
    else:
        prnums = None

    gh = Github(os.environ["GH_TOKEN"])
    repo = gh.get_repo(REPO_NAME)

    # gather up all PRs with the `tracking` label
    tracking_prs = _get_tracking_prs(repo)

    # filter on PR numbers, if provided
    if prnums is not None:
        prs = []
        for tpr in tracking_prs:
            if tpr.number in prnums:
                prs.append(tpr)
    else:
        prs = tracking_prs

    print(f"Found {len(prs)} with the 'tracking' label")

    # grab the full project board state once so we don't have to hammer the API
    # over and over
    board = _get_full_board(repo)

    # for each PR, we examine the changes to find files used for the submission
    # this is where the mapping is made between the PR and the submission files
    for pr in prs:
        print(f"Processing PR #{pr.number}")

        submission = Submission(pr, gh)
        submission.execute_state(board=board, states=states)


if __name__ == "__main__":
    main()
