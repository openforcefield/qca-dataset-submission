
#!/usr/bin/env python

import os
import glob
import json
import traceback
from collections import defaultdict
from datetime import datetime

from github import Github

REPO_NAME = "openforcefield/qca-dataset-submission"
DATASET_GLOB = "dataset*.json*"
COMPUTE_GLOB = "compute*.json*"

PRIORITIES = {'priority-low': 0, 'priority-normal': 1, 'priority-high': 2}


class Submission:
    """A submission, corresponding to a single PR, possibly multiple datasets.

    A submission has a lifecycle with well-defined states.
    This class represents the current state of a submission,
    and provides the machinery for execution of lifecycle processes based on that state.

    All lifecycle state is stored on Github in the original PR for the submission,
    mapped onto states in the "Datset Tracking" project board.

    """

    def __init__(self, pr, ghapi, repo=None, priority=1, computetag='openff'):
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
        priority : int
            Priority to use for the dataset if set by method calls;
            one of 0, 1, or 2, in increasing-priority order.
        computetag : str
            Compute tag to use for the dataset if set by method calls;
            tasks with a given compute tag will only be computed by managers
            configured to service that tag.

        """
        self.pr = pr
        self.ghapi = ghapi
        self.priority = priority
        self.computetag = computetag

        if repo is None:
            self.repo = ghapi.get_repo(REPO_NAME)
        else:
            self.repo = repo

        self.datasets = self._gather_datasets()
        self.computes = self._gather_computes()

    def _gather_datasets(self):
        files = self.pr.get_files()
        datasets = list(filter(
            lambda x: glob.fnmatch.fnmatch(os.path.basename(x), DATASET_GLOB),
            map(lambda x: x.filename, files)))

        # we only want files that actually exist
        # it can rarely be the case that a PR features changes to a path that is a file deletion
        datasets = [ds for ds in datasets if os.path.exists(ds)]

        return datasets

    def _gather_computes(self):
        files = self.pr.get_files()
        computes = list(filter(
            lambda x: glob.fnmatch.fnmatch(os.path.basename(x), COMPUTE_GLOB),
            map(lambda x: x.filename, files)))

        # we only want files that actually exist
        # it can rarely be the case that a PR features changes to a path that is a file deletion
        computes = [cs for cs in computes if os.path.exists(cs)]

        return computes

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

    def execute_state(self, board=None, states=None,
                      reset_errors=False, set_priority=False,
                      set_computetag=False):
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
            return self.execute_errorcycle(pr_card, pr_state,
                    reset_errors=reset_errors, set_priority=set_priority,
                    set_computetag=set_computetag)
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

        for compute in self.computes:
            print(f"Processing compute '{compute}'")
            ct = Compute(compute, self, self.ghapi)
            results.append(ct.execute_queued_submit())

        new_state = self.resolve_new_state(results)
        if new_state is not None:
            self.evolve_state(pr_card, pr_state, new_state)

    def execute_errorcycle(self, pr_card, pr_state,
                           reset_errors=False,
                           set_priority=False,
                           set_computetag=False):
        """Error cycle each dataset

        """
        results = []
        for dataset in self.datasets:
            print(f"Processing dataset '{dataset}'")
            ds = DataSet(dataset, self, self.ghapi,
                         priority=self.priority, computetag=self.computetag)
            results.append(ds.execute_errorcycle(reset_errors=reset_errors,
                                                 set_priority=set_priority,
                                                 set_computetag=set_computetag))

        for compute in self.computes:
            print(f"Processing compute '{compute}'")
            ct = Compute(compute, self, self.ghapi,
                         priority=self.priority, computetag=self.computetag)
            results.append(ct.execute_errorcycle(reset_errors=reset_errors,
                                                 set_priority=set_priority,
                                                 set_computetag=set_computetag))

        new_state = self.resolve_new_state(results)
        if new_state is not None:
            self.evolve_state(pr_card, pr_state, new_state)

        if new_state == "Archived/Complete":
            for dataset in self.datasets:
                ds = DataSet(dataset, self, self.ghapi)
                ds.comment_archived_complete()


class SubmittableBase:
    def __init__(self, submittable, submission, ghapi, repo=None,
                 priority=1, computetag='openff'):
        """Create new Submittable instance linking a submission dataset to its PR.

        Parameters
        ----------
        submittable : path-like
            Path to submission file.
        submission : Submission
            Submission instance corresponding to the dataset submission.
        ghapi : github.Github
            An authenticated Github Python API client object.
        repo : str
            Github repo where datasets are tracked.
        priority : int
            Priority to use for the dataset if set by method calls;
            one of 0, 1, or 2, in increasing-priority order.
        computetag : str
            Compute tag to use for the dataset if set by method calls;
            tasks with a given compute tag will only be computed by managers
            configured to service that tag.

        """
        self.submittable = submittable
        self.submission = submission
        self.pr = submission.pr
        self.ghapi = ghapi
        self.priority = priority
        self.computetag = computetag

        if repo is None:
            self.repo = ghapi.get_repo(REPO_NAME)
        else:
            self.repo = repo

    def _parse_spec(self):
        spec = self._load_submittable()

        dataset_name = spec["dataset_name"]

        if "type" in spec:
            dataset_type = spec["type"]
        elif "dataset_type" in spec:
            dataset_type = spec["dataset_type"]

        dataset_specs = spec.get("qc_specifications", None)

        return dataset_name, dataset_type, dataset_specs

    def _load_submittable(self):
        from openff.qcsubmit.serializers import deserialize
        spec = deserialize(self.submittable)

        return spec

    def _get_qca_client(self):
        import qcportal as ptl

        client = ptl.FractalClient(
            username=os.environ["QCA_USER"], password=os.environ["QCA_KEY"]
        )

        return client

    def _get_meta(self):
        import pandas as pd

        datehr = datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC")
        dataset_name, dataset_type, dataset_specs = self._parse_spec()

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
        """Submit, perhaps with some retry logic.

        """
        from openff.qcsubmit.serializers import deserialize

        client = self._get_qca_client()

        # load dataset into QCSubmit class
        ds = self._load_submittable()
        dataset_qcs = create_dataset(ds)

        try:
            # submit to QCArchive
            output = self.submit(dataset_qcs, client)
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

    def execute_errorcycle(self,
                           reset_errors=False,
                           set_priority=False,
                           set_computetag=False):
        """Obtain complete, incomplete, error stats for submittable and report.

        For suspected random errors, we perform restarts.

        If submittable complete, recommend state "Archived/Complete".

        """
        client = self._get_qca_client()

        dataset_name, dataset_type, dataset_specs = self._parse_spec()
        ds = client.get_collection(dataset_type, dataset_name)

        if dataset_type.lower() == "TorsionDriveDataset".lower():
            complete = self._errorcycle_torsiondrive(ds, client, dataset_specs,
                    reset_errors=reset_errors, set_priority=set_priority,
                    set_computetag=set_computetag)

        elif dataset_type.lower() == "OptimizationDataset".lower():
            complete = self._errorcycle_optimization(ds, client, dataset_specs,
                    reset_errors=reset_errors, set_priority=set_priority,
                    set_computetag=set_computetag)

        elif dataset_type.lower() == "GridOptimizationDataset".lower():
            complete = self._errorcycle_gridopt(ds, client, dataset_specs,
                    reset_errors=reset_errors, set_priority=set_priority,
                    set_computetag=set_computetag)

        elif dataset_type.lower() == "Dataset".lower():
            complete = self._errorcycle_dataset(ds, client, dataset_specs,
                    reset_errors=reset_errors, set_priority=set_priority,
                    set_computetag=set_computetag)

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

    def _errorcycle_torsiondrive(self, ds, client, dataset_specs,
            reset_errors=False, set_priority=False, set_computetag=False):
        import management as mgt

        tdrs, df_tdr = self._errorcycle_torsiondrive_get_tdr_errors(ds, client, dataset_specs)
        opts, df_tdr_opt = self._errorcycle_torsiondrive_get_tdr_opt_errors(ds, client, dataset_specs)

        if reset_errors:
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
            if reset_errors:
                self._errorcycle_restart_optimizations(opts, client)
                self._errorcycle_restart_torsiondrives(tdrs, client)
            if set_priority:
                self._set_priority_optimizations(opts, client)
                self._set_priority_torsiondrives(tdrs, client)
            if set_computetag:
                self._set_computetag_optimizations(opts, client)
                self._set_computetag_torsiondrives(tdrs, client)
            complete = False

        return complete

    def _errorcycle_torsiondrive_get_tdr_errors(self, ds, client, dataset_specs):
        import pandas as pd
        import management as mgt

        if dataset_specs is None:
            dataset_specs = ds.list_specifications().index.tolist()

        # gather torsiondrive results
        results = defaultdict(dict)
        all_tds = list()
        for spec in dataset_specs:
            tdrs = mgt.get_torsiondrives(ds, spec, client)
            all_tds.extend(tdrs)

            for status in ["COMPLETE", "RUNNING", "INCOMPLETE", "ERROR"]:
                results[spec][status] = len(
                    [tdr for tdr in tdrs if tdr.status == status]
                )

        df = pd.DataFrame(results).transpose()
        df.index.name = "specification"
        return all_tds, df

    def _errorcycle_torsiondrive_get_tdr_opt_errors(self, ds, client, dataset_specs):
        import pandas as pd
        import management as mgt

        if dataset_specs is None:
            dataset_specs = ds.list_specifications().index.tolist()

        # gather torsiondrive optimization results
        results = defaultdict(dict)
        all_opts = list()
        for spec in dataset_specs:
            opts = mgt.merge(mgt.get_torsiondrive_optimizations(ds, spec, client))
            all_opts.extend(opts)

            for status in ["COMPLETE", "INCOMPLETE", "ERROR"]:
                results[spec][status] = len(
                    [opt for opt in opts if opt.status == status]
                )

        df = pd.DataFrame(results).transpose()
        df.index.name = "specification"
        return all_opts, df

    def _errorcycle_torsiondrive_report(self, df_tdr, df_tdr_opt, opt_error_counts):

        if len(opt_error_counts) > 60000:
            opt_error_counts = opt_error_counts[:60000]
            opt_error_counts += "\n--- Too many errors; truncated here ---\n"

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

        <details>
        <summary><b>Tracebacks</b> (<i>click to expand</i>)</summary>
        <!-- have to be followed by an empty line! -->

        ```
        {opt_error_counts}
        ```
        </details>

        ----------
        {self._version_info_report()}

        """

        # postprocess due to raw spacing above
        comment = "\n".join([substr.strip() for substr in comment.split("\n")])

        # submit comment
        self.pr.create_issue_comment(comment)

    def _errorcycle_restart_results(self, res, client):
        import management as mgt

        mgt.restart_results(res, client)

        # handle cases where an Result has status INCOMPLETE, but data is there
        mgt.regenerate_results(res, client)

    def _errorcycle_restart_optimizations(self, opts, client):
        import management as mgt

        # TODO: add some nuance for the types of optimzations
        # we will *not* restart, such as SCF convergence issues
        mgt.restart_optimizations(opts, client)

        # handle cases where an Optimization has status INCOMPLETE, but data is there
        mgt.regenerate_optimizations(opts, client)

    def _errorcycle_restart_torsiondrives(self, tdrs, client):
        import management as mgt

        mgt.restart_torsiondrives(tdrs, client)

    def _set_priority_results(self, results, client):
        import management as mgt
        mgt.reprioritize_results(results, client, self.priority)

    def _set_priority_optimizations(self, opts, client):
        import management as mgt
        mgt.reprioritize_optimizations(opts, client, self.priority)

    def _set_priority_torsiondrives(self, tdrs, client):
        # TODO: no way to reprioritize services at the moment
        pass

    def _set_computetag_results(self, results, client):
        import management as mgt
        mgt.retag_results(results, client, self.computetag)

    def _set_computetag_optimizations(self, opts, client):
        import management as mgt
        mgt.retag_optimizations(opts, client, self.computetag)

    def _set_computetag_torsiondrives(self, tdrs, client):
        # TODO: no way to retag services at the moment
        pass

    def _errorcycle_optimization(self, ds, client, dataset_specs,
            reset_errors=False, set_priority=False, set_computetag=False):
        import management as mgt

        opts, df_opt = self._errorcycle_optimization_get_opt_errors(ds, client, dataset_specs)

        if reset_errors:
            opt_error_counts = mgt.count_unique_optimization_error_messages(
                opts, full=True, pretty_print=True, tolerate_missing=True
            )

            self._errorcycle_optimization_report(df_opt, opt_error_counts)

        if df_opt[["INCOMPLETE", "ERROR"]].sum().sum() == 0:
            complete = True
        else:
            if reset_errors:
                # restart errored optimizations
                self._errorcycle_restart_optimizations(opts, client)
            if set_priority:
                self._set_priority_optimizations(opts, client)
            if set_computetag:
                self._set_computetag_optimizations(opts, client)
            complete = False

        return complete

    def _errorcycle_optimization_get_opt_errors(self, ds, client, dataset_specs):
        import pandas as pd
        import management as mgt

        if dataset_specs is None:
            dataset_specs = ds.list_specifications().index.tolist()

        # gather optimization results
        results = defaultdict(dict)
        all_opts = list()
        for spec in dataset_specs:
            opts = mgt.get_optimizations(ds, spec, client)
            all_opts.extend(opts)

            for status in ["COMPLETE", "INCOMPLETE", "ERROR"]:
                results[spec][status] = len(
                    [opt for opt in opts if opt.status == status]
                )

        df = pd.DataFrame(results).transpose()
        df.index.name = "specification"
        return all_opts, df

    def _errorcycle_optimization_report(self, df_opt, opt_error_counts):

        if len(opt_error_counts) > 60000:
            opt_error_counts = opt_error_counts[:60000]
            opt_error_counts += "\n--- Too many errors; truncated here ---\n"

        comment = f"""
        ## Lifecycle - Error Cycling Report

        {self._get_meta().to_markdown()}

        All errored tasks will be restarted.
        Errored states prior to restart reported below.

        ### `OptimizationRecord` current status

        {df_opt.to_markdown()}

        #### `OptimizationRecord` Error Tracebacks:

        <details>
        <summary><b>Tracebacks</b> (<i>click to expand</i>)</summary>
        <!-- have to be followed by an empty line! -->

        ```
        {opt_error_counts}
        ```
        </details>

        ----------
        {self._version_info_report()}

        """

        # postprocess due to raw spacing above
        comment = "\n".join([substr.strip() for substr in comment.split("\n")])

        # submit comment
        self.pr.create_issue_comment(comment)

    def _errorcycle_dataset(self, ds, client, dataset_specs,
            reset_errors=False, set_priority=False, set_computetag=False):
        import management as mgt
        res, df_res = self._errorcycle_dataset_get_result_errors(ds, client, dataset_specs)

        if reset_errors:
            res_error_counts = mgt.count_unique_result_error_messages(
                res, full=True, pretty_print=True, tolerate_missing=True
            )

            self._errorcycle_dataset_report(df_res, res_error_counts)

        if df_res[["INCOMPLETE", "ERROR"]].sum().sum() == 0:
            complete = True
        else:
            if reset_errors:
                # restart errored tasks
                self._errorcycle_restart_results(res, client)
            if set_priority:
                self._set_priority_results(res, client)
            if set_computetag:
                self._set_computetag_results(res, client)

            complete = False

        return complete

    def _errorcycle_dataset_get_result_errors(self, ds, client, dataset_specs):
        import pandas as pd
        import management as mgt

        if dataset_specs is None:
            dataset_specs = ds.list_specifications().index.tolist()

        # gather results
        results = defaultdict(dict)
        all_res = list()
        for spec in dataset_specs:
            res = mgt.get_results(ds, spec, client)
            all_res.extend(res)

            for status in ["COMPLETE", "INCOMPLETE", "ERROR"]:
                results[spec][status] = len(
                    [r for r in res if r.status == status]
                )

        df = pd.DataFrame(results).transpose()
        df.index.name = "specification"
        return all_res, df

    def _errorcycle_dataset_report(self, df_res, res_error_counts):

        if len(res_error_counts) > 60000:
            res_error_counts = res_error_counts[:60000]
            res_error_counts += "\n--- Too many errors; truncated here ---\n"

        comment = f"""
        ## Lifecycle - Error Cycling Report

        {self._get_meta().to_markdown()}

        All errored tasks will be restarted.
        Errored states prior to restart reported below.

        ### `ResultRecord` current status

        {df_res.to_markdown()}

        #### `ResultRecord` Error Tracebacks:

        <details>
        <summary><b>Tracebacks</b> (<i>click to expand</i>)</summary>
        <!-- have to be followed by an empty line! -->

        ```
        {res_error_counts}
        ```
        </details>

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


class DataSet(SubmittableBase):
    """A dataset submitted to QCArchive.

    A dataset has a lifecycle with well-defined states.
    The state of a dataset is the state of its submission PR.

    """
    def submit(self, dataset_qcs, client):
        return dataset_qcs.submit(client=client, processes=1)


class Compute(SubmittableBase):
    """Supplemental compute submitted to QCArchive.

    """
    def submit(self, dataset_qcs, client):
        return dataset_qcs.submit(client=client, ignore_errors=True, processes=1)


def create_dataset(dataset_data):
    from openff.qcsubmit.datasets import BasicDataset, OptimizationDataset, TorsiondriveDataset

    datasets = {
        "dataset": BasicDataset,
        "optimizationdataset": OptimizationDataset,
        "torsiondrivedataset": TorsiondriveDataset,
    }

    if "type" in dataset_data:
        dataset_type = dataset_data["type"]
    elif "dataset_type" in dataset_data:
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
    packages = ["openff.qcsubmit", "openff.toolkit", "basis_set_exchange", "qcelemental"]
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
    parser.add_argument(
        "--set-priority",
        action='store_true',
        help="Triggers priority (re)setting based on Github PR label",
    )
    parser.add_argument(
        "--set-computetag",
        action='store_true',
        help="Triggers compute tag (re)setting based on Github PR label",
    )
    parser.add_argument(
        "--reset-errors",
        action='store_true',
        help="Whether to reset errored cases",
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

        # if we are setting priority, check for priority label(s) on PR
        # take highest one and set priority downstream
        # if no priority label(s), DO NOT set priority at all for this PR
        if args.set_priority:
            labels =  set(map(lambda x: x.name, pr.labels))
            priorities = set(PRIORITIES.keys()) & labels

            if not priorities:
                set_priority = False
                selected_priority = 1   # need something, but should have no effect due to `set_priority=False`
            else:
                set_priority = True
                selected_priority = 0
                for priority in priorities:
                    selected_priority = max(selected_priority, PRIORITIES[priority])

                print("Setting priority to '{}'".format(selected_priority))
        else:
            set_priority = False
            selected_priority = 1   # need something, but should have no effect due to `set_priority=False`

        if args.set_computetag:
            labels =  set(map(lambda x: x.name, pr.labels))
            computetags = [l[len('compute-'):] for l in labels if l.startswith('compute-')]

            if not computetags:
                set_computetag = False
                selected_computetag = 'openff'   # need something, but should have no effect due to `set_computetag=False`
            else:
                # if multiple compute tags on the PR, we choose the first one lexically
                set_computetag = True
                selected_computetag = sorted(computetags)[0]

                print("Setting computetag to '{}'".format(selected_computetag))
        else:
            set_computetag = False
            selected_computetag = 'openff'   # need something, but should have no effect due to `set_computetag=False`

        submission = Submission(pr, gh, priority=selected_priority, computetag=selected_computetag)
        submission.execute_state(board=board,
                                 states=states,
                                 reset_errors=args.reset_errors,
                                 set_priority=set_priority,
                                 set_computetag=set_computetag)

if __name__ == "__main__":
    main()
