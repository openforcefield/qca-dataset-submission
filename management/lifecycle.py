
#!/usr/bin/env python

import os
import glob
import json
import re
import time
import typing
from pprint import pformat
import traceback
from itertools import chain
from collections import defaultdict, Counter
from datetime import datetime

from github import Github

if typing.TYPE_CHECKING:
    import qcelemental

QCFRACTAL_URL = "https://api.qcarchive.molssi.org:443/"

REPO_NAME = "openforcefield/qca-dataset-submission"
DATASET_GLOB = "dataset*.json*"
SCAFFOLD_GLOB = "scaffold*.json*"
COMPUTE_GLOB = "compute*.json*"

PRIORITIES = {'priority-low': 0, 'priority-normal': 1, 'priority-high': 2}

DATASET_TYPES = {
        'dataset': 'singlepoint',
        'optimizationdataset': 'optimization',
        'torsiondrivedataset': 'torsiondrive'}

# matches tags with a trailing _mw-###, allowing splitting of datasets by
# molecular weight based on these tags
SPLIT_TAG = re.compile(r"_mw(-\d+)+$")
# matches just the end of a SPLIT_TAG, capturing the final number
SPLIT_TAG_END = re.compile(r"-(\d+)$")


def parse_tags(compute_tag) -> tuple[list[float], str]:
    """Parses a compute tag matching ``SPLIT_TAG`` into a sequence of molecular
    weights. Also returns the base component of the tag"""
    tag = compute_tag
    ret = list()
    while (m := SPLIT_TAG_END.search(tag)) is not None:
        ret.append(float(m[1]))
        tag = tag[: m.start(1) - 1]

    # don't change the tag if it just happens to end with _mw, SPLIT_TAG only
    # matches if there are -### following it
    return list(reversed(ret)), tag.removesuffix("_mw") if len(ret) > 0 else tag


def try_get_molecule(entry) -> typing.Optional["qcelemental.models.Molecule"]:
    """Try to extract a qcelemental Molecule from multiple fields in ``entry``."""
    molecule = (
        # this should work for singlepoints
        getattr(entry, "molecule", None)
        # this for optimizations
        or getattr(entry, "initial_molecule", None)
    )

    if molecule:
        return molecule

    # and this for torsiondrives
    if (mols := getattr(entry, "initial_molecules")) is None or len(mols) < 1:
        return None

    return mols[0]


def partition_records(
    ds, bins, include_complete=False
) -> dict[int, list[int]]:
    """Split up the records in ``ds`` based on the molecular weights (in Da) in
    ``bins``.

    ``include_complete`` is intended mostly (only?) for testing; usually you
    wouldn't want to retag complete records (and it might be an error in
    qcportal since complete records don't have tags), but for tests it's nice
    to be able to call this on a finished dataset
    """
    import numpy as np

    ds.fetch_entries()
    masses, qca_ids = list(), list()
    for entry_name, _, rec in ds.iterate_records():
        if rec.status == "complete" and not include_complete:
            continue
        entry = ds.get_entry(entry_name)
        if (mol := try_get_molecule(entry)) is None:
            print(f"failed to get molecule from {entry_name}")
            continue

        masses.append(sum(mol.masses))
        qca_ids.append(rec.id)

    # TODO: Change so that numpy ints aren't used in the first place
    # For some reason QCPortal has an issue with numpy ints and we
    # must use python ints
    qca_ids = np.array(qca_ids)
    bin_indices = np.digitize(masses, bins)
    return {
            i: [int(x) for x in qca_ids[np.where(bin_indices == i)]]
            for i in range(len(bins) + 1)
    }


def set_mw_compute_tags(client, ds, compute_tag, include_complete=False):
    bins, base_tag = parse_tags(compute_tag)
    if len(bins) == 0:
        print(f"Failed to parse molecular weight compute tags from {compute_tag}")
        # TODO should this fall back on the normal setter then? we'd need to
        # pass more arguments to make that possible
        return

    records = partition_records(ds, bins, include_complete=include_complete)
    for bin_, record_ids in records.items():
        # the largest index may be 1 past len(bins) so just call this large
        suffix = int(bins[bin_]) if bin_ < len(bins) else "large"
        new_tag = f"{base_tag}-{suffix}"
        client.modify_records(record_ids, new_tag=new_tag)


def update_compute_tags(client, dataset, specification_names, new_tag, include_complete=False):
    """Update the compute tags in ``dataset`` to ``new_tag``, unless the new
    tag matches the ``SPLIT_TAG`` pattern, in which case the dataset will be
    split up and tagged separately based on molecular weight. For example,
    ``compute-openff_mw-100-200-300`` will cause the creation of four tags:
    ``compute-openff-100`` for MW < 100, ``compute-openff-200`` for 100 <= MW <
    200, ``compute-openff-300`` for 200 <= MW < 300, and
    ``compute-openff-large`` for anything larger than 300 Da.

    Note that ``set_mw_compute_tags`` does not need access to the specification
    names because it calls ``PortalClient.modify_records``, which, despite the
    identical name, is a separate method from ``BaseDataset.modify_records``.
    The client version used by ``set_mw_compute_tags`` relies on the record IDs
    to specify records instead of the specification name.

    Note also that ``set_mw_compute_tags`` should only be called when
    ``SPLIT_TAG`` matches the tag. It will print a warning and return early,
    updating no tags, if this is not the case.
    """
    if SPLIT_TAG.search(new_tag) is None:
        dataset.modify_records(
            specification_names=specification_names,
            new_tag=new_tag,
        )
    else:
        set_mw_compute_tags(client, dataset, new_tag, include_complete=include_complete)

def _get_labels(pr):
    return [ label.name for label in pr.get_labels() ]

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
            lambda x: any(glob.fnmatch.fnmatch(os.path.basename(x), match) for match in [DATASET_GLOB, SCAFFOLD_GLOB]),
            map(lambda x: x.filename, files)
        ))

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
                if int(card.pr_number) == int(pr.number):
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
        import projectsv2

        # we're going to give up on evolving boards entirely for now
        if board is None:
            board = projectsv2._get_full_board()
        # look for the card
        pr_state = None
        pr_card = None
        for col_name, cards in board.items():
            for card in cards:
                if int(card.pr_number) == int(self.pr.number):
                    pr_card = card
                    pr_state = col_name

        # pr_card, pr_state = self._get_board_card_state(board, self.pr)
        # if card not on board, then it starts in the Backlog
        # skip this completely, we can't do it
        # if pr_state is None:
            # pr_state = self.set_backlog()

            # reload board, since we just added this card
            # board = _get_full_board(self.repo)
            # pr_card, pr_state = self._get_board_card_state(board, self.pr)

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

            Project boards are not working as expected.
            However, please consider this queued for submission.

            """

            # postprocess due to raw spacing above
            comment = "\n".join([substr.strip() for substr in comment.split("\n")])

            # submit comment
            self.pr.create_issue_comment(comment)

        #     self.evolve_state(pr_card, pr_state, "Queued for Submission")

        #     return {"new_state": "Queued for Submission"}
        # else:
        #     return {"new state": "Backlog"}

    def execute_queued_submit(self, pr_card, pr_state):
        """Submit datasets, perhaps with some retry logic.

        """
        results = []
        for dataset in self.datasets:
            if "scaffold" not in dataset:
                print(f"Processing dataset '{dataset}'")
                ds = DataSet(dataset, self, self.ghapi)
                results.append(ds.execute_queued_submit())

        for compute in self.computes:
            print(f"Processing compute '{compute}'")
            ct = Compute(compute, self, self.ghapi)
            results.append(ct.execute_queued_submit())

        new_state = self.resolve_new_state(results)
        # if new_state is not None:
        #     self.evolve_state(pr_card, pr_state, new_state)

        # comment status on PR
        if new_state != pr_state:
            self.pr.create_issue_comment(f"## Current status - {new_state}\n\n Consider manually moving this.")

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
        # if new_state is not None:
        #     self.evolve_state(pr_card, pr_state, new_state)

        # comment status on PR
        if new_state != pr_state:
            self.pr.create_issue_comment(f"## Current status - {new_state}\n\n Consider manually moving this.")

        if new_state == "Archived/Complete":
            for dataset in self.datasets:
                ds = DataSet(dataset, self, self.ghapi)
                ds.comment_archived_complete()

    def execute_requires_scientific_review(self, pr_card, pr_state):
        # add `scientific-review` label
        # remove `end-of-life`, `complete` label if present
        labels =  _get_labels(self.pr)

        add_label = "scientific-review"

        if add_label not in labels:
            self.pr.add_to_labels(add_label)

        for label in ("end-of-life", "complete"):
            if label in labels:
                self.pr.remove_from_labels(label)

    def execute_end_of_life(self, pr_card, pr_state):
        # add `end-of-life` label
        # remove `scientific-review`, `complete` label if present
        labels =  _get_labels(self.pr)

        add_label = "end-of-life"

        if add_label not in labels:
            self.pr.add_to_labels(add_label)

        for label in ("scientific-review", "complete"):
            if label in labels:
                self.pr.remove_from_labels(label)

    def execute_archived_complete(self, pr_card, pr_state):
        # add `complete` label
        # remove `scientific-review`, `end-of-life` label if present
        labels =  _get_labels(self.pr)

        add_label = "complete"

        if add_label not in labels:
            self.pr.add_to_labels(add_label)

        for label in ("scientific-review", "end-of-life"):
            if label in labels:
                self.pr.remove_from_labels(label)


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

        if "dataset_name" in spec: # with dataset*.json from QCSubmit
            dataset_name = spec["dataset_name"]
            if "type" in spec:
                dataset_type = DATASET_TYPES[spec["type"].lower()]
            elif "dataset_type" in spec:
                dataset_type = DATASET_TYPES[spec["dataset_type"].lower()]
            dataset_specs = spec.get("qc_specifications", None)
        else: # with scaffold.json
            dataset_name = spec["metadata"]["name"]
            dataset_type = spec["metadata"]["dataset_type"]
            dataset_specs = None # Will be pulled from ds from qcportal call anyway

        return dataset_name, dataset_type, dataset_specs

    def _load_submittable(self):
        from openff.qcsubmit.serializers import deserialize
        spec = deserialize(self.submittable) # Will function with scaffold too

        return spec

    def _get_qca_client(self):
        import qcportal as ptl

        client = ptl.PortalClient(
            address=QCFRACTAL_URL,
            username=os.environ["QCA_USER"],
            password=os.environ["QCA_KEY"]
        )

        return client

    def _get_meta(self):
        import pandas as pd

        datehr = datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC")
        dataset_name, dataset_type, _ = self._parse_spec()

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
        ds = client.get_dataset(dataset_type, dataset_name)

        if dataset_type == "torsiondrive":
            complete = self._errorcycle_torsiondrive(
                    ds, client, dataset_specs,
                    reset_errors=reset_errors, set_priority=set_priority,
                    set_computetag=set_computetag)

        elif dataset_type == "optimization":
            complete = self._errorcycle_dataset(
                    ds, client, dataset_specs,
                    self._errorcycle_optimization_report,
                    reset_errors=reset_errors, set_priority=set_priority,
                    set_computetag=set_computetag)

        elif dataset_type == "singlepoint":
            complete = self._errorcycle_dataset(
                    ds, client, dataset_specs,
                    self._errorcycle_dataset_report,
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

    @staticmethod
    def count_unique_error_messages(errors_in, pretty_print=False):
        errors = defaultdict(set)
    
        for id, error in errors_in.items():
            errors["\n".join([error[i] for i in ['error_type', 'error_message']])].add(id)
    
        errors = dict(errors)
    
        content = ""
        if pretty_print:
            for count, key, value in sorted([(len(value), key, value) for key, value in errors.items()], reverse=True):
                content += '-------------------------------------\n'
                content += f"count : {count}\n"
                content += '\n'
                content += f'{key}\n'
                content += '\n'
                content += 'ids : \n'
                content += f'{pformat(value, width=80, compact=True)}\n'
                content += '-------------------------------------\n'
            return content
        else:
            return errors

    def _errorcycle_torsiondrive(self, ds, client, dataset_specs,
            reset_errors=False, set_priority=False, set_computetag=False):
        import pandas as pd
        from qcportal.record_models import RecordStatusEnum

        if dataset_specs is None:
            dataset_specs = ds.specification_names

        df_status = self._errorcycle_get_status(ds, dataset_specs)

        if reset_errors:
            erred_rec_ids = []
            erred_opts = {}
            status_counts = {}
            for ds_spec in dataset_specs:
                recs = ds.iterate_records(
                        specification_names=[ds_spec], 
                        #status='error'
                        )

                # build up optimization statuses and errors, if present
                erred_opts[ds_spec] = []
                status_counts[ds_spec] = Counter({status.value.upper(): 0 for status in list(RecordStatusEnum)})
                for entry, spec, rec in recs:
                    if rec.status == 'error':
                        erred_rec_ids.append(rec.id)
                    for opt in chain.from_iterable(rec.optimizations.values()):
                        status_counts[ds_spec][opt.status.value.upper()] += 1

                        if opt.status == 'error':
                            erred_opts[ds_spec].append((opt.id, opt.error))

            # create status counts dataframe
            df_opt_status = pd.DataFrame(status_counts).transpose()
            df_opt_status = df_opt_status[['COMPLETE', 'RUNNING', 'WAITING', 'ERROR', 'CANCELLED', 'INVALID', 'DELETED']]
            df_opt_status.index.name = 'specification'

            # aggregate all errors to get single set of counts for error messages
            errors = {}
            for ds_spec in erred_opts:
                errors.update({r[0]: r[1] for r in erred_opts[ds_spec]})

            error_counts = self.count_unique_error_messages(errors, pretty_print=True)

            self._errorcycle_torsiondrive_report(df_status, df_opt_status, error_counts)

        if df_status[["WAITING", "RUNNING", "ERROR"]].sum().sum() == 0:
            complete = True
        else:
            if reset_errors:
                client.reset_records(erred_rec_ids)
            if set_priority:
                ds.modify_records(specification_names=list(dataset_specs),
                                  new_priority=self.priority)
            if set_computetag:
                update_compute_tags(
                    client=client,
                    dataset=ds,
                    specification_names=list(dataset_specs),
                    new_tag=self.computetag,
                )

            complete = False

        return complete

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

    def _errorcycle_get_status(self, ds, dataset_specs):
        import pandas as pd
        from qcportal.record_models import RecordStatusEnum

        if dataset_specs is None:
            dataset_specs = ds.specification_names
        
        status = ds.status()
        status_ = {key: {status.value.upper(): counts.get(status, 0)
                         for status in list(RecordStatusEnum)}
                   for key, counts in status.items() if key in dataset_specs.keys()}

        df = pd.DataFrame(status_).transpose()
        df = df[['COMPLETE', 'RUNNING', 'WAITING', 'ERROR', 'CANCELLED', 'INVALID', 'DELETED']]
        df.index.name = 'specification'

        return df

    def _errorcycle_optimization_report(self, df_status, opt_error_counts):

        if len(opt_error_counts) > 60000:
            opt_error_counts = opt_error_counts[:60000]
            opt_error_counts += "\n--- Too many errors; truncated here ---\n"

        comment = f"""
        ## Lifecycle - Error Cycling Report

        {self._get_meta().to_markdown()}

        All errored tasks will be restarted.
        Errored states prior to restart reported below.

        ### `OptimizationRecord` current status

        {df_status.to_markdown()}

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

    def _errorcycle_dataset(self, ds, client, dataset_specs, report_method,
            reset_errors=False, set_priority=False, set_computetag=False):

        if dataset_specs is None:
            dataset_specs = ds.specification_names

        df_status = self._errorcycle_get_status(ds, dataset_specs)

        if reset_errors:
            erred_recs = ds.iterate_records(
                    specification_names=list(dataset_specs), 
                    status='error')

            errors = {r.id: r.error for entry, spec, r in erred_recs}
            error_counts = self.count_unique_error_messages(errors, pretty_print=True)

            report_method(df_status, error_counts)

        if df_status[["WAITING", "RUNNING", "ERROR"]].sum().sum() == 0:
            complete = True
        else:
            if reset_errors:
                client.reset_records(list(errors))
            if set_priority:
                ds.modify_records(specification_names=list(dataset_specs),
                                  new_priority=self.priority)
            if set_computetag:
                update_compute_tags(
                    client=client,
                    dataset=ds,
                    specification_names=list(dataset_specs),
                    new_tag=self.computetag,
                )
            complete = False

        return complete

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

    def submit(self, dataset_qcs, client):
        return dataset_qcs.submit(client=client, ignore_errors=True)


class DataSet(SubmittableBase):
    """A dataset submitted to QCArchive.

    A dataset has a lifecycle with well-defined states.
    The state of a dataset is the state of its submission PR.

    """
    ...


class Compute(SubmittableBase):
    """Supplemental compute submitted to QCArchive.

    """
    ...


def create_dataset(dataset_data):
    from openff.qcsubmit.datasets import BasicDataset, OptimizationDataset, TorsiondriveDataset

    datasets = {
        "dataset": BasicDataset,
        "optimizationdataset": OptimizationDataset,
        "torsiondrivedataset": TorsiondriveDataset,
    }

    if "type" in dataset_data:
        dataset_type = dataset_data["type"].lower()
    elif "dataset_type" in dataset_data:
        dataset_type = dataset_data["dataset_type"].lower()

    dataset_class = datasets.get(dataset_type, None)
    if dataset_class is not None:
        return dataset_class.parse_obj(dataset_data)
    else:
        raise RuntimeError(f"The dataset type {dataset_type} is not supported.")


# def _get_full_board(repo):
#     proj = [proj for proj in repo.get_projects() if proj.name == "Dataset Tracking"][0]
#     board = {col.name: [card for card in col.get_cards()] for col in proj.get_columns()}
#
#     # attach pr number to each card; we do this *once* here to avoid too many API calls,
#     # exhausting our limit
#     for col, cards in board.items():
#         for card in cards:
#             card.pr_number = card.get_content().number
#
#     return board


def _get_tracking_prs(repo):
    prs = [
        pr
        for pr in repo.get_pulls(state="all")
        if "tracking" in _get_labels(pr)
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
    import gc

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
    #board = _get_full_board(repo)
    import projectsv2
    board = projectsv2._get_full_board()

    # for each PR, we examine the changes to find files used for the submission
    # this is where the mapping is made between the PR and the submission files
    for pr in prs:

        print(f"Processing PR #{pr.number}")

        # if we are setting priority, check for priority label(s) on PR
        # take highest one and set priority downstream
        # if no priority label(s), DO NOT set priority at all for this PR
        if args.set_priority:
            labels =  _get_labels(pr)
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
            labels =  _get_labels(pr)
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

        gc.collect()

if __name__ == "__main__":
    main()
