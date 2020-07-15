#!/usr/bin/env python

import os
import json
from collections import defaultdict
from datetime import datetime

from github import Github
import qcportal as ptl
import pandas as pd

import management as mgt

REPO_NAME = 'openforcefield/qca-dataset-submission'


class DataSet:
    """A dataset submitted to QCArchive.
    
    A dataset has a lifecycle with well-defined states.
    This class represents the current state of a dataset,
    and provides the machinery for execution of lifecycle processes based on that state.

    All lifecycle state is stored on Github in the original PR for the submission,
    mapped onto states in the "Datset Tracking" project board.
    
    """
    
    def __init__(self, directory, pr, ghapi, repo=None):
        """Create new DataSet instance linking a submission directory to its PR.

        Parameters
        ----------
        directory : path-like
            Path to dataset submission artifacts directory.
        pr : github.PullRequest
            PullRequest corresponding to the dataset submission.
        ghapi : github.Github
            An authenticated Github Python API client object.
        repo : str
            Github repo where datasets are tracked.

        """
        self.directory = os.path.abspath(directory)
        self.pr = pr
        self.ghapi = ghapi

        if repo is None:
            self.repo = ghapi.get_repo(REPO_NAME)
        else:
            self.repo = repo

    def _parse_spec(self):
        flagpost_file = 'spec.json'
        filepath = os.path.join(self.directory, flagpost_file)

        if os.path.exists(filepath):
            with open(filepath, 'r') as f:
                spec = json.load(f)

        dataset_name = spec['dataset_name']
        dataset_type = spec['dataset_type']

        return dataset_name, dataset_type

    def execute_state(self, board=None):
        """Based on current state of the PR, perform appropriate actions.

        """
        if board is None:
            board = _get_full_board(self.repo)

        # get lifecycle state, if it exists
        pr_state = None
        for state, cardconts in board.items():
            for cardcont in cardconts:
                if cardcont.number == self.pr.number:
                    pr_state = state
                    break

        # if card not on board, then it starts in the Backlog
        if pr_state is None:
            pr_state = self.set_backlog()
        
        if pr_state == "Backlog":
            self.execute_backlog()
        elif pr_state == "Queued for Submission":
            self.execute_queued()
        elif pr_state == "Error Cycling":
            return self.execute_errorcycle()
        elif pr_state == "Requires Scientific Review":
            self.execute_requires_scientific_review()
        elif pr_state == "End of Life":
            self.execute_end_of_life()
        elif pr_state == "Archived/Complete":
            self.execute_archived_complete()

    def _get_column(self, column):
        proj = [proj for proj in self.repo.get_projects()
                if proj.name == 'Dataset Tracking'][0]

        cols = list(proj.get_columns())
        return [col for col in cols if col.name == column][0]

    def set_backlog(self):
        backlog = self._get_column('Backlog')
        backlog.create_card(content_id=self.pr.id,
                            content_type='PullRequest')

    def execute_backlog(self):
        pass

    def execute_queued(self):
        """

        """

        result = self._submit_dataset()
        pass

    def _submit_dataset(self):
        """Use QCSubmit components to submit dataset to QCArchive for compute.

        """
        pass

    def execute_errorcycle(self, restart=False):

        dataset_name, dataset_type = self._parse_spec()

        client = ptl.FractalClient()
        ds = client.get_collection(dataset_type, dataset_name)
        
        if dataset_type == 'TorsionDriveDataset':
            return self._errorcycle_torsiondrive(ds, client)

        elif dataset_type == 'OptimizationDataset':
            return self._errorcycle_optimization(ds, client)

        elif dataset_type == 'GridOptimizationDataset':
            return self._errorcycle_gridopt(ds, client)

        elif dataset_type == 'Dataset':
            return self._errorcycle_dataset(ds, client)


    def _errorcycle_torsiondrive(self, ds, client):
        df_tdr = self._errorcycle_torsiondrive_get_tdr_errors(ds, client)
        df_tdr_opt = self._errorcycle_torsiondrive_get_tdr_opt_errors(ds, client)

        self._errorcycle_torsiondrive_report(df_tdr, df_tdr_opt)

        # restart errored torsiondrives and optimizations


    def _errorcycle_torsiondrive_get_tdr_errors(self, ds, client):
        # gather torsiondrive results
        results = defaultdict(dict)
        for spec in ds.list_specifications().index.tolist():
            tdrs = mgt.get_torsiondrives(ds, spec, client)

            for status in ['COMPLETE', 'RUNNING', 'ERROR']:
                results[spec][status] =  len(
                        [tdr for tdr in tdrs if tdr.status == status])
        
        df = pd.DataFrame(results).transpose()
        df.index.name = 'specification'
        return df

    def _errorcycle_torsiondrive_get_tdr_opt_errors(self, ds, client):
        # gather torsiondrive optimization results
        results = defaultdict(dict)
        for spec in ds.list_specifications().index.tolist():
            opts = mgt.merge(
                    mgt.get_torsiondrive_optimizations(ds, spec, client))

            for status in ['COMPLETE', 'INCOMPLETE', 'ERROR']:
                results[spec][status] =  len(
                        [opt for opt in opts if opt.status == status])
        
        df = pd.DataFrame(results).transpose()
        df.index.name = 'specification'
        return df

    def _errorcycle_torsiondrive_report(self, df_tdr, df_tdr_opt):
        datehr = datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC")

        dataset_name, dataset_type = self._parse_spec()

        meta = {'**Dataset Name**': dataset_name,
                '**Dataset Type**': dataset_type,
                '**UTC Date**': datehr}

        meta = pd.DataFrame(pd.Series(meta, name=""))

        # TODO: add error message samples for optimizations
        comment = f"""
        ## Error Cycling Report 

        {meta.to_markdown()}

        ### `TorsionDriveRecord` current status

        {df_tdr.to_markdown()}

        ### `OptimizationRecord` current status

        {df_tdr_opt.to_markdown()}

        """

        # postprocess due to raw spacing above
        comment = "\n".join([substr.strip() for substr in comment.split('\n')])

        # submit comment
        self.pr.create_issue_comment(comment)

    def _errorcycle_optimization(self, ds):
        pass

    def execute_requires_scientific_review(self):
        # check state
        ## exit early if already past this state
        pass

    def execute_end_of_life(self):
        # check state
        ## exit early if already past this state
        pass

    def execute_archived_complete(self):
        # check state
        ## exit early if already past this state
        pass


def _get_full_board(repo):
    proj = [proj for proj in repo.get_projects() if proj.name == 'Dataset Tracking'][0]
    board = {col.name: [card.get_content() for card in col.get_cards()] 
             for col in proj.get_columns()}
    return board

def _get_tracking_prs(repo):
    prs = [pr for pr in repo.get_pulls(state='all') 
       if 'tracking' in list(map(lambda x: x.name, pr.labels))]
    return prs

def _get_datadirs(files):
    """Given a file list from a PR, get the unique directories containing a
    flagpost file."""

    flagpost_file = 'spec.json'

    dirs = []
    for filepath in files:
        if flagpost_file in filepath.filename:
            dirs.append(os.path.dirname(filepath.filename))

    return sorted(set(dirs))

def main():
    """Map PRs tagged with 'tracking' into corresponding data directories.

    """
    import argparse
    
    parser = argparse.ArgumentParser(description='Process PRs according to dataset lifecycle')
    parser.add_argument('states', type=str, nargs='*',
                        help='states to process; if not provided, all states processed')
    
    args = parser.parse_args()

    gh = Github(os.environ['GH_TOKEN'])
    repo = gh.get_repo(REPO_NAME)

    # gather up all PRs with the `tracking` tag
    prs = _get_tracking_prs(repo)

    # grab the full project board state once so we don't have to hammer the API
    # over and over
    board = _get_full_board(repo)

    # if states provided, use it to filter our board
    # to only those states
    if args.states:
        board = {key: value for key, value in board.items() if key in args.states}

    # for each PR, we examine the changes to determine the directory for the submission
    # this is where the mapping is made between the PR and the submission files
    for pr in prs:
        
        files = pr.get_files()
        datadirs = _get_datadirs(files)

        # execute lifecycle process based on current state
        # TODO: adde excessive stdout logging for actions
        for datadir in datadirs:
            ds = DataSet(datadir, pr, gh)
            ds.execute_state(board)


if __name__ == '__main__':
    main()
