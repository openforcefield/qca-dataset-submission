
#!/usr/bin/env python

"""Lifecycle management for QCArchive datasets using GraphQL interface"""

import os
import pathlib
import requests
import textwrap


DATASET_GLOB = "dataset*.json*"
COMPUTE_GLOB = "compute*.json*"

PRIORITIES = {'priority-low': 0, 'priority-normal': 1, 'priority-high': 2}


def _post_query(query, variables=None):
    """Post a query to the GitHub GraphQL API"""
    headers = {"Authorization": f"Bearer {os.environ['GH_TOKEN']}"}
    json_dict = {"query": query}
    if variables:
        json_dict["variables"] = variables
    response = requests.post('https://api.github.com/graphql', json=json_dict, headers=headers)

    data = response.json()
    return data


class PullRequest:
    """
    A single pull request on a repository.

    Parameters
    ----------
    repo : Repo
        The repository where the PR is located
    id : str
        The node ID of the PR
    number : int
        The PR number
    title : str
        The PR title
    url : str
        The URL of the PR
    merged : bool, optional
        Whether the PR has been merged
    """
    def __init__(self, repo, id: str, number: int, title: str, url: str, merged=None):
        self.repo = repo
        self.id = id
        self.number = number
        self.title = title
        self.url = url
        self.merged = merged

    @classmethod
    def from_node_item(cls, repo, node_item):
        return cls(
            repo,
            node_item["id"],
            node_item["number"],
            node_item.get("title"),
            node_item.get("url"),
            node_item.get("merged"),
        )
    
    def get_label_names(self) -> list[str]:
        query = """
        query {
          repository(owner: "$owner", name: "$name") {
            pullRequest(number: $number) {
              labels(first: 10) {
                nodes {
                  name
                }
              }
            }
          }
        }
        """
        variables = {"owner": self.repo.owner, "name": self.repo.name, "number": self.number}
        data = _post_query(query, variables)
        label_names = []
        for node in data["data"]["repository"]["pullRequest"]["labels"]["nodes"]:
            label_names.append(node["name"])
        return label_names
    
    def add_to_labels(self, label: str):
        label_id = self.repo.get_label_id(label)
        query = """
        mutation {
          addLabelsToLabelable(input: {labelableId: "$id", labelIds: ["$label_id"]}) {
            labelable {
              id
            }
          }
        }
        """
        variables = {"id": self.id, "label_id": label_id}
        return _post_query(query, variables)
    
    def remove_from_labels(self, label: str):
        label_id = self.repo.get_label_id(label)
        query = """
        mutation {
          removeLabelsFromLabelable(input: {labelableId: "$id", labelIds: ["$label_id"]}) {
            labelable {
              id
            }
          }
        }
        """
        variables = {"id": self.id, "label_id": label_id}
        return _post_query(query, variables)
    
    def add_issue_comment(self, body: str):
        query = """
        mutation {
          addComment(input: {subjectId: "$id", body: "$body"}) {
            commentEdge {
              node {
                id
              }
            }
          }
        }
        """
        variables = {"id": self.id, "body": body}
        return _post_query(query, variables)
    
    def get_file_paths(self) -> list[pathlib.Path]:
        query = """
        query {
          repository(owner: "$owner", name: "$name") {
            pullRequest(number: $number) {
              files(first: 100) {
                nodes {
                  path
                }
              }
            }
          }
        }
        """
        variables = {"owner": self.repo.owner, "name": self.repo.name, "number": self.number}
        data = _post_query(query, variables)
        files = []
        for node in data["data"]["repository"]["pullRequest"]["files"]["nodes"]:
            files.append(pathlib.Path(node["path"]))
        return files


class Repo:
    def __init__(
        self,
        name: str = "qca-dataset-submission",
        owner: str = "openforcefield",
    ):
        self.name = name
        self.owner = owner
        self.repo_name = f"{owner}/{name}"

    def get_label_id(self, label: str):
        query = """
        query {
          repository(owner: "$owner", name: "$name") {
            label(name: "$label") {
              id
            }
          }
        }
        """
        variables = {"owner": self.owner, "name": self.name, "label": label}
        data = _post_query(query, variables)
        return data["data"]["repository"]["label"]["id"]
    
    def get_tracking_pull_requests(self) -> list[PullRequest]:
        """Get pull requests with the 'tracking' label"""

        query = """
        query {
          repository(owner: "$owner", name: "$name") {
            pullRequests(first: 100, labels: ["tracking"], after: $cursor) {
              pageInfo {
                hasNextPage
                endCursor
              }
              nodes {
                id
                number
                title
                url
              }
            }
          }
        }
        """

        variables = {"owner": self.owner, "name": self.name, "cursor": None}
        data = _post_query(query, variables)
        has_next_page = data["data"]["repository"]["pullRequests"]["pageInfo"]["hasNextPage"]

        prs = []
        for node_item in data["data"]["repository"]["pullRequests"]["nodes"]:
            pr = PullRequest.from_node_item(self, node_item)
            prs.append(pr)
        
        while has_next_page:
            cursor = data["data"]["repository"]["pullRequests"]["pageInfo"]["endCursor"]
            variables["cursor"] = cursor
            data = _post_query(query, variables)
            has_next_page = data["data"]["repository"]["pullRequests"]["pageInfo"]["hasNextPage"]
            for node_item in data["data"]["repository"]["pullRequests"]["nodes"]:
                pr = PullRequest.from_node_item(self, node_item)
                prs.append(pr)
        return prs
    
    def get_pull_request(self, number: int) -> PullRequest:
        """
        Get a pull request by number
        
        Parameters
        ----------
        number : int
            The PR number
            
        Returns
        -------
        PullRequest
        
        """
        query = """
        query {
          repository(owner: "$owner", name: "$name") {
            pullRequest(number: $number) {
              id
              title
              url
              merged
            }
          }
        }
        """
        variables = {"owner": self.owner, "name": self.name, "number": number}
        data = _post_query(query, variables)
        return PullRequest.from_node_item(self, data["data"]["repository"]["pullRequest"])


class ProjectV2PRCard:
    """
    A single card on a project board, corresponding to a single PR.

    Parameters
    ----------
    project : Project
        The project board where the card is located
    column : ProjectV2Column
        The column where the card is located
    card_node_id : str
        The node ID of the card
    card_url : str
        The URL of the card
    card_name : str
        The name of the card
    number : int
        The PR number
    """
    def __init__(self, project, column, card_node_id, card_url, card_name, number):
        self.project = project
        self.card_node_id = card_node_id
        self.card_url = card_url
        self.card_name = card_name
        self.column = column
        self.number = number

    def get_item(self) -> PullRequest:
        """Retrieve the PR associated with the card"""
        repo = self.project.repo
        return repo.get_pull_request(self.number)

class ProjectV2Column:
    """
    A single column on a project board.

    Parameters
    ----------
    project : Project
        The project board where the column is located
    column_node_id : str
        The node ID of the column
    column_name : str
        The name of the column


    Attributes
    ----------
    cards : list[ProjectV2PRCard]
        The cards in the column
    """
    def __init__(self, project, column_node_id, column_name):
        self.project = project
        self.column_node_id = column_node_id
        self.column_name = column_name
        self.cards = list()

    def add_card(self, item: PullRequest):
        """Add a card to the top of the specified column"""
        query = """
        mutation {
          addProjectCard(input: {contentId: "$content_id", projectColumnId: "$column_id"}) {
            cardEdge {
              node {
                id
                content {
                    __typename
                    ... on Issue {
                      title
                      url
                    }
                    ... on PullRequest {
                      title
                      url
                    }
                  }
              }
            }
          }
        }
        """
        variables = {
            "content_id": item.id,
            "column_id": self.column_node_id
        }
        data = _post_query(query, variables)
        return self._add_card_to_self_from_content(
            data["data"]["addProjectCard"]["cardEdge"]["node"]
        )
        

    def _add_card_to_self(self, card_node_id, card_url, card_name, card_number):
        """Updates self with card information"""
        card = ProjectV2PRCard(self.project, self, card_node_id, card_url, card_name, card_number)
        self.cards.append(card)
        self.project.cards_by_id[card_node_id] = card
        return card
    
    def _add_card_to_self_from_content(self, content):
        """Updates self with card information from content"""
        card_id = content['id']
        card_name = content['content']['title']
        card_url = content['content']['url']
        card_number = content['content'].get("number")
        return self._add_card_to_self(
            card_id, card_url, card_name, card_number
        )


class Project:
    """
    A single project board, corresponding to a single repository.

    Many assumptions are made as this is created for the Dataset Tracking
    board. This is not a general-purpose project board class.

    Parameters
    ----------
    repo : Repo
        The repository where the project board is located
    project_node_id : str
        The node ID of the project board
    """
    @classmethod
    def from_repo(cls, repo: Repo, project_number: int = 2):
        query = """
        query {
          organization(login: "$owner") {
          projectV2(number: $project_number) {
            id
            }
          }
        }
        """
        variables = {
            "owner": repo.owner,
            "project_number": project_number
        }
        data = _post_query(query, variables)
        project_node_id = data["data"]["repository"]["project"]["id"]
        return cls(repo, project_node_id)
    

    def _get_item_card(self, item: PullRequest):
        """
        Retrieve the card associated with an issue or PR. Currently only PRs supported
        """
        for card in self.cards_by_id.values():
            if card.number == item.number:
                return card
    

    def __init__(self, repo, node_id: str):
        self.repo = repo
        self.project_node_id = node_id
        self._reinitialize()


    def _reinitialize(self):
        self.columns_by_name = {}
        self.columns_by_id = {}
        self.cards_by_id = {}

        # set up project board
        project_data = self._get_project_data()
        # this is the card item
        for node_item in project_data:
            for field in node_item['fieldValues']['nodes']:
                if "name" in field:  # this is the column item
                    column_name = field['name']
                    column_node_id = field['id']
                    column = self.__create_or_retrieve_column(column_name, column_node_id)
                    column._add_card_to_self_from_content(node_item)

    def _create_or_retrieve_column(
        self,
        column_name: str,
        column_node_id: str,
    ):
        if column_name in self.columns_by_name:
            assert column_node_id in self.columns_by_id
            return self.columns_by_name[column_name]
        column = ProjectV2Column(self, column_node_id, column_name)
        self.columns_by_name[column_name] = column
        self.columns_by_id[column_node_id] = column
        return column
    
    

    

    def move_card_to_column(self, card, column: str):
        """Moves card to the top of the specified column"""
        if isinstance(card, str):
            card = self.cards_by_id[card]
        
        query = """
        mutation {
          moveProjectCard(input: {cardId: "$card_id", columnId: "$column_id"}) {
            cardEdge {
              node {
                id
              }
            }
          }
        }
        """
        variables = {
            "card_id": card.card_node_id,
            "column_id": self.columns_by_name[column].column_node_id
        }
        return _post_query(query, variables)



    def _get_project_data(self):
        query = """
        query {
          node(id: "$project_node_id") {
            ... on ProjectV2 {
              items(first: 100, after: $cursor) {
                nodes {
                  id
                  content {
                    __typename
                    ... on Issue {
                      title
                      url
                    }
                    ... on PullRequest {
                      title
                      url
                    }
                  }
                  fieldValues(first: 10) {
                    nodes {
                      ... on ProjectV2ItemFieldSingleSelectValue {
                        name
                        id
                      }
                    }
                  }
                }
                pageInfo {
                  hasNextPage
                  endCursor
                }
              }
            }
          }
        }
        """
        variables = {"project_node_id": self.project_node_id, "cursor": None}
        data = _post_query(query, variables)

        output_data = list(data['data']['node']['items']['nodes'])
        has_next_page = data["data"]["node"]["items"]["pageInfo"]["hasNextPage"]
        while has_next_page:
            cursor = data["data"]["node"]["items"]["pageInfo"]["endCursor"]
            variables["cursor"] = cursor
            data = _post_query(query, variables)
            output_data.extend(data['data']['node']['items']['nodes'])
            has_next_page = data["data"]["node"]["items"]["pageInfo"]["hasNextPage"]
        return output_data



class Submission:
    """A submission, corresponding to a single PR, possibly multiple datasets.

    A submission has a lifecycle with well-defined states.
    This class represents the current state of a submission,
    and provides the machinery for execution of lifecycle processes based on that state.

    All lifecycle state is stored on Github in the original PR for the submission,
    mapped onto states in the "Dataset Tracking" project board.

    Parameters
    ----------
    project : Project
        The project board where the submission is tracked
    item : PullRequest
        The PR corresponding to the submission
    priority : int, optional
        Priority to use for the dataset if set by method calls;
        one of 0, 1, or 2, in increasing-priority order.
    computetag : str, optional
        Compute tag to use for the dataset if set by method calls;
        tasks with a given compute tag will only be computed by managers
        configured to service that tag.

    """
    def __init__(self, project, item, priority=1, computetag="openff"):
        self.project = project
        self.item = item
        self.priority = priority
        self.computetag = computetag

        files = item.get_file_paths()
        self.computes = self._get_computes(files)
        self.datasets = self._get_datasets(files)

    
    def _get_datasets(self, files):
        datasets = []
        for file in files:
            if file.exists() and file.is_file() and file.match(DATASET_GLOB):
                datasets.append(str(file.resolve()))
        return datasets

    def _get_computes(self, files):
        computes = []
        for file in files:
            if file.exists() and file.is_file() and file.match(COMPUTE_GLOB):
                computes.append(str(file.resolve()))
        return computes



    def execute_state(self, states=None):
        card = self.project._get_item_card(self.item)
        # if card not on board, then it starts in the Backlog
        if card is None:
            column = self.project.columns_by_name["Backlog"]
            card = column.add_card(self.item)

            # reload board, since we just added this card
            self.project._reinitialize()
        
        # exit early if states specified, and this PR is not
        # in one of those
        if states is not None:
            if card.column.column_name not in states:
                return
            
        ACTIONS = {
            "Backlog": self.execute_backlog,
            "Queued for Submission": self.execute_queued_submit,
            "Error Cycling": self.execute_errorcycle,
            "Requires Scientific Review": self.execute_requires_scientific_review,
            "End of Life": self.execute_end_of_life,
            "Archived/Complete": self.execute_archived_complete,
        }
        if card.column.column_name in ACTIONS:
            return ACTIONS[card.column.column_name](card)
        
    
    def execute_backlog(self, card):
        item = card.get_item()
        if item.merged:
            comment = textwrap.dedent(
            """\
            ## Lifecycle - Backlog

            Merged dataset moved from "Backlog" to "Queued for Submission".

            """
            )
            item.add_issue_comment(comment)
            self.project.move_card_to_column(card, "Queued for Submission")

    def resolve_new_state(self, dataset_results) -> str:
        """If new state agreed upon by dataset results, that state is returned.
        Otherwise, returns `None`.
        """
        # get unique states recommended by datasets for this PR
        # may not always be the same, say, if e.g. submission fails for one
        # of many datasets in this submission
        new_state = None
        new_card_state = set(res["new_state"] for res in dataset_results)

        # if all datasets agree on the new card state, we change to that state
        if len(new_card_state) == 1:
            new_state = list(new_card_state)[0]
        return new_state

    def execute_queued_submit(self, card):
        from submittable import DataSet, Compute

        results = []
        for dataset in self.datasets:
            print(f"Processing dataset '{dataset}'")
            ds = DataSet(dataset, self)
            results.append(ds.submit())

        for compute in self.computes:
            print(f"Processing compute '{compute}'")
            ct = Compute(compute, self)
            results.append(ct.submit())
        
        new_state = self.resolve_new_state(results)
        if new_state is not None:
            self.project.move_card_to_column(card, new_state)

    def execute_errorcycle(self, card, reset_errors=False,
                           set_priority=False,
                           set_computetag=False):
        from submittable import DataSet, Compute

        results = []
        for dataset in self.datasets:
            print(f"Processing dataset '{dataset}'")
            ds = DataSet(
                dataset,
                self,
                priority=self.priority,
                computetag=self.computetag
            )
            results.append(ds.execute_errorcycle(reset_errors=reset_errors,
                                                 set_priority=set_priority,
                                                 set_computetag=set_computetag))

        for compute in self.computes:
            print(f"Processing compute '{compute}'")
            ct = Compute(
                compute,
                self,
                priority=self.priority,
                computetag=self.computetag
            )
            results.append(ct.execute_errorcycle(reset_errors=reset_errors,
                                                 set_priority=set_priority,
                                                 set_computetag=set_computetag))
        
        new_state = self.resolve_new_state(results)
        if new_state is not None:
            self.project.move_card_to_column(card, new_state)
        
        if new_state == "Archived/Complete":
            for dataset in self.datasets:
                ds = DataSet(dataset, self)
                ds.comment_archived_complete()

    def execute_requires_scientific_review(self, card):
        # add `scientific-review` label
        # remove `end-of-life`, `complete` label if present
        labels =  self.item.get_label_names()

        add_label = "scientific-review"

        if add_label not in labels:
            self.item.add_to_labels(add_label)

        for label in ("end-of-life", "complete"):
            if label in labels:
                self.item.remove_from_labels(label)

    
    def execute_end_of_life(self, card):
        # add `end-of-life` label
        # remove `scientific-review`, `complete` label if present
        labels =  self.item.get_label_names()

        add_label = "end-of-life"

        if add_label not in labels:
            self.item.add_to_labels(add_label)

        for label in ("scientific-review", "complete"):
            if label in labels:
                self.item.remove_from_labels(label)

    def execute_archived_complete(self, card):
        # add `complete` label
        # remove `scientific-review`, `end-of-life` label if present
        labels = self.item.get_label_names()

        add_label = "complete"

        if add_label not in labels:
            self.item.add_to_labels(add_label)

        for label in ("scientific-review", "end-of-life"):
            if label in labels:
                self.item.remove_from_labels(label)


    

def main():
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
    states = args.states if args.states else None
    prnums = args.prs if args.prs else None
    
    repo = Repo()

    # gather up all PRs with the `tracking` label
    prs = repo.get_tracking_pull_requests()
    if prnums is not None:
        prs = [
            pr.number
            for pr in prs
            if pr.number in prnums
        ]
    
    print(f"Found {len(prs)} with the 'tracking' label")
    # grab the full project board state once so we don't have to hammer the API
    project = Project.from_repo(repo, project_number=2)

    # for each PR, we examine the changes to find files used for the submission
    # this is where the mapping is made between the PR and the submission files
    for pr in prs:
        print(f"Processing PR #{pr.number}")

        labels = pr.get_label_names()

        # if we are setting priority, check for priority label(s) on PR
        # take highest one and set priority downstream
        # if no priority label(s), DO NOT set priority at all for this PR
        set_priority = False
        selected_priority = 1
        if args.set_priority:
            priorities = set(PRIORITIES.keys()).intersection(labels)

            if priorities:
                set_priority = True
                selected_priority = max([PRIORITIES[p] for p in priorities])
            print(f"Setting priority to {selected_priority} for PR #{pr.number}")
        
        set_computetag = False
        selected_computetag = "openff"
        if args.set_computetag:
            prefix = "compute-"
            n_prefix = len(prefix)
            computetags = [
                label[n_prefix:]
                for label in labels
                if label.startswith(prefix)
            ]

            if computetags:
                set_computetag = True
                # if multiple compute tags on the PR, we choose the first one lexically
                selected_computetag = sorted(computetags)[0]
            print(f"Setting compute tag to {selected_computetag} for PR #{pr.number}")
    
        submission = Submission(
            project,
            pr,
            priority=selected_priority,
            computetag=selected_computetag
        )
        submission.execute_state(
            states=states,
            reset_errors=args.reset_errors,
            set_priority=set_priority,
            set_computetag=set_computetag
        )
        