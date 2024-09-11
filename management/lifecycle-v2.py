
#!/usr/bin/env python

"""Lifecycle management for QCArchive datasets using GraphQL interface"""

import json
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
    def __repr__(self):
        return f"PullRequest({self.repo.repo_name}, {self.number}, {self.title})"

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
        """Get the names of the labels on the PR"""
        query = """
        query($owner: String!, $name: String!, $number: Int!) {
          repository(owner: $owner, name: $name) {
            pullRequest(number: $number) {
              labels(first: 100) {
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
        """Add a label to the PR by name

        Note: labels must already exist in the repository
        """
        label_id = self.repo.get_label_id(label)
        query = """
        mutation($id: ID!, $label_id: ID!) {
          addLabelsToLabelable(input: { labelableId: $id, labelIds: [$label_id] }) {
            labelable {
                labels {
                    nodes { name }
                }
            }
          }
        }
        """
        variables = {"id": self.id, "label_id": label_id}
        return _post_query(query, variables)
    
    def remove_from_labels(self, label: str):
        """Remove a label from the PR by name
        
        Note: labels must already exist in the repository
        """
        label_id = self.repo.get_label_id(label)
        query = """
        mutation($id: ID!, $label_id: ID!) {
          removeLabelsFromLabelable(input: { labelableId: $id, labelIds: [$label_id]}) {
            labelable {
              labels {
                nodes { name }
              }
            }
          }
        }
        """
        variables = {"id": self.id, "label_id": label_id}
        return _post_query(query, variables)
    
    def add_issue_comment(self, body: str):
        """Add a comment to the PR"""
        query = """
        mutation($id: ID!, $body: String!) {
          addComment(input: { subjectId: $id, body: $body }) {
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
        query($owner: String!, $name: String!, $number: Int!) {
          repository(owner: $owner, name: $name) {
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
    """A single repository on GitHub.
    
    Parameters
    ----------
    name : str, optional
        The name of the repository
    owner : str, optional
        The owner of the repository
    
    """
    def __init__(
        self,
        name: str = "qca-dataset-submission",
        owner: str = "openforcefield",
    ):
        self.name = name
        self.owner = owner
        self.repo_name = f"{owner}/{name}"

    def get_label_id(self, label: str) -> str:
        """Get the node ID of a label"""
        query = """
        query($owner: String!, $name: String!, $label: String!) {
          repository(owner: $owner, name: $name) {
            label(name: $label) {
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

        query_base = """
        query($owner: String!, $name: String! %s) {
          repository(owner: $owner, name: $name) {
            pullRequests(first: 100, labels: ["tracking"] %s) {
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
        query = query_base % ("", "")
        variables = {"owner": self.owner, "name": self.name}
        data = _post_query(query, variables)
        has_next_page = data["data"]["repository"]["pullRequests"]["pageInfo"]["hasNextPage"]

        prs = []
        for node_item in data["data"]["repository"]["pullRequests"]["nodes"]:
            pr = PullRequest.from_node_item(self, node_item)
            prs.append(pr)
        
        while has_next_page:
            query = query_base % (", $cursor: String", ", after: $cursor")
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
        query($owner: String!, $name: String!, $number: Int!) {
          repository(owner: $owner, name: $name) {
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
    def __repr__(self):
        return f"ProjectV2PRCard({self.project}, {self.column}, {self.card_node_id}, {self.card_name}, {self.number})"
    
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
    column_option_id : str
        The option ID of the column
    column_name : str
        The name of the column


    Attributes
    ----------
    cards : list[ProjectV2PRCard]
        The cards in the column
    """
    def __repr__(self):
        return f"ProjectV2Column({self.project}, {self.column_option_id}, {self.column_name})"

    def __init__(self, project, column_option_id, column_name):
        self.project = project
        self.column_option_id = column_option_id
        self.column_name = column_name
        self.cards = list()

    def add_card(self, item: PullRequest):
        """Add a card to the top of the specified column"""
        add_card_query = """
            mutation($project_id: ID!, $content_id: ID!) {
                addProjectV2ItemById(input: { projectId: $project_id, contentId: $content_id }) {
                    item {
                        id
                    }
                }
            }
        """

        variables = {
            "project_id": self.project.project_node_id,
            "content_id": item.id
        }

        data = _post_query(add_card_query, variables)
        card_id = data["data"]["addProjectV2ItemById"]["item"]["id"]

        card = self._add_card_to_self(
            card_id,
            item.url,
            item.title,
            item.number
        )
        self.project.move_card_to_column(card, self.column_name)


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
        query($owner: String!, $project_number: Int!) {
          organization(login: $owner) {
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
        project_node_id = data["data"]["organization"]["projectV2"]["id"]
        return cls(repo, project_node_id)
    
    def add_item_to_column(self, item: PullRequest, column: str):
        if isinstance(column, str):
            column = self.columns_by_name[column]

        return column.add_card(item)
    

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
        # The _column_status_id is the fieldId needed to label the Status
        self._column_status_id = ""
        self._reinitialize()
        


    def _reinitialize(self):
        self.columns_by_name = {}
        self.columns_by_id = {}
        self.cards_by_id = {}

        # set up project board
        # first get all columns and create them initially
        column_data = self._get_all_columns()
        for item in column_data:
            if item.get("name") == "Status":
                self._column_status_id = item["id"]

                for option in item["options"]:
                    self._create_or_retrieve_column(option["name"], option["id"])

        project_data = self._get_project_data()
        # this is the card item
        for node_item in project_data:
            for field in node_item['fieldValues']['nodes']:
                if "name" in field and field["field"]["name"] == "Status":  # this is the column item
                    column_name = field['name']
                    column_option_id = field['optionId']
                    column = self._create_or_retrieve_column(column_name, column_option_id)
                    column._add_card_to_self_from_content(node_item)


    def _create_or_retrieve_column(
        self,
        column_name: str,
        column_option_id: str,
    ):
        if column_name in self.columns_by_name:
            assert column_option_id in self.columns_by_id
            return self.columns_by_name[column_name]
        column = ProjectV2Column(self, column_option_id, column_name)
        self.columns_by_name[column_name] = column
        self.columns_by_id[column_option_id] = column
        return column
    

    def move_card_to_column(self, card, column: str):
        """Moves card to the top of the specified column"""
        if isinstance(card, str):
            card = self.cards_by_id[card]

        query = """
        mutation($card_id: ID!, $column_status_id: ID!, $project_id: ID!, $new_column_id: String!) {
          updateProjectV2ItemFieldValue(input: {
            itemId: $card_id,
            fieldId: $column_status_id,
            projectId: $project_id,
            value: { singleSelectOptionId: $new_column_id }
        }) {
            projectV2Item { id }
          }
        }
        """
        column = self.columns_by_name[column]
        variables = {
            "card_id": card.card_node_id,
            "column_status_id": self._column_status_id,
            "project_id": self.project_node_id,
            "new_column_id": column.column_option_id
        }
        return _post_query(query, variables)

    def _get_all_columns(self):
        """Get all columns, even if they're empty with no cards"""

        # 100 should be more. We can include pagination later if necessary

        query = """
        query($project_node_id: ID!) {
          node(id: $project_node_id) {
            ... on ProjectV2 {
              fields(first: 100) {
                nodes {
                  ... on ProjectV2SingleSelectField {
                    name
                    id
                    options {
                      id
                      name
                    }
                  }
                }
              }
            }
          }
        }
        """
        variables = {"project_node_id": self.project_node_id}
        data = _post_query(query, variables)
        return data["data"]["node"]["fields"]["nodes"]

    def _get_project_data(self):
        query_base = """
        query($project_node_id: ID! %s) {
          node(id: $project_node_id) {
            ... on ProjectV2 {
              items(first: 100 %s) {
                nodes {
                  id
                  content {
                    __typename
                    ... on Issue {
                      title
                      url
                      number
                    }
                    ... on PullRequest {
                      title
                      url
                      number
                    }
                  }
                  fieldValues(first: 10) {
                    nodes {
                      ... on ProjectV2ItemFieldSingleSelectValue {
                        field {
                          ... on ProjectV2SingleSelectField {
                            name
                            id
                          }
                        }
                        name
                        optionId
                        
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
        query = query_base % ("", "")
        variables = {"project_node_id": self.project_node_id}
        data = _post_query(query, variables)

        output_data = list(data['data']['node']['items']['nodes'])
        has_next_page = data["data"]["node"]["items"]["pageInfo"]["hasNextPage"]
        while has_next_page:
            query = query_base % (", $cursor: String", ", after: $cursor")
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
        """Process a PR in the 'Queued for Submission' state"""
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
        """Process a PR in the 'Error Cycling' state"""
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
        """Process a PR in the 'Requires Scientific Review' state"""
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
        """Process a PR in the 'End of Life' state"""
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
        """Process a PR in the 'Archived/Complete' state"""
        # add `complete` label
        # remove `scientific-review`, `end-of-life` label if present
        labels = self.item.get_label_names()

        add_label = "complete"

        if add_label not in labels:
            self.item.add_to_labels(add_label)

        for label in ("scientific-review", "end-of-life"):
            if label in labels:
                self.item.remove_from_labels(label)


def run_tests():
    repo = Repo()

    # gather up all PRs with the `tracking` label
    prs = repo.get_tracking_pull_requests()
    print(f"Found {len(prs)} with the 'tracking' label")
    print(prs)

    print("Creating project")
    project = Project.from_repo(repo, project_number=2)
    print(f"Project {project.project_node_id}")
    print(f"Columns: {project.columns_by_name.keys()}")
    for column_name, column in project.columns_by_name.items():
        print("===")
        print(column_name, column.column_option_id)
        for card in column.cards:
            print(card.card_name, card.card_node_id)
    
    print("===")
    print("Cards")
    for card in project.cards_by_id.values():
        print(card)
    print("")

    # try test on latest
    pr = sorted(prs, key=lambda pr: pr.number)[-1]
    print(f"Processing PR #{pr.number}")
    labels = pr.get_label_names()

    card = project._get_item_card(pr)
    previous_column = None
    if card:
        previous_column = card.column.column_name
    print(f"PR #{pr.number} is in column {previous_column}")

    # print files
    files = pr.get_file_paths()
    print("Files:")
    print(files)
    assert len(files) > 0

    # temporarily move to "Backlog"
    if card:
        print("Moving card to backlog")
        data = project.move_card_to_column(card, "Backlog")
    else:
        print(f"Adding card to backlog")
        data = project.add_item_to_column(pr, "Backlog")
        print(data)
    project._reinitialize()
    card = project._get_item_card(pr)
    assert card.column.column_name == "Backlog"

    # move back to original column
    if previous_column:
        project.move_card_to_column(card, previous_column)
        project._reinitialize()
        card = project._get_item_card(pr)
        assert card.column.column_name == previous_column


    # temporarily add label
    data = pr.add_to_labels("test-label")
    print(data)
    label_names = pr.get_label_names()
    assert "test-label" in label_names

    data = pr.remove_from_labels("test-label")
    print(data)
    labels = pr.get_label_names()
    assert "test-label" not in labels

    # add comment
    pr.add_issue_comment("This is a test comment from running CI. Please ignore")

    

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
    parser.add_argument(
        "--run-tests",
        action='store_true',
        help="Ignores everything else and runs a test function",
    )

    args = parser.parse_args()
    if args.run_tests:
        run_tests()
        return


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
        gc.collect()



if __name__ == "__main__":
    main()
