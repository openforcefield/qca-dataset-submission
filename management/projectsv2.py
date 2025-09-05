import os
import requests


class ProjectV2Project(dict):
    def __init__(self, project_node_id):
        self.project_node_id = project_node_id
        data = self._get_project_data(project_node_id)
        print("Got Data")
#        print(data) # I think this causes execute_submission() to occur here for some reason
        self.columns = dict()
        print("Made it")
        for item in data['data']['node']['items']['nodes']:
            for field in item['fieldValues']['nodes']:
                if 'name' in field:
                    column_name = field['name']
                    if column_name not in self.columns:
                        self.columns[column_name] = ProjectV2Column(self,
                                                                    field['id'],
                                                                    column_name
                                                                    )
                    self.columns[column_name].cards.append(ProjectV2PRCard(self,
                                                                           self.columns[column_name],
                                                                           item['id'],
                                                                           item['content']['url'],
                                                                           item['content']['title']
                                                                           ))
                    print(column_name, field['id'])

    def _get_project_data(self, project_node_id):
        query = """
        query {
          node(id: "%s") {
            ... on ProjectV2 {
              items(first: 100) {
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
        """ % project_node_id

        headers = {"Authorization": f"Bearer {os.environ['GH_TOKEN']}"}
        response = requests.post('https://api.github.com/graphql', json={'query': query}, headers=headers)

        data = response.json()

        return data


class ProjectV2Column:
    def __init__(self, project, column_node_id, column_name):
        self.project = project
        self.column_node_id = column_node_id
        self.column_name = column_name
        self.cards = list()


class ProjectV2PRCard:
    def __init__(self, project, column, card_node_id, card_url, card_name):
        self.project = project
        self.card_node_id = card_node_id
        self.card_url = card_url
        self.card_name = card_name


    def move(position, column):
        pass


def _get_full_board():
    proj = ProjectV2Project("PVT_kwDOARrkss4Am84U")
    board = {col.column_name: [card for card in col.cards] for col in proj.columns.values()}

    for col, cards in board.items():
        for card in cards:
            card.pr_number = card.card_url.split('/')[-1]

    return board
