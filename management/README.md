## GitHub Automation

For the GitHub automation to read and update the Dataset Tracking project board, it's necessary to make a Personal Access Token with the "resource owner" set to "openforcefield" available to the automation with at least 
* org:project:(read and write) access
* repo:pull request:(read and write) access

This token is accessed via the `QCA_DATASET_SUBMISSION_PAT` secret.

These tokens must be renewed periodically as GH only permits them a max lifetime of 1 year.

## QCArchive Credentials
There two sets of credentials that we use for adding / changing our data on QCArchive, one set for GitHub Actions
and another for OpenFF staff. These credentials can be changed at any time by logging in with the current credentials and running:
```python
client.change_user_password()
```
Where a new password will be printed. Be sure to update 1Password and GitHub secrets accordingly.
