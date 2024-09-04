For the GitHub automation to read and update the Dataset Tracking project board, it's necessary to make a token available to the automation with at least 
* org:project:(read and write) access
* repo:pull request:(read and write) access

This token is accessed via the `QCA_DATASET_SUBMISSION_PAT` secret.

These tokens must be renewed periodically as GH only gives them a max lifetime of 1 year.
