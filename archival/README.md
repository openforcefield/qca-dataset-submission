# Archival Support for Released Force Field Datasets

The quantum chemical (QC) records used to train each respective OpenFF force field from Sage 2.0.0 and later have been grouped into a single dataset on the QCArchive servers (https://api.qcarchive.molssi.org:443/) run by the Molecular Sciences Software Institute (MolSSI). This document provides a protocol to create archival storage of a force field dataset on Zenodo.

QCArchive stores quantum chemical (QC) data in a relational database in the SQLite format. Although a user would access it through the QCPortal python API, there is an option to create a SQLite "dataset view" which can be downloaded locally to explore, but cannot be further altered. We leverage this feature to store our complete training datasets in the common SQLite format on Zenodo for redundant archival.

To archive a dataset follow this protocol:

## Table of Contents

- [Create Dataset Views](#create-dataset-views)
- [Create a Docker Image](#create-the-docker-image)
- [Create a Zenodo Record](#create-a-zenodo-record)

## Create Dataset Views

```python
import os
from qcportal import PortalClient
client = PortalClient(
    "https://api.qcarchive.molssi.org:443/", 
    username="your_username",
    password="your_password",
    cache_dir=".",
)
os.makedirs("views", exist_ok=True)

ds_opt = client.get_dataset('optimization', 'OpenFF SMIRNOFF Sage 2.1.0')
internal_job_opt = ds_opt.create_view(
    description="Full Sage 2.1.0 Optimization dataset", 
    provenance={}, 
    include=['**'], 
    exclude=["wavefunction"], 
    include_children=True
)
ds_opt.download_view(
    destination_path=f"views/ds_{ds_opt.id}_{ds_opt.name.replace(' ','-')}_optimization_view.sqlite"
)

ds_td = client.get_dataset('torsiondrive', 'OpenFF SMIRNOFF Sage 2.1.0')
internal_job_td = ds_td.create_view(
    description="Full Sage 2.1.0 Torsiondrive dataset", 
    provenance={}, 
    include=['**'], 
    exclude=["trajectory", "compute_history", "outputs", "native_files", "wavefunction"],
    include_children=True
)
ds_td.download_view(
    destination_path=f"views/ds_{ds_td.id}_{ds_td.name.replace(' ','-')}_torsiondrive_view.sqlite"
)
```

These views will be included in the Zenodo entry. In order to access these files in the Docker image, put them in a directory called "views", the name of this directory can be adjusted in the `docker run` command.

## Docker Image

### Create the Docker Image

To create the Docker image the following files are needed:

- Dockerfile
- dataset_handling.ipynb
- environment.yaml

Compile the Docker image with:

`docker buildx build -t docker_handle_dataset_views .`

If using a M* MAC, consider using the flag `--platform=linux/amd64`. In this case, we create an image named, "docker_handle_dataset_views".

### Save the Docker Image to a File

Saving the Docker image to a file will allow it to be included in a Zenodo submission.

`docker save -o docker_handle_dataset_views.tar docker_handle_dataset_views:latest`

Optionally compress the file further:

`gzip docker_handle_dataset_views.tar`

### Loading a Compressed Docker Image File

Load the Docker image with:

`docker load -i docker_handle_dataset_views.tar.gz`

### Ensure the Docker Image Can Run and Load the Dataset

A user can run the Docker image to spawn the Jupyter notebook.

```bash
mkdir views; mv *sqlite views/
mkdir outputs
docker run -p 8888:8888 -v ./views:/workspace/views -v ./outputs:/workspace/outputs docker_handle_dataset_views
```

The `-p` flag exposed the port `8888` inside the Docker image to the port by the same name externally. 
The `-v` flag exposes a directory (in this case `./views`, so put your dataset views there) to a directory inside the Docker image so that the Jupyter notebook and access them.
The ./outputs directory provides another shared directory that can be useful to pass output files.
If using a M* MAC, the flag `--platform=linux/amd64` could be required.

Entering the URL that starts with `http://127.0.0.1:8888...` in a internet browser should lead to a JupyterLab instance.

Run the Jupyter notebook to ensure that:
 - The notebook runs completely without errors
 - The expected output files are produced by the notebook in the Docker container
 - The expected output can be accessed from the local `outputs` directory outside of the container

## Create a Zenodo Record

See our [Confluence page](https://openforcefield.atlassian.net/wiki/spaces/OFFO/pages/83951665/Zenodo) on the subject for more information.

- Title: "QC Fitting Datasets for *Force Field Dataset Name Here*"
- Authors: All dataset curators, generators, and associated PIs.
- Description: 
    - Copy from dataset description (details in the QCArchive datasets that comprise this overall dataset are removed and a referral is given to the repository)
    - Include instructions for running the Docker
    - Include descriptive statistics
- License: Creative Commons by 4.0
- Funding: Author dependent
- Journal: Publication Information (if the paper isn't out yet, skip this section)

Manifest:

- Optimization dataset view
- Torsiondrive dataset view
- Compressed Docker image file
- README.txt