# Force field Release Dataset Archival

The quantum chemical (QC) records used to train each respective OpenFF force field from Sage 2.0.0 and later have been grouped into a single dataset on the QCArchive servers (https://api.qcarchive.molssi.org:443/) run by the Molecular Sciences Software Institute (MolSSI). This document provides a protocol to create archival storage of a force field dataset on Zenodo.

QCArchive stores quantum chemical (QC) in a relational database in the SQLite format. Although a user would access it through the QCPortal python API, there is an option to create a SQLite "dataset view" which can be downloaded locally to explore, but cannot be further altered. We leverage this feature to store our complete training datasets in the common SQLite format on Zenodo for redundant archival.

To archive a dataset follow this protocol:

## Table of Contents

- [Create Dataset Views](#create-dataset-views)
- [Create a Docker File](#create-a-docker-file)
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
    destination_path=f"views/ds_{ds.id}_{ds.name.replace(' ','-')}_optimization_view.sqlite"
)

ds_td = client.get_dataset('torsiondrive', 'OpenFF SMIRNOFF Sage 2.1.0')
internal_job_td = ds_td.create_view(
    description="Full Sage 2.1.0 Torsiondrive dataset", 
    provenance={}, 
    include=['**'], 
    exclude=["wavefunction"], 
    include_children=True
)
ds_td.download_view(
    destination_path=f"views/ds_{ds.id}_{ds.name.replace(' ','-')}_torsiondrive_view.sqlite"
)
```

## Create a Docker File

```docker
```

Where the conda environment yaml could contain:
```yaml
name: qca
channels:
  - conda-forge/label/libint_dev
  - conda-forge
dependencies:
  - python=3.11
  - pip
  - qcfractalcompute>=0.61
  - qcengine
  - qcelemental
```
Commands to compile docker image

## Create a Zenodo Record

See our [Confluence page](https://openforcefield.atlassian.net/wiki/spaces/OFFO/pages/83951665/Zenodo) on the subject for more information.

- Title: Same as force field
- Authors: All submitters, curators, and generators
- Description: Copy from release notes
- License: Creative Commons by 4.0
- Funding: Same as publication
- Keywords: The same as the publication
- Related Identifiers: GitHub repository link
- Journal: Publication Information
