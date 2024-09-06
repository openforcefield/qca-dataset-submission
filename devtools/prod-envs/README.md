# Production Environments

This file contains production environments for OpenFF QCArchive manager deployments.

Once environments are uploaded, they can be created from anywhere using the following line:
```bash
mamba env create -n <env_name> -f https://raw.githubusercontent.com/openforcefield/qca-dataset-submission/use-ghcr/devtools/prod-envs/qcarchive-worker-openff-psi4.yaml 
```

To update an environment the current environment needs to be removed:
```bash
mamba env remove -n <env_name>
```
It should be noted that this does not remove the packages, and creating a new environment should be very quick
as most required packages should already be installed.

To modify an environment, open a pull request to this repository and modify these files directly. Newly started workers will immediately begin using the updated yamls. 

To make new docker images, run the [Deployment - QCArchive Prod Docker Environments action](https://github.com/openforcefield/qca-dataset-submission/actions/workflows/prod-envs-docker.yml). The new images will be served at https://github.com/openforcefield/qca-dataset-submission/pkgs/container/qca-dataset-submission.  