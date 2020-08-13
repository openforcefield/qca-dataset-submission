# Production Environments

This file contains production environments for OpenFF QCArchive manager deployments.

Once environments are uploaded, they can be created from anywhere using the following line:
```bash
conda env create openforcefield/<env_name> -n <env_name>
```

To update an environment the current environment needs to be removed:
```bash
conda env remove -n <env_name>
```
It should be noted that this does not remove the packages, and creating a new environment should be very quick
as most required packages should already be installed.

To upload an environment:
```bash
anaconda -t <api_token> upload --user openforcefield <env_file.yaml>
```

