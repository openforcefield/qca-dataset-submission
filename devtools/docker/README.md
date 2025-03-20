# QCFractal Dockerfiles

QCFractal Dockerfiles in this directory correspond to images provided on the [GitHub Container Registry (GHCR) associated with this repo](https://github.com/openforcefield/qca-dataset-submission/pkgs/container/qca-dataset-submission/versions).


## `qcarchive_worker_openff`

This Dockerfile builds a container intended to be used as a compute worker.

Its entrypoint launches a compute manager based on the configuration file, which must be provided at `/etc/qcfractal-manager/manager.yaml`.


## Updating worker images 

* Do the following on the `master` branch:
    * Change the appropriate conda yaml in [this folder](https://github.com/openforcefield/qca-dataset-submission/tree/master/devtools/prod-envs). The filenames map to the images in the [container registry](https://github.com/openforcefield/qca-dataset-submission/pkgs/container/qca-dataset-submission/versions).
    * Uncomment only the env_names corresponding to the images you want to rebuild in the [docker image building action](https://github.com/openforcefield/qca-dataset-submission/blob/master/.github/workflows/prod-envs-docker.yml).
* Run the container building action [here](https://github.com/openforcefield/qca-dataset-submission/actions/workflows/prod-envs-docker.yml)
* IF you're deploying this on NRP
    * Shut down your previous deployments using this image with ex. `kubectl delete`
    * Update your `deployment.yaml` (NOT the template on GitHub, but your local one) to force use of the new image (change the `-image` line to `- image: "ghcr.io/openforcefield/qca-dataset-submission:qcarchive-worker-openff-xtb-latest@sha256:3bae3b16ec0b34d0bc0b9359961bea715eff972952d54364166d6d38bdf2ce57"`)
    * Re-launch your deployment with ex. `kubectl apply`
 

### Debugging
* To interactively run the docker image locally: `docker run -it --pull always --entrypoint /bin/bash ghcr.io/openforcefield/qca-dataset-submission:qcarchive-worker-openff-<ENV_NAME>-latest`
    * The `--entrypoint` is needed because using the default entrypoint will try to spin up a compute manager and immediately die since it doesn't have access to the `manager.yaml`.
    * The `--pull always` is needed because otherwise it will be really lazy about fetching the latest version of an image.
    * If it still doesn't fetch the latest/desired image, you can try forcing it by specifying a specific container hash, which you can get from the [container registry page](https://github.com/openforcefield/qca-dataset-submission/pkgs/container/qca-dataset-submission/versions) ex: `docker run -it --entrypoint /bin/bash ghcr.io/openforcefield/qca-dataset-submission:qcarchive-worker-openff-xtb-latest:sha256:ed56d81f6a48e195832f7b44e74f6042b787e9535ba3af7e7d44be7f5d89b7af`
* Ensure that the env files are really being pulled from `master` [here](https://github.com/openforcefield/qca-dataset-submission/blob/master/devtools/docker/qcarchive-worker-openff/Dockerfile#L5)

