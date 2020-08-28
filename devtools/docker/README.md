# QCFractal Dockerfiles

QCFractal Dockerfiles in this directory correspond to images provided on [Docker Hub](https://cloud.docker.com/u/openff/repository/list).


## `qcarchive_worker_openff`

This Dockerfile builds a container intended to be used as a compute worker.
It contains QCFractal as well as tools used by the [OpenFF](https://openforcefield.org/) workflow:

* [Psi4](http://www.psicode.org), [dftd3](https://github.com/loriab/dftd3), and [MP2D](https://github.com/Chandemonium/MP2D>)
* [RDKit](https://www.rdkit.org)
* [geomeTRIC](https://github.com/leeping/geomeTRIC)

including the following MM tools:
* [openforcefield](https://github.com/openforcefield/openforcefield)
* [openforcefields](https://github.com/openforcefield/openforcefields)
* [openmm](https://github.com/openmm/openmm)
* [openmmforcefields](https://github.com/openmm/openmmforcefields)

Its entrypoint launches a compute manager based on the configuration file, which must be provided at `/etc/qcfractal-manager/manager.yaml`.
