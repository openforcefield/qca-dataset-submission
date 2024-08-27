FROM condaforge/miniforge3

SHELL ["/bin/bash", "-c"]

RUN mamba env create -n qcfractal -f https://raw.githubusercontent.com/openforcefield/qca-dataset-submission/master/devtools/prod-envs/${ENV_NAME}.yaml

RUN groupadd -g 999 qcfractal && \
    useradd -m -r -u 999 -g qcfractal qcfractal

USER qcfractal

ENV PATH /opt/conda/bin:$PATH
RUN echo "source activate qcfractal" > ~/.bashrc

ENTRYPOINT /bin/bash -c "source activate qcfractal && qcfractal-compute-manager --config /etc/qcfractal-manager/manager.yaml --verbose"
