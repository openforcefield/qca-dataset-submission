FROM condaforge/miniforge3

SHELL ["/bin/bash", "-c"]

RUN mamba install anaconda-client
RUN mamba env create -n qcfractal openforcefield/${ENV_NAME}

RUN groupadd -g 999 qcfractal && \
    useradd -m -r -u 999 -g qcfractal qcfractal

USER qcfractal

ENV PATH /opt/conda/bin:$PATH
RUN echo "source activate qcfractal" > ~/.bashrc

ENTRYPOINT /bin/bash -c "source activate qcfractal && qcfractal-compute-manager --config /etc/qcfractal-manager/manager.yaml --verbose"
