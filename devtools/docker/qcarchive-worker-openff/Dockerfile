FROM condaforge/miniforge3

SHELL ["/bin/bash", "-c"]

RUN mamba env create -n qcfractal -f https://raw.githubusercontent.com/openforcefield/qca-dataset-submission/master/devtools/prod-envs/${ENV_NAME}.yaml && \
    if echo "$ENV_NAME" | grep -q "qcarchive-worker-openff-psi4-basis-set-update"; then \
        mamba run -n qcfractal pip install basis_set_exchange && \
        mamba run -n qcfractal python -c "import basis_set_exchange as bse, os; \
            conda_prefix = os.environ.get('CONDA_PREFIX', '/opt/conda/envs/qcfractal'); \
            path = os.path.join(conda_prefix, 'share', 'psi4', 'basis', 'aug-cc-pvtz-dk.gbs'); \
            data = bse.get_basis('aug-cc-pVTZ-dk', fmt='psi4', header=True, optimize_general=True); \
            open(path, 'w').write(data)"; \
    fi

RUN groupadd -g 999 qcfractal && \
    useradd -m -r -u 999 -g qcfractal qcfractal

USER qcfractal

ENV PATH=/opt/conda/bin:$PATH
RUN echo "source activate qcfractal" > ~/.bashrc

ENTRYPOINT /bin/bash -c "source activate qcfractal && qcfractal-compute-manager --config /etc/qcfractal-manager/manager.yaml --verbose"