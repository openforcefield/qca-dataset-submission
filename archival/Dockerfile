FROM condaforge/miniforge3

SHELL ["/bin/bash", "-c"]

ADD environment.yaml /tmp/env.yaml
ADD dataset_handling.ipynb /workspace/dataset_handling.ipynb

RUN mamba env create -n qca -f /tmp/env.yaml && \
    mamba shell init --shell bash && \
    eval "$(mamba shell hook --shell bash)"

WORKDIR /workspace

ENV CONDA_DEFAULT_ENV=qca
ENV CONDA_PREFIX=/opt/conda/envs/qca
ENV PATH=/opt/conda/envs/qca/bin:$PATH

# Expose Jupyter port
EXPOSE 8888

# Start Jupyter Lab with the notebook
CMD ["jupyter", "lab", "dataset_handling.ipynb", "--ip=0.0.0.0", "--no-browser", "--allow-root"]