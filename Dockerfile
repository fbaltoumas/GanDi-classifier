FROM condaforge/miniforge3:latest

LABEL description="GanDi-classifier: contig classifier for the Global Anaerobic Digestion (GanDi) database"

# External bioinformatics dependencies (skani, prodigal, prodigal-gv, diamond) via bioconda,
# plus git/pip needed to install gandi-classifier from source below.
RUN mamba install -y -c bioconda -c conda-forge \
        skani \
        prodigal \
        prodigal-gv \
        diamond \
        git \
        pip \
    && mamba clean -afy

# TODO: once published to PyPI, replace the git clone + pip install below with:
#   RUN pip install --no-cache-dir gandi
RUN git clone --depth 1 https://github.com/fbaltoumas/GanDi-classifier.git /opt/GanDi-classifier \
    && pip install --no-cache-dir /opt/GanDi-classifier \
    && rm -rf /opt/GanDi-classifier

WORKDIR /data

ENTRYPOINT ["gandi-classifier"]
CMD ["--help"]
