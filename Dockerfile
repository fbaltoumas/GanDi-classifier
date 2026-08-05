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

# Building this image on a VM whose hypervisor masks CPU features (e.g. Hyper-V
# "Processor Compatibility Mode") hits the same numpy baseline RuntimeError at
# build time, since the build runs on the host's actual (masked) CPU. Run the
# CPU-baseline check before installing gandi-classifier so the correct numpy
# build is already in place when pip resolves its dependencies.
RUN git clone --depth 1 https://github.com/fbaltoumas/GanDi-classifier.git /opt/GanDi-classifier \
    && python3 /opt/GanDi-classifier/scripts/check_cpu_baseline.py \
    && pip install --no-cache-dir /opt/GanDi-classifier \
    && rm -rf /opt/GanDi-classifier

WORKDIR /data

ENTRYPOINT ["gandi-classifier"]
CMD ["--help"]
