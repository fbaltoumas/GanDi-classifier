FROM condaforge/miniforge3:latest

LABEL description="GanDi-classifier: contig classifier for the Global Anaerobic Digestion (GanDi) database"

# External bioinformatics dependencies (skani, prodigal, prodigal-gv, diamond) via bioconda,
# git/pip needed to install gandi-classifier from source below, and C/C++ compilers in case
# scripts/check_cpu_baseline.py needs to rebuild numpy from source (see comment below) --
# numpy's build requires both.
RUN mamba install -y -c bioconda -c conda-forge \
        skani \
        prodigal \
        prodigal-gv \
        diamond \
        git \
        pip \
        c-compiler \
        cxx-compiler \
    && mamba clean -afy

# Building this image on a VM whose hypervisor masks CPU features (e.g. Hyper-V
# "Processor Compatibility Mode") hits the same numpy baseline RuntimeError at
# build time, since the build runs on the host's actual (masked) CPU. Run the
# CPU-baseline check before installing gandi-classifier so the correct numpy
# build is already in place when pip resolves its dependencies.
#
# TMPDIR is redirected off of /tmp: if /tmp is mounted noexec on the build host
# (common on hardened VMs), meson's compiler sanity check -- which compiles and
# then runs a tiny test binary under TMPDIR -- fails with a generic "Executables
# created by c compiler ... are not runnable" error during the numpy source
# rebuild below. /opt is on the same writable filesystem we're already using.
ENV TMPDIR=/opt/build-tmp
RUN mkdir -p "$TMPDIR"

RUN git clone --depth 1 https://github.com/fbaltoumas/GanDi-classifier.git /opt/GanDi-classifier \
    && python3 /opt/GanDi-classifier/scripts/check_cpu_baseline.py \
    && pip install --no-cache-dir /opt/GanDi-classifier \
    && rm -rf /opt/GanDi-classifier "$TMPDIR"

WORKDIR /data

ENTRYPOINT ["gandi-classifier"]
CMD ["--help"]
