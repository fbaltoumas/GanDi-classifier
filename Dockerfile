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

# C/C++ compilers in case scripts/check_cpu_baseline.py needs to rebuild numpy from
# source (see comment below). Deliberately the base image's own system gcc/g++ via
# apt, not conda-forge's c-compiler/cxx-compiler packages: those bundle their own
# sysroot/glibc built by conda-forge's own CI, which is a plausible explanation for
# meson's compiler sanity check failing with "Executables created by c compiler
# ... are not runnable" on a VM whose hypervisor masks CPU features (glibc's
# CPU-dispatched routines, e.g. memcpy/strlen IFUNC resolvers, could behave
# differently between that bundled glibc and the CPU than the base image's own
# glibc) -- this was tried after ruling out compiler-default -march flags and a
# /tmp noexec mount as the cause. The system compiler is linked against the same
# glibc already used everywhere else in this image, avoiding that mismatch.
# CC/CXX are set explicitly so pip's build backend picks these up regardless of
# what else ends up on PATH.
RUN apt-get update \
    && apt-get install -y --no-install-recommends build-essential \
    && rm -rf /var/lib/apt/lists/*
ENV CC=gcc
ENV CXX=g++

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
