# GanDi-classifier
Contig classifier for the Global Anaerobic Digestion (GanDi) database

## Table of Contents

- [Dependencies](#dependencies)
  - [External executables](#external-executables)
  - [Python packages](#python-packages)
- [Installation](#installation)
  - [a) Create a conda environment](#a-create-a-conda-environment)
  - [b) Install the external dependencies through conda](#b-install-the-external-dependencies-through-conda)
  - [c) Install the `gandi` package](#c-install-the-gandi-package)
  - [Windows users](#windows-users)
  - [Docker](#docker)
- [Usage](#usage)
  - [Overview](#overview)
  - [Database download](#database-download)
  - [Classification](#classification)
  - [Examples](#examples)
- [Troubleshooting](#troubleshooting)
  - [`RuntimeError: NumPy was built with baseline optimizations... but your machine doesn't support...`](#runtimeerror-numpy-was-built-with-baseline-optimizations-but-your-machine-doesnt-support)
  - [`Illegal instruction (core dumped)` after importing polars](#illegal-instruction-core-dumped-after-importing-polars)

## Dependencies

`gandi` requires **Python 3.11+**.

### External executables

`gandi classifier` shells out to `skani`, `prodigal`/`prodigal-gv`, and `diamond`; `gandi download-db` shells out to `skani` and `diamond`. These are not Python packages and are not installed automatically — they need to be installed separately and be available on `PATH`. See [Installation](#installation) for how to install them via conda.

- **[skani](https://github.com/bluenote-1215/skani)** (ANI workflow)

  Alternatively to conda, download a prebuilt binary from the [skani releases page](https://github.com/bluenote-1215/skani/releases), or install via `cargo install skani` if you have a Rust toolchain.

- **[prodigal](https://github.com/hyattpd/Prodigal)** (AAI workflow)

  On Debian/Ubuntu, it's also available via `apt`:
  ```bash
  sudo apt install prodigal
  ```

- **[prodigal-gv](https://github.com/apcamargo/prodigal-gv)** (AAI workflow, used automatically by `gandi classifier` when `-c viruses` is given)

- **[diamond](https://github.com/bbuchfink/diamond)** (AAI workflow)

  Alternatively to conda, download a prebuilt binary from the [diamond releases page](https://github.com/bbuchfink/diamond/releases).

### Python packages

- `numpy>=2.4.1`
- `polars[rtcompat]>=1.37.1`
- `pyfastx>=2.3.0`
- `biopython>=1.86`

These are installed automatically by `pip install .` — no separate step is needed for them.

## Installation

### a) Create a conda environment

```bash
conda create -n gandi python=3.11
conda activate gandi
```

### b) Install the external dependencies through conda

```bash
conda install -c bioconda -c conda-forge skani prodigal prodigal-gv diamond
```

### c) Install the `gandi` package

The package isn't published on PyPI yet, so install it from a git clone:

```bash
git clone https://github.com/fbaltoumas/GanDi-classifier.git
cd GanDi-classifier
pip install .
```

For an editable install (picks up code changes without reinstalling), use `pip install -e .` instead. This also installs the Python packages listed in [Dependencies](#python-packages).

### Windows users

`gandi` and its dependencies (`skani`, `prodigal`, `prodigal-gv`, `diamond`) are not supported natively on Windows — `bioconda`, which provides these tools, only supports Linux and macOS. To run this tool on Windows:

1. Install and configure the Windows Subsystem for Linux (WSL). See the [official Microsoft WSL installation guide](https://learn.microsoft.com/en-us/windows/wsl/install).
2. Once inside your WSL Linux environment, proceed as normal — follow the steps above exactly as you would on native Linux.

### Docker

A `Dockerfile` is provided that bundles `gandi` together with all of its external dependencies (`skani`, `prodigal`, `prodigal-gv`, `diamond`) via bioconda, so no separate dependency installation is needed.

Build the image:

```bash
docker build -t gandi .
```

Run it (mount a local directory to `/data` so input/output/database files are accessible from the host). First download the reference database, then classify a genome against it:

```bash
docker run --rm -v "$(pwd)":/data gandi download-db -o /data/gandi_db

docker run --rm -v "$(pwd)":/data gandi classifier \
    -i /data/input.fasta -o /data/output_dir -w full \
    -c plasmids -d /data/gandi_db
```

Running the image with no arguments shows the top-level help text listing both subcommands (`docker run --rm gandi`).

Note: the `Dockerfile` currently installs the `gandi` package via `git clone` + `pip install .` from source, since the package isn't published on PyPI yet. Once it is, the corresponding `RUN` step can be swapped for a plain `pip install gandi`.

## Usage

### Overview

`gandi` is a single command with two subcommands:

- `gandi classifier` — classifies contigs/genomes against a GanDi reference database (ANI/AAI workflows).
- `gandi download-db` — downloads and builds the GanDi reference database.

```
$ gandi -h
usage: gandi [-h] <subcommand> ...

gandi: contig/genome classifier and reference-database downloader for the
Global Anaerobic Digestion (GanDi) database.

Classify a genome against the database:
  gandi classifier -i genome.fasta -o output_dir -c plasmids -d /path/to/database

Download and build the reference database:
  gandi download-db -o /path/to/database

positional arguments:
  <subcommand>
    classifier   Classify contigs/genomes against a GanDi reference database
                 (ANI/AAI workflows).
    download-db  Download and build the GanDi reference database.

options:
  -h, --help     show this help message and exit
```

### Database download

Download and build the reference database (only needs to be done once):

```bash
gandi download-db -o /path/to/database
```

This requires `skani` and `diamond` to already be on `PATH` (see [Installation](#installation)).

### Classification

Classify a genome or set of contigs against the database:

```bash
gandi classifier -i input.fasta -o output_dir -w full \
    -c plasmids -d /path/to/database
```

Run `gandi --help` for an overview of both subcommands, or `gandi classifier --help` / `gandi download-db --help` for their full lists of options.

### Examples

The `examples/` directory contains sample inputs for each category and workflow. Assuming a database has already been downloaded to `/path/to/database`:

MAG input, once per workflow type:

```bash
gandi classifier -i examples/mag_input.fna -o mag_full_output -c mags -d /path/to/database -w full
gandi classifier -i examples/mag_input.fna -o mag_ani_output -c mags -d /path/to/database -w ani
gandi classifier -i examples/mag_input.fna -o mag_aai_output -c mags -d /path/to/database -w aai
```

Isolate genome, against the MAGs database:

```bash
gandi classifier -i examples/isolate_genome_input.fna -o isolate_output -c mags -d /path/to/database
```

Single virus genome:

```bash
gandi classifier -i examples/single_virus.fna -o single_virus_output -c viruses -d /path/to/database
```

Multiple virus genomes in one FASTA file, using `--multiple_genomes` so each sequence is treated as a separate genome:

```bash
gandi classifier -i examples/multiple_viruses.fna -o multiple_viruses_output -c viruses -d /path/to/database --multiple_genomes
```

Single plasmid:

```bash
gandi classifier -i examples/single_plasmid.fna -o single_plasmid_output -c plasmids -d /path/to/database
```

## Troubleshooting

### `RuntimeError: NumPy was built with baseline optimizations... but your machine doesn't support...`

This means the CPU numpy is running on is missing instructions (SSSE3/SSE4.1/SSE4.2/POPCNT) that numpy's official PyPI wheels require. This is not a bug in `gandi` — it's most commonly seen on virtual machines where the hypervisor masks CPU features from the guest for live-migration compatibility (e.g. Hyper-V's "Processor Compatibility Mode"), even though the physical CPU fully supports them.

Fixes, in order of preference:

1. **If you manage the VM/hypervisor**: disable processor/CPU compatibility mode on the VM so its real CPU features are exposed to the guest. See [Configure processor compatibility mode in Hyper-V](https://learn.microsoft.com/en-us/windows-server/virtualization/hyper-v/configure-processor-compatibility-mode) (official Microsoft docs). This is the only fix that also restores full numpy performance.

2. **If you can't change the VM's CPU configuration**, rebuild numpy from source with a lower SIMD baseline. If you installed from a git clone, a helper script is included that detects whether this is actually needed and does it for you:
   ```bash
   python3 scripts/check_cpu_baseline.py
   ```
   Run it before `pip install .` on a fresh setup, or afterwards to repair an already-broken install (no need to reinstall the `gandi` package itself afterwards — just re-run `gandi`).

   If you don't have the repo cloned (e.g. installed via `pip install gandi`), run the equivalent command directly instead:
   ```bash
   pip install --force-reinstall numpy --no-binary numpy -Csetup-args=-Dcpu-baseline="none"
   ```
   This build is slower for numpy-heavy operations than the stock wheel, but functionally correct on any x86-64 CPU.

### `Illegal instruction (core dumped)` after importing polars

Same root cause as above (a CPU/VM missing instructions polars' default build requires — in this case a wider set, including AVX/AVX2/FMA/BMI), but a different failure mode: polars crashes the whole process outright (`SIGILL`) rather than raising a catchable Python exception, so this can't be turned into a friendly error message — it has to be fixed at the install level.

`gandi` already depends on `polars[rtcompat]` rather than plain `polars`, which installs an extra, broadly-compatible runtime (`polars-runtime-compat`) alongside the normal fast one; polars picks whichever one actually matches the CPU it's running on at import time, automatically and with no performance cost on capable CPUs. If you still hit this (e.g. an existing environment installed before this was added, or a version of `gandi` predating it), fix it directly:
```bash
pip install --force-reinstall "polars[rtcompat]"
```
