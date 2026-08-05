# GanDi-classifier
Contig classifier for the Global Anaerobic Digestion (GanDi) database

## Installation

```bash
pip install gandi
```

This installs the `gandi-classifier` command on your `PATH`.

### Installing from source

```bash
git clone https://github.com/fbaltoumas/GanDi-classifier.git
cd GanDi-classifier
pip install .
```

For an editable install (picks up code changes without reinstalling), use `pip install -e .` instead.

`gandi-classifier` also relies on the external tools listed below being installed separately and available on `PATH` — see [Installing dependencies](#installing-dependencies).

## Installing dependencies

`gandi-classifier` shells out to `skani`, `prodigal`/`prodigal-gv`, and `diamond`. These are not Python packages and are not installed by `pip install gandi` — they need to be installed separately and be available on `PATH`.

The easiest way is via [bioconda](https://bioconda.github.io/), which covers all four tools:

```bash
conda install -c bioconda -c conda-forge skani prodigal prodigal-gv diamond
```

Or install each one individually:

- **[skani](https://github.com/bluenote-1215/skani)** (ANI workflow)
  ```bash
  conda install -c bioconda skani
  ```
  Alternatively, download a prebuilt binary from the [skani releases page](https://github.com/bluenote-1215/skani/releases), or install via `cargo install skani` if you have a Rust toolchain.

- **[prodigal](https://github.com/hyattpd/Prodigal)** (AAI workflow)
  ```bash
  conda install -c bioconda prodigal
  ```
  On Debian/Ubuntu, it's also available via `apt`:
  ```bash
  sudo apt install prodigal
  ```

- **[prodigal-gv](https://github.com/apcamargo/prodigal-gv)** (AAI workflow, used with `--viral`)
  ```bash
  conda install -c bioconda prodigal-gv
  ```

- **[diamond](https://github.com/bbuchfink/diamond)** (AAI workflow)
  ```bash
  conda install -c bioconda diamond
  ```
  Alternatively, download a prebuilt binary from the [diamond releases page](https://github.com/bbuchfink/diamond/releases).

## Windows users

`gandi-classifier` and its dependencies (`skani`, `prodigal`, `prodigal-gv`, `diamond`) are not supported natively on Windows — `bioconda`, which provides these tools, only supports Linux and macOS. To run this tool on Windows:

1. Install and configure the Windows Subsystem for Linux (WSL). See the [official Microsoft WSL installation guide](https://learn.microsoft.com/en-us/windows/wsl/install).
2. Once inside your WSL Linux environment, proceed as normal — follow the [Installation](#installation) and [Installing dependencies](#installing-dependencies) sections above exactly as you would on native Linux.

## Docker

A `Dockerfile` is provided that bundles `gandi-classifier` together with all of its external dependencies (`skani`, `prodigal`, `prodigal-gv`, `diamond`) via bioconda, so no separate dependency installation is needed.

Build the image:

```bash
docker build -t gandi-classifier .
```

Run it (mount a local directory to `/data` so input/output files are accessible from the host):

```bash
docker run --rm -v "$(pwd)":/data gandi-classifier \
    -i input.fasta -o output_dir -w full \
    --skani_database /data/skani.db \
    --diamond_database_dmnd /data/diamond.dmnd
```

Running the image with no arguments shows the help text (`docker run --rm gandi-classifier`).

Note: the `Dockerfile` currently installs `gandi-classifier` via `git clone` + `pip install .` from source, since the package isn't published on PyPI yet. Once it is, the corresponding `RUN` step can be swapped for a plain `pip install gandi`.

## Usage

```bash
gandi-classifier -i input.fasta -o output_dir -w full \
    --skani_database /path/to/skani.db \
    --diamond_database_dmnd /path/to/diamond.dmnd
```

Run `gandi-classifier --help` for the full list of options.
