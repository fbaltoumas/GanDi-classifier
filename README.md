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

## Usage

```bash
gandi-classifier -i input.fasta -o output_dir -w full \
    --skani_database /path/to/skani.db \
    --diamond_database_dmnd /path/to/diamond.dmnd
```

Run `gandi-classifier --help` for the full list of options.
