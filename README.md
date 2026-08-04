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

`gandi-classifier` also relies on the following external tools being installed separately and available on `PATH` (e.g. via conda/bioconda):
- [`skani`](https://github.com/bluenote-1215/skani) (ANI workflow)
- [`prodigal`](https://github.com/hyattpd/Prodigal) / [`prodigal-gv`](https://github.com/apcamargo/prodigal-gv) (AAI workflow)
- [`diamond`](https://github.com/bbuchfink/diamond) (AAI workflow)

## Usage

```bash
gandi-classifier -i input.fasta -o output_dir -w full \
    --skani_database /path/to/skani.db \
    --diamond_database_dmnd /path/to/diamond.dmnd
```

Run `gandi-classifier --help` for the full list of options.
