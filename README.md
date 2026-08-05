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
    --diamond_database /data/diamond.dmnd
```

Running the image with no arguments shows the help text (`docker run --rm gandi-classifier`).

Note: the `Dockerfile` currently installs `gandi-classifier` via `git clone` + `pip install .` from source, since the package isn't published on PyPI yet. Once it is, the corresponding `RUN` step can be swapped for a plain `pip install gandi`.

## Usage

```bash
gandi-classifier -i input.fasta -o output_dir -w full \
    --skani_database /path/to/skani.db \
    --diamond_database /path/to/diamond.dmnd
```

Run `gandi-classifier --help` for the full list of options.

## Troubleshooting

### `RuntimeError: NumPy was built with baseline optimizations... but your machine doesn't support...`

This means the CPU numpy is running on is missing instructions (SSSE3/SSE4.1/SSE4.2/POPCNT) that numpy's official PyPI wheels require. This is not a bug in `gandi-classifier` — it's most commonly seen on virtual machines where the hypervisor masks CPU features from the guest for live-migration compatibility (e.g. Hyper-V's "Processor Compatibility Mode"), even though the physical CPU fully supports them.

Fixes, in order of preference:

1. **If you manage the VM/hypervisor**: disable processor/CPU compatibility mode on the VM so its real CPU features are exposed to the guest. See [Configure processor compatibility mode in Hyper-V](https://learn.microsoft.com/en-us/windows-server/virtualization/hyper-v/configure-processor-compatibility-mode) (official Microsoft docs). This is the only fix that also restores full numpy performance.

2. **If you can't change the VM's CPU configuration**, rebuild numpy from source with a lower SIMD baseline. If you installed from a git clone, a helper script is included that detects whether this is actually needed and does it for you:
   ```bash
   python3 scripts/check_cpu_baseline.py
   ```
   Run it before `pip install .` on a fresh setup, or afterwards to repair an already-broken install (no need to reinstall the `gandi` package itself afterwards — just re-run `gandi-classifier`).

   If you don't have the repo cloned (e.g. installed via `pip install gandi`), run the equivalent command directly instead:
   ```bash
   pip install --force-reinstall numpy --no-binary numpy -Csetup-args=-Dcpu-baseline="none"
   ```
   This build is slower for numpy-heavy operations than the stock wheel, but functionally correct on any x86-64 CPU.

### `Illegal instruction (core dumped)` after importing polars

Same root cause as above (a CPU/VM missing instructions polars' default build requires — in this case a wider set, including AVX/AVX2/FMA/BMI), but a different failure mode: polars crashes the whole process outright (`SIGILL`) rather than raising a catchable Python exception, so this can't be turned into a friendly error message — it has to be fixed at the install level.

`gandi-classifier` already depends on `polars[rtcompat]` rather than plain `polars`, which installs an extra, broadly-compatible runtime (`polars-runtime-compat`) alongside the normal fast one; polars picks whichever one actually matches the CPU it's running on at import time, automatically and with no performance cost on capable CPUs. If you still hit this (e.g. an existing environment installed before this was added, or a version of `gandi` predating it), fix it directly:
```bash
pip install --force-reinstall "polars[rtcompat]"
```
