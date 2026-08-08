import os
import subprocess as sp
import argparse as ap
import http.client
import logging
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from shutil import which, unpack_archive, rmtree
from functools import partial
from urllib.error import URLError
from urllib.request import urlopen, Request

logger = logging.getLogger(__name__)


URL = "https://www.gandi-db.org/datasets/gandi-classifier-db.tar.gz"


class DatabaseDownloader:
    def __init__(self, url: str, output_file: Path, # type: ignore
                 chunk_size: int = 32768, report_interval: float = 0.5, threads: int = 1,
                 max_retries: int = 5, retry_backoff: float = 1.0) -> None:
        self.url = url
        self.output_file = output_file
        self.chunk_size = chunk_size
        self.report_interval = report_interval
        self.threads = max(1, threads)
        self.max_retries = max_retries
        self.retry_backoff = retry_backoff
        self._downloaded = 0
        self._lock = threading.Lock()

    @staticmethod
    def _format_size(num_bytes):
        size = float(num_bytes)
        for unit in ["B", "KB", "MB", "GB"]:
            if size < 1024:
                return f"{size:.1f}{unit}"
            size /= 1024
        return f"{size:.1f}TB"

    @classmethod
    def _progress_line(cls, downloaded, total_size, elapsed):
        speed = downloaded / elapsed if elapsed > 0 else 0
        if total_size:
            pct = downloaded / total_size * 100
            return f"  {pct:5.1f}%  {cls._format_size(downloaded)} / {cls._format_size(total_size)}  ({cls._format_size(speed)}/s)"
        return f"  {cls._format_size(downloaded)} downloaded  ({cls._format_size(speed)}/s)"

    def _probe_range_support(self):
        request = Request(self.url, headers={"User-Agent": "Mozilla/5.0", "Range": "bytes=0-0"})
        with urlopen(request) as response:
            response.read(1)
            if response.status == 206:
                content_range = response.headers.get("Content-Range")
                if content_range and "/" in content_range:
                    total_size = int(content_range.rsplit("/", 1)[-1])
                    return total_size, True
            content_length = response.headers.get("Content-Length")
            total_size = int(content_length) if content_length is not None else None
            return total_size, False

    @staticmethod
    def _split_ranges(total_size, num_parts):
        part_size = total_size // num_parts
        ranges = []
        start = 0
        for i in range(num_parts):
            end = total_size - 1 if i == num_parts - 1 else start + part_size - 1
            ranges.append((start, end))
            start = end + 1
        return ranges

    def _download_range(self, start, end):
        # end is None means "download from start to EOF, no Range header" (server doesn't
        # support ranges, or this is the single-connection fallback). That case can't be
        # resumed on failure, since there's no byte offset to resume from -- it has to restart.
        resumable = end is not None
        position = start
        attempt = 0
        while True:
            headers = {"User-Agent": "Mozilla/5.0"}
            if resumable:
                headers["Range"] = f"bytes={position}-{end}"
            request = Request(self.url, headers=headers)
            try:
                with urlopen(request) as response, open(self.output_file, "r+b") as fl:
                    fl.seek(position)
                    for chunk in iter(partial(response.read, self.chunk_size), b""):
                        fl.write(chunk)
                        position += len(chunk)
                        with self._lock:
                            self._downloaded += len(chunk)
                return
            except (OSError, URLError, http.client.HTTPException) as e:
                if resumable and position > end:
                    # already received every byte of this range; the error came from a
                    # spurious extra read after the peer had nothing left to send.
                    return
                attempt += 1
                if attempt > self.max_retries:
                    raise RuntimeError(
                        f"Failed to download bytes {start}-{end if resumable else 'EOF'} "
                        f"after {self.max_retries} retries: {e}"
                    ) from e
                if not resumable:
                    # can't resume without range support: undo the counted progress and restart
                    with self._lock:
                        self._downloaded -= (position - start)
                    position = start
                time.sleep(self.retry_backoff * attempt)

    def _report_progress(self, total_size, start_time, stop_event):
        # The live, in-place-updating (\r) progress meter is deliberately plain print(),
        # not logging: it's a terminal UI element, not a discrete log event, and running
        # it through a logging.Formatter would break the in-place update (each call would
        # get its own timestamp/level prefix and newline).
        while not stop_event.is_set():
            stop_event.wait(self.report_interval)
            with self._lock:
                downloaded = self._downloaded
            print(f"\r{self._progress_line(downloaded, total_size, time.monotonic() - start_time)}", end="", flush=True)

    def download(self):
        total_size, range_support = self._probe_range_support()
        use_parallel = self.threads > 1 and total_size is not None and range_support
        if self.threads > 1 and not use_parallel:
            logger.warning("Server does not support parallel (range) downloads; falling back to a single connection.")

        logger.info(f"Downloading to {self.output_file}" + (f" ({self._format_size(total_size)})" if total_size else "") + "...")
        if use_parallel:
            logger.info(f"Using {self.threads} parallel connections.")

        self._downloaded = 0
        start_time = time.monotonic()
        stop_event = threading.Event()
        reporter = threading.Thread(target=self._report_progress, args=(total_size, start_time, stop_event), daemon=True)
        reporter.start()

        try:
            if use_parallel:
                with open(self.output_file, "wb") as fl:
                    fl.truncate(total_size)
                ranges = self._split_ranges(total_size, self.threads)
                errors = []
                with ThreadPoolExecutor(max_workers=self.threads) as executor:
                    futures = [executor.submit(self._download_range, start, end) for start, end in ranges]
                    for future in as_completed(futures):
                        try:
                            future.result()
                        except Exception as e:
                            errors.append(str(e))
                if errors:
                    raise RuntimeError("Parallel download failed:\n" + "\n".join(errors))
            else:
                open(self.output_file, "wb").close()
                self._download_range(0, None)
        finally:
            stop_event.set()
            reporter.join()

        print(f"\r{self._progress_line(self._downloaded, total_size, time.monotonic() - start_time)}")
        logger.info(f"Download complete: {self.output_file}")


def cmd_arguments():
    parser = ap.ArgumentParser()
    parser.add_argument("-o", "--output", required=False, default=None, help="Path to store databases (present working directory, i.e. '.' by default)")
    parser.add_argument("-t", "--threads", required=False, default=0, type=int, help='Number of CPU threads to use. Default: 0 (use all CPUs)')
    parser.add_argument("--quiet", required=False, action='store_true', help='Suppress informational logging output (only warnings and errors are shown).')

    return parser.parse_args()


def get_exe_location(exe_name: str) -> str:
    exe_location = which(exe_name)
    if exe_location is None:
        raise RuntimeError(
                        f"Could not find the '{exe_name}' executable. "
                        f"Make sure it is installed and available on your PATH."
                    )
    return exe_location



def skani_sketch(input_path: Path, output_path: Path, multifasta: bool = False, threads: int = 1) -> None:
    logger.info(f"Building skani sketch database: {input_path} -> {output_path} (threads={threads})")
    skani_exe = get_exe_location("skani")
    input_path = Path(input_path)
    output_path = Path(output_path)
    if output_path.exists():
        if output_path.is_dir():
            rmtree(output_path)
        else:
            output_path.unlink()
    cmd = [
        skani_exe,
        "sketch",
        "-t", str(threads),
        "-o", str(output_path)
    ]

    list_file = None
    if multifasta is True:
        cmd.append("-i")
        cmd.append(str(input_path))
    elif input_path.is_dir():
        # A directory of many individual genome files (e.g. mOTUs/MAGs). Passing each
        # file as its own argv entry can blow past the OS's ARG_MAX with large
        # directories, so list them in a temp file and use skani's "-l" instead.
        genome_files = sorted(p.resolve() for p in input_path.iterdir() if p.is_file())
        if not genome_files:
            raise RuntimeError(f"No genome files found in directory '{input_path}'.")
        fd, list_file = tempfile.mkstemp(prefix="skani_genome_list_", suffix=".txt")
        with os.fdopen(fd, "w") as f:
            f.write("\n".join(str(p) for p in genome_files) + "\n")
        cmd.extend(["-l", list_file])
    else:
        cmd.append(str(input_path))

    try:
        sp.run(cmd, check=True, capture_output=True, text=True)
    except FileNotFoundError as e:
        raise RuntimeError(
            f"Could not find the '{skani_exe}' executable. "
            f"Make sure it is installed and available on your PATH."
        ) from e
    except sp.CalledProcessError as e:
        raise RuntimeError(
            f"'{skani_exe}' failed on input '{input_path}' "
            f"Command: {' '.join([str(c) for c in cmd])}\n"
            f"Stderr: {e.stderr.strip() if e.stderr else '(no stderr output)'}"
        ) from e
    finally:
        if list_file is not None:
            Path(list_file).unlink(missing_ok=True)

def diamond_makedb(input_path: Path, output_path: Path, threads: int = 1) -> None:
    logger.info(f"Building diamond database: {input_path} -> {output_path} (threads={threads})")
    diamond_exe = get_exe_location("diamond")
    output_path = Path(output_path)
    if output_path.exists():
        if output_path.is_dir():
            rmtree(output_path)
        else:
            output_path.unlink()
    # diamond makedb --threads 4 --db plasmids.dmnd --in plasmids.faa
    cmd = [
        diamond_exe,
        "makedb",
        "--threads", str(threads),
        "--db", str(output_path),
        "--in", str(input_path)
    ]
    try:
        sp.run(cmd, check=True, capture_output=True, text=True)
    except FileNotFoundError as e:
        raise RuntimeError(
            f"Could not find the '{diamond_exe}' executable. "
            f"Make sure it is installed and available on your PATH."
        ) from e
    except sp.CalledProcessError as e:
        raise RuntimeError(
            f"'{diamond_exe}' failed on input '{input_path}' "
            f"Command: {' '.join([str(c) for c in cmd])}\n"
            f"Stderr: {e.stderr.strip() if e.stderr else '(no stderr output)'}"
        ) from e


def main():
    args = cmd_arguments()
    logging.basicConfig(
        level=logging.WARNING if args.quiet else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )

    output_dir = Path("./")


    # step 0. First, check if the required dependencies are installed (skani, diamond). If they are not, report the error and exit
    logger.info("=== Step 0: Checking required executables (skani, diamond) ===")
    get_exe_location("skani")
    get_exe_location("diamond")
    logger.info("Found skani and diamond on PATH.")

    # step 1. Check output path...
    logger.info("=== Step 1: Preparing output directory ===")
    if args.output is not None:
        output_dir = Path(args.output)

    if not output_dir.exists():
        logger.info(f"Creating path {output_dir} to store databases...")
        output_dir.mkdir()
    else:
        if not output_dir.is_dir():
            logger.error(f"Output path {output_dir} exists and is not a directory. Exiting...")
            exit()
        if output_dir != Path("."):
            logger.info("Path already exists, continuing...")

    # step 1b. Resolve thread count
    threads = args.threads
    if args.threads <= 0:
        if hasattr(os, "sched_getaffinity"):
            threads = len(os.sched_getaffinity(0))
        else:
            threads = os.cpu_count() or 1

    # step 2. Download and extract data
    logger.info("=== Step 2: Downloading and extracting the database archive ===")
    output_file = Path(f"{output_dir}/gandi-classifier-db.tar.gz")
    download_threads = threads if threads < 8 else 8
    logger.info(f"Using {threads} thread(s) for skani/diamond, capped at {download_threads} for the download itself.")
    downloader = DatabaseDownloader(URL, output_file, threads=download_threads)
    downloader.download()
    logger.info(f"Extracting {output_file}...")
    unpack_archive(output_file, output_dir, "gztar")
    logger.info("Extraction complete.")
    # step 3. delete the (large) downloaded archive now that it's extracted
    logger.info("=== Step 3: Removing the downloaded archive ===")
    logger.info(f"Deleting {output_file}...")
    output_file.unlink()

    skani_dir = Path(f"{output_dir}/database/skani")
    diamond_dir = Path(f"{output_dir}/database/diamond")

    # step 4. create skani databases
    logger.info("=== Step 4: Building skani sketch databases ===")
    # 4.1 plasmids
    logger.info("--- Step 4.1: Plasmids ---")
    # 4.1.1 unzip
    logger.info(f"Extracting {skani_dir / 'plasmids.tar.gz'}...")
    unpack_archive(skani_dir / "plasmids.tar.gz", skani_dir, "gztar")
    # 4.1.2 sketch
    skani_sketch(skani_dir / "plasmids.fna",
                 skani_dir / "plasmids", multifasta=True,
                 threads=threads)
    # 4.1.3 delete intermediate data
    logger.info("Cleaning up plasmid intermediate files...")
    (skani_dir / "plasmids.tar.gz").unlink()
    (skani_dir / "plasmids.fna").unlink()
    # 4.2 viruses
    logger.info("--- Step 4.2: Viruses ---")
    # 4.2.1 unzip
    logger.info(f"Extracting {skani_dir / 'viruses.tar.gz'}...")
    unpack_archive(skani_dir / "viruses.tar.gz", skani_dir, "gztar")
    # 4.2.2 sketch
    skani_sketch(skani_dir / "viruses.fna",
                 skani_dir / "viruses", multifasta=True,
                 threads=threads)
    # 4.2.3 delete intermediate data
    logger.info("Cleaning up virus intermediate files...")
    (skani_dir / "viruses.tar.gz").unlink()
    (skani_dir / "viruses.fna").unlink()
    # 4.3 mOTUs
    logger.info("--- Step 4.3: mOTUs/MAGs ---")
    # 4.3.1 unzip
    logger.info(f"Extracting {skani_dir / 'motu_genomes.tar.gz'}...")
    unpack_archive(skani_dir / "motu_genomes.tar.gz", skani_dir, "gztar")
    # 4.3.2 sketch
    skani_sketch(skani_dir / "motu_genomes",
                 skani_dir / "mags", multifasta=False,
                 threads=threads)
    # 4.3.3 delete intermediate data
    logger.info("Cleaning up mOTU/MAG intermediate files...")
    (skani_dir / "motu_genomes.tar.gz").unlink()
    rmtree(skani_dir / "motu_genomes")
    logger.info("skani databases complete.")

    # step 5. create diamond databases
    logger.info("=== Step 5: Building diamond databases ===")
    # 5.1 plasmids
    logger.info("--- Step 5.1: Plasmids ---")
    # 5.1.1 unzip
    logger.info(f"Extracting {diamond_dir / 'plasmids.tar.gz'}...")
    unpack_archive(diamond_dir / "plasmids.tar.gz", diamond_dir, "gztar")
    # 5.1.2 sketch
    diamond_makedb(diamond_dir / "plasmids.faa",
                 diamond_dir / "plasmids.dmnd",
                 threads=threads)
    # 5.1.3 delete intermediate data
    logger.info("Cleaning up plasmid intermediate files...")
    (diamond_dir / "plasmids.tar.gz").unlink()
    (diamond_dir / "plasmids.faa").unlink()
    # 5.2 viruses
    logger.info("--- Step 5.2: Viruses ---")
    # 5.2.1 unzip
    logger.info(f"Extracting {diamond_dir / 'viruses.tar.gz'}...")
    unpack_archive(diamond_dir / "viruses.tar.gz", diamond_dir, "gztar")
    # 5.2.2 sketch
    diamond_makedb(diamond_dir / "viruses.faa",
                 diamond_dir / "viruses.dmnd",
                 threads=threads)
    # 5.2.3 delete intermediate data
    logger.info("Cleaning up virus intermediate files...")
    (diamond_dir / "viruses.tar.gz").unlink()
    (diamond_dir / "viruses.faa").unlink()
    # 5.3 mOTUs
    logger.info("--- Step 5.3: mOTUs/MAGs ---")
    # 5.3.1 unzip
    logger.info(f"Extracting {diamond_dir / 'mags.tar.gz'}...")
    unpack_archive(diamond_dir / "mags.tar.gz", diamond_dir, "gztar")
    # 5.3.2 sketch
    diamond_makedb(diamond_dir / "mags.faa",
                 diamond_dir / "mags.dmnd",
                 threads=threads)
    # 5.3.3 delete intermediate data
    logger.info("Cleaning up mOTU/MAG intermediate files...")
    (diamond_dir / "mags.tar.gz").unlink()
    (diamond_dir / "mags.faa").unlink()
    logger.info("diamond databases complete.")
    logger.info(f"=== All done! Databases are ready in {output_dir} ===")

if __name__ == "__main__":
    main()
    