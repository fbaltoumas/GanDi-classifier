import os
import subprocess as sp
import argparse as ap
from pathlib import Path
from shutil import which, copy2

# 3rd-party
import pyfastx

def cmd_arguments():
    parser = ap.ArgumentParser()
    parser.add_argument("-o", "--output", required=False, default=None, help="Path to store databases (present working directory, i.e. '.' by default)")
    parser.add_argument("-t", "--threads", required=False, default=0, help='Number of CPU threads to use. Default: 0 (use all CPUs)')

    return parser.parse_args()


def get_exe_location(exe_name: str) -> str:
    exe_location = which(exe_name)
    if exe_location is None:
        raise RuntimeError(
                        f"Could not find the '{exe_name}' executable. "
                        f"Make sure it is installed and available on your PATH."
                    )
    return exe_location


def main():
    args = cmd_arguments()
    output_dir = Path("./")


    # step 0. First, check if the required dependencies are installed (skani, diamond). If they are not, report the error and exit
    get_exe_location("skani")
    get_exe_location("diamond")

    # step 1. Checking output path...
    if args.output is not None:
        output_dir = Path(args.output)

    if not output_dir.exists():
        print(f"Creating path {output_dir} to store databases...")
        output_dir.mkdir()
    else:
        if output_dir != Path("."):
            print("Path already exists, continuing...")

    # step 2. Downloading data
    

if __name__ == "__main__":
    main()
    