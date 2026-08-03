import multiprocessing as mp
import subprocess as sp
from pathlib import Path
import argparse as ap

# 3rd party modules
import pandas as pd
import pyfastx


def get_exe_location(exe_name: str) -> str:
    cmd = ["which", exe_name]
    proc = sp.Popen(cmd, stdout=sp.PIPE)
    out = proc.communicate()
    exe_location = out[0].decode('utf-8').rstrip()
    if exe_location == "" or exe_location is None:
        raise RuntimeError(
                        f"Could not find the '{exe_name}' executable. "
                        f"Make sure it is installed and available on your PATH."
                    )
    return exe_location


def cmd_arguments():
    parser = ap.ArgumentParser()
    general = parser.add_argument_group("Input/output and run type options")
    general.add_argument("-i", "--input", required=True, help='Input genome(s) file in FASTA format.')
    general.add_argument("-o", "--output", required=True, help="Output prefix")
    general.add_argument("-w", "--workflow", required=False, default="full",help="Workflow type. Can be one of 'ani' (genome-based), 'aai' (proteome-based), or 'full' (both). Default is 'full'")
    general.add_argument("-t", "--threads", required=False, default=1, type=int, help="Number of CPU threads to use. Default: 1. Use 0 to get all CPU threads")
    skani = parser.add_argument_group("ANI (skani) options")
    skani.add_argument("--multiple_genomes", action='store_true', required=False, help='If set, treat each sequencei in the input FASTA as a separate genome')
    skani.add_argument("--skani_database", required=False, default=None, help="Path to an skani-compatible sketch database.")
    skani.add_argument("--skani_exe", required=False, help="Path to the skani executable.")
    protein = parser.add_argument_group("AAI (prodigal + diamond) options")
    protein.add_argument("--metagenome", required=False, action='store_true', help="If set, run prodigal in 'meta' mode. Used when dealing with multiple genomes or metagenomes")
    protein.add_argument("--viral", required=False, action='store_true', help="If set, treat sequences as viruses")
    protein.add_argument("--diamond_database", required=False, help='Path to the diamond protein database')
    annot = parser.add_argument_group("Post-processing and annotation options")
    annot.add_argument("--metadata", required=False, help="Path to the database metadata table (tab-delimited)")
    return parser.parse_args()





class ProdigalJob:
    def __init__(self, input_file: Path,
                 output_prefix: Path, is_virus: bool = False,
                 mode: str = 'meta') -> None:
        self.input_file = input_file
        self.output_prefix = output_prefix
        self.is_virus = is_virus
        self.prodigal_exe = "prodigal"
        if self.is_virus is True:
            self.prodigal_exe = "prodigal-gv"
        self.mode = mode.lower()
        if self.mode not in ['single', 'meta']:
            self.mode = 'meta'
    
    def run(self):
        prodigal_output = f"{self.output_prefix}.gff"
        prodigal_genes = f"{self.output_prefix}.fna"
        prodigal_proteins = f"{self.output_prefix}.faa"
        cmd = [self.prodigal_exe,
               '-i', str(self.input_file),
               "-f", "gff",
               "-o", prodigal_output,
               "-a", prodigal_proteins,
               "-d", prodigal_genes,
               "-p", self.mode]
        try:
            sp.run(cmd, check=True, capture_output=True, text=True)
        except FileNotFoundError as e:
            raise RuntimeError(
                f"Could not find the '{self.prodigal_exe}' executable. "
                f"Make sure it is installed and available on your PATH."
            ) from e
        except sp.CalledProcessError as e:
            raise RuntimeError(
                f"'{self.prodigal_exe}' failed on input '{self.input_file}' "
                f"with exit code {e.returncode}.\n"
                f"Command: {' '.join(cmd)}\n"
                f"Stderr: {e.stderr.strip() if e.stderr else '(no stderr output)'}"
            ) from e


class SkaniJob:
    def __init__(self, input_file: Path,
                 output_prefix: Path,
                 database: Path,
                 multifasta_individual_sequences: bool = False,
                 cpus: int = 1, skani_exe: str = None) -> None:
        self.input_file = input_file
        self.output_prefix = output_prefix
        self.database = database
        self.qi = False
        if multifasta_individual_sequences is True:
            self.qi = True
        self.cpus = cpus
        # exe
        self.skani_exe = skani_exe
        # processed output:
        self.ani : pd.DataFrame = None
    
    def run(self):
        q_arg = "-q"
        output = f"{self.output_prefix}.skani"
        if self.qi is True:
            q_arg = "--pi"
        cmd = [self.skani_exe,
               q_arg,
               str(self.input_file),
               "-o",
               output,
               "-d",
               str(self.database),
               "-t",
               str(self.cpus)
        ]
        try:
            sp.run(cmd, check=True, capture_output=True, text=True)
        except FileNotFoundError as e:
            raise RuntimeError(
                f"Could not find the '{self.skani_exe}' executable. "
                f"Make sure it is installed and available on your PATH."
            ) from e
        except sp.CalledProcessError as e:
            raise RuntimeError(
                f"'{self.skani_exe}' failed on input '{self.input_file}' "
                f"against database '{self.database}' with exit code {e.returncode}.\n"
                f"Command: {' '.join(cmd)}\n"
                f"Stderr: {e.stderr.strip() if e.stderr else '(no stderr output)'}"
            ) from e

    def process_output(self):
        raw_output = pd.read_csv(f"{self.output_prefix}.skani", sep="\t")
        proc_output = raw_output[['Query_name', 'Ref_name', 'ANI', 'Align_fraction_query', 'Align_fraction_ref']]
        proc_output.columns ['query', 'hit', 'ani', 'qcov', 'tcov']
        self.ani = proc_output



def workflow_skani(input_file:Path,
                   output_path: Path,
                   skani_database: Path,
                   skani_exe: Path,
                   multifasta: bool = False, 
                   cpus: int = 1) -> bool:
    """"""
    skani_job = SkaniJob(
        input_file,
        output_path,
        skani_database,
        multifasta,
        cpus,
        skani_exe
    )
    skani_job.run()
    skani_job.process_output()
    return skani_job.ani



class DiamondJob:
    def __init__(self, input_file:Path,
                 output_prefix:Path,
                 database: Path,
                 cpus: int = 1, diamond_exe: str = None) -> None:
        self.input_file = input_file
        self.output_prefix  = output_prefix
        self.database = database
        self.cpus = cpus
        self.diamond_exe = diamond_exe

    def run(self):
        output = f"{self.output_prefix}.blout"

        cmd = [
            self.diamond_exe,
            "blastp",
            "--masking",
            "none",
            "-k",
            "1000",
            "-e",
            "1e-3",

            "--query",
            self.input_file,
            "--db",
            self.database,
            "--out",
            output,
            "--outfmt",
            "6"
        ]
        try:
            sp.run(cmd, check=True, capture_output=True, text=True)
        except FileNotFoundError as e:
            raise RuntimeError(
                f"Could not find the '{self.diamond_exe}' executable. "
                f"Make sure it is installed and available on your PATH."
            ) from e
        except sp.CalledProcessError as e:
            raise RuntimeError(
                f"'{self.diamond_exe}' failed on input '{self.input_file}' "
                f"against database '{self.database}' with exit code {e.returncode}.\n"
                f"Command: {' '.join(cmd)}\n"
                f"Stderr: {e.stderr.strip() if e.stderr else '(no stderr output)'}"
            ) from e





if __name__ == "__main__":
    args = cmd_arguments()
    print(type(args))

    # check if input file xists
    input_file = Path(args.input)
    if not input_file.exists():
        print(f"Incorrect input file name: {input_file}")
        exit()

    # check if output path exists, and if not, create it
    output_path = Path(args.output)
    if not output_path.exists():
        output_path.mkdir()
    else:
        print("Path already exists, continuing...")

    # check cpu count
    cpus = args.threads
    if args.threads <=0:
        cpus = mp.cpu_count()

    # set workflow
    workflow = args.workflow.lower()
    if workflow not in ['full', 'ani', 'aai']:
        workflow = 'ani'


    # now, based on workflow type, run the appropriate function
    if workflow in ['ani', 'full']:
        """"""
        if args.skani_database is None:
            print("No path to a sketch database given. Exiting...")
            exit()
        skani_db = args.skani_database
        skani_exe = args.skani_exe
        if skani_exe is None:
            skani_exe = get_exe_location("skani")

        multifasta = args.multiple_genomes

        skani_job = SkaniJob(
            input_file,
            output_path,
            skani_db,
            multifasta,
            cpus,
            skani_exe
        )
        skani_job.run()
        skani_job.process_output()


    elif workflow in ['aai', 'full']:
        """"""
    else:
        print("Non-recognized workflow type. Exiting...")
        exit(1)