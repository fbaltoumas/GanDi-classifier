import multiprocessing as mp
import subprocess as sp
from pathlib import Path
from shutil import which, copy2
import argparse as ap

# 3rd party modules
import pandas as pd
import pyfastx


def get_exe_location(exe_name: str) -> str:
    exe_location = which(exe_name)
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
    skani.add_argument("--multiple_genomes", action='store_true', required=False, help='If set, treat each sequence in the input FASTA as a separate genome')
    skani.add_argument("--skani_database", required=False, default=None, help="Path to an skani-compatible sketch database.")
    protein = parser.add_argument_group("AAI (prodigal + diamond) options")
    protein.add_argument("--metagenome", required=False, action='store_true', help="If set, run prodigal in 'meta' mode. Used when dealing with multiple genomes or metagenomes")
    protein.add_argument("--viral", required=False, action='store_true', help="If set, treat sequences as viruses")
    protein.add_argument("--diamond_database_dmnd", required=False, help='Path to the diamond protein database DMND file')

    annot = parser.add_argument_group("Post-processing and annotation options")
    annot.add_argument("--metadata", required=False, help="Path to the database metadata table (tab-delimited)")
    return parser.parse_args()


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
               "search",
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
        proc_output.columns = ['query', 'hit', 'ani', 'qcov', 'tcov']
        proc_output = proc_output.sort_values(by=['ani', 'qcov', 'tcov'], ascending=[False,False,False], ignore_index=True)
        self.ani = proc_output


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
            "6",
            # define columns in outfmt 6
            "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore",
            "qcovhsp", "scovhsp", #query and target coverage of the HSP
            "--query-cover",
            "50",
            "--subject-cover",
            "50"
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


class GenomeAAICalculator:
    def __init__(self,
                 diamond_output_file: Path) -> None:
        self.diamond_result = pd.read_csv(str(diamond_output_file), sep="\t", header=None,
                                          names=['query', 'hit', 'pid', 'aln_length', 'mismatches', 'gap_opens',
                                                 'q_start', 'q_end', 's_start', 's_end', 'evalue', 'bit_score',
                                                 'qcovhsp', 'scovhsp'])
        self.aai_result = None

    def calculate_aai(self, pid_cutoff: float = 30.0):
        diamond_result = self.diamond_result
        diamond_result['query_genome'] = diamond_result['query'].str.rsplit('_', n=1).str[0]
        diamond_result['hit_genome'] = diamond_result['hit'].str.rsplit('_', n=1).str[0]

        # sort table by bitscore and evalue first
        diamond_result_sorted = diamond_result.sort_values(by=['bit_score', 'evalue'], ascending=[False,True])
        # group by query and hit genome, and take the first row (best hit) for each query-hit genome pair
        best_hits_per_genome = diamond_result_sorted.groupby(['query', 'hit_genome'], as_index=False).first()
        # filter by percentage identity to get ortholog hits
        best_hits_filtered = best_hits_per_genome[best_hits_per_genome['pid'] >= pid_cutoff]

        aai_result = (
            best_hits_filtered.groupby(["query_genome", "hit_genome"], as_index=False).agg(
                aai = ("pid", "mean"),
                std = ("pid", "std"),
                ortholog_hits = ("pid", "count"),
                raw_bitscore = ("bit_score", "sum"),
                qcov = ("qcovhsp", "mean"),
                tcov = ("scovhsp", "mean")
            )
        ).rename(columns={"query_genome":'query', 'hit_genome':'hit'})
        aai_result['norm_bitscore'] = aai_result['raw_bitscore'] / aai_result['ortholog_hits']

        # reorder columns
        aai_result = aai_result[['query', 'hit', 'aai', 'std', 'ortholog_hits', 'qcov', 'tcov', 'raw_bitscore', 'norm_bitscore']]
        self.aai_result = aai_result


if __name__ == "__main__":
    args = cmd_arguments()

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

    # copy input to work path
    # also check that the input file is a valid FASTA file
    work_input_file = output_path / input_file.name
    copy2(input_file, work_input_file)
    try:
        fasta_input = pyfastx.Fasta(str(work_input_file))
    except Exception as e:
        print(f"Input file is not a valid FASTA file: {input_file}")
        exit()

    # check cpu count
    cpus = args.threads
    if args.threads <=0:
        cpus = mp.cpu_count()

    # set workflow
    workflow = args.workflow.lower()
    if workflow not in ['full', 'ani', 'aai']:
        workflow = 'full'


    # now, based on workflow type, run the appropriate function
    if workflow in ['ani', 'full']:
        if args.skani_database is None:
            print("No path to a sketch database given. Exiting...")
            exit()
        skani_db = args.skani_database

        skani_exe = get_exe_location("skani")

        multifasta = args.multiple_genomes
        output_prefix_raw = f"{output_path}/skani_search"
        skani_job = SkaniJob(
            work_input_file,
            output_prefix_raw,
            skani_db,
            multifasta,
            cpus,
            skani_exe
        )
        skani_job.run()
        skani_job.process_output()
        # write processed ANI result
        skani_job.ani.to_csv(f"{output_path}/ANI-result.tsv", sep="\t", index=False)


    if workflow in ['aai', 'full']:
        if args.diamond_database_dmnd is None:
            print("No path to a diamond database given. Exiting...")
            exit()
        
        # step 1: run prodigal to generate proteins
        prodigal_mode = 'single'
        if args.metagenome is True:
            prodigal_mode = 'meta'
        # but, also do a sequence length check
        if prodigal_mode =='single' and fasta_input.size < 20000:
            print("Input sequence length < 20000 bps, automatically switching to 'meta' mode for prodigal.")
            prodigal_mode = 'meta'
        output_prefix = f"{output_path}/prodigal"

        prodigal_exe = get_exe_location("prodigal")
        if args.viral is True:
            prodigal_exe = get_exe_location("prodigal-gv")
        prodigal_job = ProdigalJob(
            work_input_file,
            output_prefix,
            args.viral,
            prodigal_mode
        )
        prodigal_job.run()


        # step 2. Run diamond
        diamond_db = args.diamond_database_dmnd
        diamond_exe = get_exe_location("diamond")
        diamond_job = DiamondJob(
            f"{output_path}/prodigal.faa",
            f"{output_path}/diamond-search",
            diamond_db,
            cpus,
            diamond_exe
        )
        diamond_job.run()

        # step 3. calculate AAI
        aai_calculator = GenomeAAICalculator(
            f"{output_path}/diamond-search.blout"
        )
        aai_calculator.calculate_aai()
        aai_calculator.aai_result.to_csv(f"{output_path}/AAI-result.tsv", sep="\t", index=False)