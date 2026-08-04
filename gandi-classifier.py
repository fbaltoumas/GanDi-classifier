import os
import subprocess as sp
from pathlib import Path
from shutil import which, copy2
import argparse as ap

# 3rd party modules
import numpy as np
import pandas as pd
import pyfastx
import Bio.Phylo as Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor


def get_exe_location(exe_name: str) -> str:
    exe_location = which(exe_name)
    if exe_location is None:
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

    phylogeny = parser.add_argument_group("Phylogeny tree options")
    phylogeny.add_argument("--create_tree", required=False, action='store_true', help="If set, construct a phylogenetic tree from the ANI/AAI distances and save it in Newick format.")
    phylogeny.add_argument("--tree_method", required=False, default='nj', choices=['nj', 'upgma'], help="Method used to construct the phylogenetic tree. Can be 'nj' (neighbor-joining) or 'upgma'. Default: 'nj'")
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
        proc_output['qcov'] = proc_output['qcov'] * 100
        proc_output['tcov'] = proc_output['tcov'] * 100
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
            str(self.input_file),
            "--db",
            str(self.database),
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
        aai_result = aai_result.sort_values(by=['aai', 'qcov', 'tcov'], ascending=[False,False,False], ignore_index=True)
        self.aai_result = aai_result


class Phylogeny:
    def __init__(self, result: pd.DataFrame, metric: str, method: str = 'nj') -> None:
        self.result = result
        self.value_column = metric
        self.method = method.lower()
        if self.method not in ['nj', 'upgma']:
            self.method = 'nj'
        self.tree = None

    def _build_distance_matrix(self) -> DistanceMatrix:
        df = self.result[['query', 'hit', self.value_column]].dropna(subset=[self.value_column])
        genomes = sorted(set(df['query']) | set(df['hit']))
        if len(genomes) < 2:
            raise RuntimeError(
                f"Not enough genomes with a valid {self.value_column.upper()} value "
                f"to construct a phylogenetic tree."
            )
        genome_index = {genome: i for i, genome in enumerate(genomes)}
        n = len(genomes)
        distances = np.full((n, n), np.nan)
        np.fill_diagonal(distances, 0.0)
        for _, row in df.iterrows():
            i, j = genome_index[row['query']], genome_index[row['hit']]
            distance = (100 - row[self.value_column]) / 100
            distances[i, j] = distance
            distances[j, i] = distance

        missing = np.isnan(distances)
        n_missing_pairs = int(missing.sum() / 2)
        if n_missing_pairs > 0:
            observed = distances[~missing]
            fill_value = float(np.max(observed)) if observed.size > 0 else 1.0
            print(
                f"Phylogeny: {n_missing_pairs} genome pairs have no {self.value_column.upper()} value; "
                f"filling with the maximum observed distance ({fill_value:.4f})."
            )
            distances[missing] = fill_value

        matrix = [list(distances[i, :i + 1]) for i in range(n)]
        return DistanceMatrix(names=genomes, matrix=matrix)

    def build_tree(self):
        distance_matrix = self._build_distance_matrix()
        constructor = DistanceTreeConstructor()
        if self.method == 'upgma':
            self.tree = constructor.upgma(distance_matrix)
        else:
            self.tree = constructor.nj(distance_matrix)
        return self.tree

    def write_tree(self, output_file: Path) -> None:
        if self.tree is None:
            raise RuntimeError("Tree has not been constructed yet. Call build_tree() first.")
        Phylo.write(self.tree, str(output_file), "newick")


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
    except Exception:
        print(f"Input file is not a valid FASTA file: {input_file}")
        exit()

    # check cpu count
    cpus = args.threads
    if args.threads <=0:
        if hasattr(os, "sched_getaffinity"):
            cpus = len(os.sched_getaffinity(0))
        else:
            cpus = os.cpu_count() or 1

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

    # now, depending on the workflow, choose what result will be the final output
    if workflow == 'ani':
        result = skani_job.ani[['query', 'hit', 'ani', 'qcov', 'tcov']].rename(
            columns={'qcov': 'genome_qcov', 'tcov': 'genome_tcov'}
        )
    elif workflow == 'aai':
        result = aai_calculator.aai_result[['query', 'hit', 'aai', 'qcov', 'tcov']].rename(
            columns={'qcov': 'proteome_qcov', 'tcov': 'proteome_tcov'}
        )
    else:
        # for full workflow, we need to merge the ANI and AAI results.
        # we will keep only the ani/aai, qcov and tcov columns for each and rename them accordingly
        ani_renamed = skani_job.ani[['query', 'hit', 'ani', 'qcov', 'tcov']].rename(
            columns={'qcov': 'genome_qcov', 'tcov': 'genome_tcov'}
        )
        aai_renamed = aai_calculator.aai_result[['query', 'hit', 'aai', 'qcov', 'tcov']].rename(
            columns={'qcov': 'proteome_qcov', 'tcov': 'proteome_tcov'}
        )
        result = ani_renamed.merge(
            aai_renamed,
            on=['query', 'hit'],
            how='outer',
            indicator=True
        )
        match_counts = result['_merge'].value_counts()
        print(
            f"ANI/AAI merge: {match_counts.get('both', 0)} query-hit pairs matched in both, "
            f"{match_counts.get('left_only', 0)} ANI-only, {match_counts.get('right_only', 0)} AAI-only."
        )
        result = result.drop(columns='_merge')
        result = result.sort_values(by=['ani', 'genome_qcov', 'genome_tcov', 'aai', 'proteome_qcov', 'proteome_tcov'], ascending=[False,False,False,False,False,False], ignore_index=True)

    # if metadata is provided, we will use it to annotate the results based on the 'hit' column
    if args.metadata is not None:
        metadata_file = Path(args.metadata)
        if not metadata_file.exists():
            print(f"Metadata file does not exist: {metadata_file}. Exiting...")
            exit()
        metadata = pd.read_csv(metadata_file, sep="\t")
        if 'hit' not in metadata.columns:
            print(f"Metadata file does not contain a 'hit' column. Exiting...")
            exit()
        result = result.merge(metadata, on='hit', how='left', indicator=True)
        metadata_match_counts = result['_merge'].value_counts()
        print(
            f"Metadata merge: {metadata_match_counts.get('both', 0)} rows matched metadata, "
            f"{metadata_match_counts.get('left_only', 0)} rows had no metadata match."
        )
        result = result.drop(columns='_merge')

    # save the processed result to a file
    result.to_csv(f"{output_path}/final-result.tsv", sep="\t", index=False)
    print(f"Final result saved to: {output_path}/final-result.tsv")

    # if requested, construct phylogenetic tree(s) from the ANI/AAI distances
    if args.create_tree:
        metrics_to_build = []
        if workflow in ['ani', 'full']:
            metrics_to_build.append('ani')
        if workflow in ['aai', 'full']:
            metrics_to_build.append('aai')

        for metric in metrics_to_build:
            tree_file = output_path / f"phylogeny_{metric}.tree"
            try:
                phylogeny = Phylogeny(result, metric, args.tree_method)
                phylogeny.build_tree()
            except RuntimeError as e:
                print(f"Could not construct {metric.upper()}-based phylogenetic tree: {e}. Skipping...")
                continue
            phylogeny.write_tree(tree_file)
            print(f"{metric.upper()}-based phylogenetic tree saved to: {tree_file}")