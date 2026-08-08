import os
import subprocess as sp
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from shutil import which, copy2
import argparse as ap

# 3rd party modules
try:
    import numpy as np
except RuntimeError as e:
    if "baseline optimizations" in str(e):
        raise RuntimeError(
            "numpy failed to load because this CPU is missing instructions "
            "(SSE3/SSSE3/SSE4.1/SSE4.2/POPCNT) that the installed numpy build requires.\n"
            "This commonly happens on virtual machines where the hypervisor masks CPU "
            "features for live-migration compatibility (e.g. Hyper-V 'Processor "
            "Compatibility Mode') even though the physical CPU supports them.\n\n"
            "To fix this, either:\n"
            "  1. Disable processor/CPU compatibility mode on the VM so its real CPU "
            "features are exposed to the guest, or\n"
            "  2. Rebuild numpy from source with a lower SIMD baseline:\n"
            "       pip uninstall numpy\n"
            "       pip install numpy --no-binary numpy -Csetup-args=-Dcpu-baseline=\"none\"\n"
        ) from e
    raise
import polars as pl
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
    general.add_argument("-i", "--input", required=True, type=str, help="Input genome(s) file in FASTA format.")
    general.add_argument("-o", "--output", required=True, type=str, help="Output prefix")
    general.add_argument("-c", "--category", required=True, type=str, help="Data category to search. Must be one of these three: viruses, plasmids, or mags")
    general.add_argument("-d", "--database", required=True, type=str, help="Path to the GanDi-classifier database, e.g. '~/gandi/database/'")
    general.add_argument("-w", "--workflow", required=False, type=str, default="full", help="Workflow type. Can be one of 'ani' (genome-based), 'aai' (proteome-based), or 'full' (both). Default is 'full'")
    general.add_argument("-m", "--multiple_genomes", action='store_true', required=False, help='If set, treat each sequence in the input FASTA as a separate genome.  Enables "-qi" option in skani and "-p meta" in prodigal')
    general.add_argument("-t", "--threads", required=False, default=1, type=int, help="Number of CPU threads to use. Default: 1. Use 0 to get all CPU threads")

    phylogeny = parser.add_argument_group("Phylogeny tree options")
    phylogeny.add_argument("--create_tree", required=False, action='store_true', help="If set, construct a phylogenetic tree from the ANI/AAI distances and save it in Newick format.")
    phylogeny.add_argument("--tree_method", required=False, default='nj', choices=['nj', 'upgma'], help="Method used to construct the phylogenetic tree. Can be 'nj' (neighbor-joining) or 'upgma'. Default: 'nj'")
    return parser.parse_args()



class Database:
    def __init__(self, database_type: str, database_path: Path) -> None:
        if database_type.lower() not in ['viruses', 'plasmids', 'mags']:
            raise ValueError(
                f"Search database must be one of 'plasmids', 'viruses', or 'mags' (got '{database_type}')"
            )
        database_path = Path(database_path)
        if not database_path.exists():
            raise FileNotFoundError(
                f"{str(database_path)} does not exist"
            )
        self.database_type = database_type.lower()
        self.database_root = database_path
        self.skani_db  = Path(f"{self.database_root}/skani/{self.database_type}")
        self.diamond_db = Path(f"{self.database_root}/diamond/{self.database_type}.dmnd")
        self.metadata_db = Path(f"{self.database_root}/metadata/{self.database_type}.tsv.gz")
        for db in [self.skani_db, self.diamond_db, self.metadata_db]:
            if not db.exists():
                raise FileNotFoundError(
                    f"Database file {db} not found in the supplied path."
                )






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
        self.ani : pl.DataFrame = None

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
        raw_output = pl.read_csv(f"{self.output_prefix}.skani", separator="\t", infer_schema_length=None)
        proc_output = raw_output.select(['Query_name', 'Ref_name', 'ANI', 'Align_fraction_query', 'Align_fraction_ref'])
        proc_output = proc_output.rename({
            'Query_name': 'query', 'Ref_name': 'seq_name', 'ANI': 'ani',
            'Align_fraction_query': 'qcov', 'Align_fraction_ref': 'tcov'
        })
        proc_output = proc_output.with_columns([
            (pl.col('qcov') * 100).alias('qcov'),
            (pl.col('tcov') * 100).alias('tcov'),
        ])
        proc_output = proc_output.sort(['ani', 'qcov', 'tcov'], descending=[True, True, True], nulls_last=True)
        self.ani = proc_output


class ProdigalJob:
    def __init__(self, input_file: Path,
                 output_prefix: Path,
                 mode: str = 'meta', cpus: int = 1, prodigal_exe: str = None) -> None:
        self.input_file = input_file
        self.output_prefix = output_prefix
        self.mode = mode.lower()
        if self.mode not in ['single', 'meta']:
            self.mode = 'meta'
        self.cpus = cpus
        self.prodigal_exe = prodigal_exe

    def _run_prodigal(self, input_file: Path, output_prefix: Path) -> None:
        prodigal_output = f"{output_prefix}.gff"
        prodigal_genes = f"{output_prefix}.fna"
        prodigal_proteins = f"{output_prefix}.faa"
        cmd = [self.prodigal_exe,
               '-i', str(input_file),
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
                f"'{self.prodigal_exe}' failed on input '{input_file}' "
                f"with exit code {e.returncode}.\n"
                f"Command: {' '.join(cmd)}\n"
                f"Stderr: {e.stderr.strip() if e.stderr else '(no stderr output)'}"
            ) from e

    def run(self):
        self._run_prodigal(self.input_file, self.output_prefix)

    def run_parallel(self):
        records = list(pyfastx.Fasta(str(self.input_file)))
        num_chunks = max(1, min(self.cpus, len(records)))

        if num_chunks == 1:
            self.run()
            return

        chunks = [records[i::num_chunks] for i in range(num_chunks)]
        chunk_prefixes = [f"{self.output_prefix}.chunk{i}" for i in range(num_chunks)]
        for chunk_prefix, chunk_records in zip(chunk_prefixes, chunks):
            with open(f"{chunk_prefix}.fasta", 'w') as handle:
                for record in chunk_records:
                    handle.write(f">{record.name}\n{record.seq}\n")

        try:
            with ThreadPoolExecutor(max_workers=num_chunks) as executor:
                futures = [
                    executor.submit(self._run_prodigal, f"{chunk_prefix}.fasta", chunk_prefix)
                    for chunk_prefix in chunk_prefixes
                ]
                errors = []
                for future in as_completed(futures):
                    try:
                        future.result()
                    except RuntimeError as e:
                        errors.append(str(e))
            if errors:
                raise RuntimeError("Parallel prodigal run failed:\n" + "\n".join(errors))

            # merge chunk outputs into the final gff/fna/faa files
            for extension in ['gff', 'fna', 'faa']:
                with open(f"{self.output_prefix}.{extension}", 'w') as out_handle:
                    for chunk_prefix in chunk_prefixes:
                        with open(f"{chunk_prefix}.{extension}", 'r') as in_handle:
                            out_handle.write(in_handle.read())
        finally:
            # delete the intermediate chunk files, whether the run succeeded or failed
            for chunk_prefix in chunk_prefixes:
                for extension in ['fasta', 'gff', 'fna', 'faa']:
                    Path(f"{chunk_prefix}.{extension}").unlink(missing_ok=True)


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
        self.diamond_result = pl.read_csv(
            str(diamond_output_file), separator="\t", has_header=False, infer_schema_length=None,
            new_columns=['query', 'hit', 'pid', 'aln_length', 'mismatches', 'gap_opens',
                         'q_start', 'q_end', 's_start', 's_end', 'evalue', 'bit_score',
                         'qcovhsp', 'scovhsp']
        )
        self.aai_result = None

    def calculate_aai(self, pid_cutoff: float = 30.0):
        diamond_result = self.diamond_result.with_columns([
            pl.col('query').str.replace(r'_[^_]*$', '').alias('query_genome'),
            pl.col('hit').str.replace(r'_[^_]*$', '').alias('hit_genome'),
        ])

        # sort table by bitscore and evalue first
        diamond_result_sorted = diamond_result.sort(['bit_score', 'evalue'], descending=[True, False], nulls_last=True)
        # group by query and hit genome, and take the first row (best hit) for each query-hit genome pair
        best_hits_per_genome = diamond_result_sorted.group_by(['query', 'hit_genome'], maintain_order=True).first()
        # filter by percentage identity to get ortholog hits
        best_hits_filtered = best_hits_per_genome.filter(pl.col('pid') >= pid_cutoff)

        aai_result = best_hits_filtered.group_by(['query_genome', 'hit_genome']).agg([
            pl.col('pid').mean().alias('aai'),
            pl.col('pid').std().alias('std'),
            pl.col('pid').count().alias('ortholog_hits'),
            pl.col('bit_score').sum().alias('raw_bitscore'),
            pl.col('qcovhsp').mean().alias('qcov'),
            pl.col('scovhsp').mean().alias('tcov'),
        ])
        aai_result = aai_result.rename({'query_genome': 'query', 'hit_genome': 'seq_name'})
        aai_result = aai_result.with_columns(
            (pl.col('raw_bitscore') / pl.col('ortholog_hits')).alias('norm_bitscore')
        )

        # reorder columns
        aai_result = aai_result.select(['query', 'seq_name', 'aai', 'std', 'ortholog_hits', 'qcov', 'tcov', 'raw_bitscore', 'norm_bitscore'])
        aai_result = aai_result.sort(['aai', 'qcov', 'tcov'], descending=[True, True, True], nulls_last=True)
        self.aai_result = aai_result


class Phylogeny:
    def __init__(self, result: pl.DataFrame, metric: str, method: str = 'nj') -> None:
        self.result = result
        self.value_column = metric
        self.method = method.lower()
        if self.method not in ['nj', 'upgma']:
            self.method = 'nj'
        self.tree = None

    def _build_distance_matrix(self) -> DistanceMatrix:
        df = self.result.select(['query', 'seq_name', self.value_column]).drop_nulls(subset=[self.value_column])
        genomes = sorted(set(df['query'].to_list()) | set(df['seq_name'].to_list()))
        if len(genomes) < 2:
            raise RuntimeError(
                f"Not enough genomes with a valid {self.value_column.upper()} value "
                f"to construct a phylogenetic tree."
            )
        genome_index = {genome: i for i, genome in enumerate(genomes)}
        n = len(genomes)
        distances = np.full((n, n), np.nan)
        np.fill_diagonal(distances, 0.0)
        for row in df.iter_rows(named=True):
            i, j = genome_index[row['query']], genome_index[row['seq_name']]
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


def main():
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


    # Define database
    database = Database(args.category, args.database)


    # set workflow
    workflow = args.workflow.lower()
    if workflow not in ['full', 'ani', 'aai']:
        workflow = 'full'
   
    # now, based on workflow type, run the appropriate function
    if workflow in ['ani', 'full']:
        skani_db = database.skani_db

        skani_exe = get_exe_location("skani")
        skani_multifasta = args.multiple_genomes
        output_prefix_raw = f"{output_path}/skani_search"
        skani_job = SkaniJob(
            work_input_file,
            output_prefix_raw,
            skani_db,
            skani_multifasta,
            cpus,
            skani_exe
        )
        skani_job.run()
        skani_job.process_output()
        # write processed ANI result
        skani_job.ani.write_csv(f"{output_path}/ANI-result.tsv", separator="\t")


    if workflow in ['aai', 'full']:
        # step 1: run prodigal to generate proteins
        # first, check if multiple_genomes is true and if so, set to meta
        prodigal_mode = 'single'
        if args.multiple_genomes:
            prodigal_mode = 'meta'
        # also do a sequence length check. if below 20000, automatically set to 'meta'
        if prodigal_mode == 'single' and fasta_input.size < 20000:
            print("Input sequence length < 20000 bps, automatically switching to 'meta' mode for prodigal.")
            prodigal_mode = 'meta'
        output_prefix = f"{output_path}/prodigal"

        prodigal_exe = get_exe_location("prodigal")
        if args.category.lower() == 'viruses':
            prodigal_exe = get_exe_location("prodigal-gv")
        prodigal_job = ProdigalJob(
            work_input_file,
            output_prefix,
            prodigal_mode,
            cpus,
            prodigal_exe
        )
        if cpus >= 2:
            prodigal_job.run_parallel()
        else:
            prodigal_job.run()


        # step 2. Run diamond
        diamond_db = database.diamond_db
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
        aai_calculator.aai_result.write_csv(f"{output_path}/AAI-result.tsv", separator="\t")

    # now, depending on the workflow, choose what result will be the final output
    if workflow == 'ani':
        result = skani_job.ani.select(['query', 'seq_name', 'ani', 'qcov', 'tcov']).rename(
            {'qcov': 'genome_qcov', 'tcov': 'genome_tcov'}
        )
    elif workflow == 'aai':
        result = aai_calculator.aai_result.select(['query', 'seq_name', 'aai', 'qcov', 'tcov']).rename(
            {'qcov': 'proteome_qcov', 'tcov': 'proteome_tcov'}
        )
    else:
        # for full workflow, we need to merge the ANI and AAI results.
        # we will keep only the ani/aai, qcov and tcov columns for each and rename them accordingly
        ani_renamed = skani_job.ani.select(['query', 'seq_name', 'ani', 'qcov', 'tcov']).rename(
            {'qcov': 'genome_qcov', 'tcov': 'genome_tcov'}
        ).with_columns(pl.lit(True).alias('_ani_side'))
        aai_renamed = aai_calculator.aai_result.select(['query', 'seq_name', 'aai', 'qcov', 'tcov']).rename(
            {'qcov': 'proteome_qcov', 'tcov': 'proteome_tcov'}
        ).with_columns(pl.lit(True).alias('_aai_side'))
        result = ani_renamed.join(
            aai_renamed,
            on=['query', 'seq_name'],
            how='full',
            coalesce=True
        )
        match_counts = result.select([
            (pl.col('_ani_side').is_not_null() & pl.col('_aai_side').is_not_null()).sum().alias('both'),
            (pl.col('_ani_side').is_not_null() & pl.col('_aai_side').is_null()).sum().alias('ani_only'),
            (pl.col('_ani_side').is_null() & pl.col('_aai_side').is_not_null()).sum().alias('aai_only'),
        ]).row(0, named=True)
        both, ani_only, aai_only = match_counts['both'], match_counts['ani_only'], match_counts['aai_only']
        print(
            f"ANI/AAI merge: {both} query-hit pairs matched in both, "
            f"{ani_only} ANI-only, {aai_only} AAI-only."
        )
        result = result.drop(['_ani_side', '_aai_side'])
        result = result.sort(
            ['ani', 'genome_qcov', 'genome_tcov', 'aai', 'proteome_qcov', 'proteome_tcov'],
            descending=[True, True, True, True, True, True],
            nulls_last=True
        )

    # Annotate results based on the 'seq_name' column, by merging with the metadata df

    metadata = pl.read_csv(database.metadata_db, separator="\t", infer_schema_length=None)
    if 'seq_name' not in metadata.columns:
        print(f"Metadata file does not contain a 'seq_name' column. Exiting...")
        exit()
    metadata = metadata.with_columns(pl.lit(True).alias('_metadata_matched'))
    result = result.join(metadata, on='seq_name', how='left', coalesce=True)
    matched = int(result.select(pl.col('_metadata_matched').is_not_null().sum()).item())
    unmatched = result.height - matched
    print(
        f"Metadata merge: {matched} rows matched metadata, "
        f"{unmatched} rows had no metadata match."
    )
    result = result.drop('_metadata_matched')

    # save the processed result to a file
    result.write_csv(f"{output_path}/final-result.tsv", separator="\t")
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


if __name__ == "__main__":
    main()
