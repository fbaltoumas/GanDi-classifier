import argparse as ap
import sys

from . import __version__, classifier
from . import download_db

DESCRIPTION = (
    f"gandi {__version__}: contig/genome classifier and reference-database downloader for the\n"
    "Global Anaerobic Digestion (GanDi) database.\n"
    "\n"
    "Classify a genome against the database:\n"
    "  gandi classifier -i genome.fasta -o output_dir -c plasmids -d /path/to/database\n"
    "\n"
    "Download and build the reference database:\n"
    "  gandi download-db -o /path/to/database"
)


def _build_parser():
    parser = ap.ArgumentParser(
        prog="gandi",
        description=DESCRIPTION,
        formatter_class=ap.RawDescriptionHelpFormatter,
    )
    parser.add_argument("-v", "--version", action="version", version=f"gandi {__version__}")
    subparsers = parser.add_subparsers(dest="command", metavar="<subcommand>")
    subparsers.add_parser(
        "classifier",
        help="Classify contigs/genomes against a GanDi reference database (ANI/AAI workflows)."
    )
    subparsers.add_parser(
        "download-db",
        help="Download and build the GanDi reference database."
    )
    return parser


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    parser = _build_parser()

    if not argv:
        parser.print_help()
        sys.exit(2)

    if argv[0] in ("-h", "--help"):
        parser.print_help()
        sys.exit(0)

    if argv[0] in ("-v", "--version"):
        print(f"gandi {__version__}")
        sys.exit(0)

    command, rest = argv[0], argv[1:]
    if command == "classifier":
        classifier.main(rest, prog="gandi classifier")
    elif command == "download-db":
        download_db.main(rest, prog="gandi download-db")
    else:
        parser.print_help()
        print(f"\nUnknown subcommand: {command!r}", file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()
