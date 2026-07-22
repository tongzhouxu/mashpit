#!/usr/bin/env python3

import argparse
import sys

from mashpit import __version__
from mashpit import build
from mashpit import gui
from mashpit import query


def positive_int(value):
    value = int(value)
    if value < 1:
        raise argparse.ArgumentTypeError("value must be at least 1")
    return value


def nonnegative_float(value):
    value = float(value)
    if value < 0:
        raise argparse.ArgumentTypeError("value must be non-negative")
    return value


def commandToArgs(commandline):
    parser = argparse.ArgumentParser(
        description="A lightweight, sketch-based genomic surveillance platform."
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s " + __version__,
    )

    subparsers = parser.add_subparsers(dest="command", metavar="command")
    subparsers.required = True

    subparser_build = subparsers.add_parser(
        "build",
        help="Build a Mashpit database",
    )
    subparser_query = subparsers.add_parser(
        "query",
        help="Query for the most similar isolates",
    )
    subparser_gui = subparsers.add_parser(
        "gui",
        help="Launch the Streamlit GUI",
    )

    # Build arguments
    subparser_build.add_argument(
        "type",
        choices=["taxon", "accession"],
        help="database type",
    )
    subparser_build.add_argument(
        "name",
        help="database name",
    )
    subparser_build.add_argument(
        "--quiet",
        action="store_true",
        help="disable console logs",
    )
    subparser_build.add_argument(
        "--number",
        type=positive_int,
        default=1000,
        help="maximum number of hashes for sourmash (default: 1000)",
    )
    subparser_build.add_argument(
        "--ksize",
        type=positive_int,
        default=31,
        help="k-mer size for sourmash (default: 31)",
    )
    subparser_build.add_argument(
        "--species",
        help="NCBI Pathogen Detection taxon name; required for taxon builds",
    )
    subparser_build.add_argument(
        "--email",
        help="Entrez email; required for accession builds",
    )
    subparser_build.add_argument(
        "--key",
        help="NCBI API key",
    )
    subparser_build.add_argument(
        "--pd_version",
        help="Pathogen Detection release accession; default is the latest release",
    )
    subparser_build.add_argument(
        "--list",
        help="path to a file containing NCBI BioSample accessions",
    )
    subparser_build.add_argument(
        "--radius",
        type=nonnegative_float,
        default=20.0,
        help=(
            "maximum SNP-tree distance from an isolate to its representative "
            "for taxon builds (default: 20)"
        ),
    )
    subparser_build.add_argument(
        "--download-attempts",
        dest="download_attempts",
        type=positive_int,
        default=3,
        help="maximum assembly download attempts (default: 3)",
    )
    subparser_build.add_argument(
        "--download-batch-size",
        dest="download_batch_size",
        type=positive_int,
        default=500,
        help="number of assemblies per NCBI Datasets batch (default: 500)",
    )
    subparser_build.add_argument(
        "--retry-delay",
        dest="retry_delay",
        type=nonnegative_float,
        default=5.0,
        help="seconds between assembly download attempts (default: 5)",
    )
    subparser_build.add_argument(
        "--max-reselection-rounds",
        dest="max_reselection_rounds",
        type=positive_int,
        default=5,
        help=(
            "maximum representative reselection rounds after unavailable "
            "assemblies are excluded (default: 5)"
        ),
    )
    subparser_build.set_defaults(func=build.build)

    # Query arguments
    subparser_query.add_argument(
        "sample",
        help="path to the query assembly",
    )
    subparser_query.add_argument(
        "database",
        help="path to the database folder",
    )
    subparser_query.add_argument(
        "--number",
        type=positive_int,
        default=200,
        help="number of isolates in the query output (default: 200)",
    )
    subparser_query.add_argument(
        "--threshold",
        type=float,
        default=0.85,
        help="minimum Jaccard similarity for mashtree (default: 0.85)",
    )
    subparser_query.add_argument(
        "--annotation",
        help="metadata field used to annotate mashtree tips",
    )
    subparser_query.add_argument(
        "--tie-tolerance-hashes",
        dest="tie_tolerance_hashes",
        type=positive_int,
        default=2,
        help=(
            "sketch hashes of similarity slack used to flag a cluster as "
            "statistically tied with the top hit in the cluster candidates "
            "table, for taxon databases (default: 2)"
        ),
    )
    subparser_query.set_defaults(func=query.query)

    # GUI arguments
    subparser_gui.add_argument(
        "--port",
        type=positive_int,
        help="port for the Streamlit server (default: Streamlit's own default, 8501)",
    )
    subparser_gui.set_defaults(func=gui.gui)

    return parser.parse_args(commandline)


def main():
    if len(sys.argv) == 1:
        print(
            "Subcommand is required to run mashpit. "
            "Use -h or --help to show help information.\n"
        )
        print("Subcommand options:\n")
        return

    args = commandToArgs(sys.argv[1:])
    args.func(args)


if __name__ == "__main__":
    main()
