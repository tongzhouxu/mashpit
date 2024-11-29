#!/usr/bin/env python3
import sys
import argparse

from mashpit import __version__
from mashpit import build
from mashpit import update
from mashpit import query
from mashpit import webserver


def commandToArgs(commandline):
    parser = argparse.ArgumentParser(
        description="A sketch-based surveillance platform."
    )
    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s " + __version__
    )

    # subparsers for mashpit
    subparsers = parser.add_subparsers(metavar="command")
    subparsers.required = True
    subparser_build = subparsers.add_parser("build", help="Build mashpit database")
    subparser_update = subparsers.add_parser("update", help="Update a database")
    subparser_query = subparsers.add_parser(
        "query", help="Query for the most similar isolates"
    )
    subparser_webserver = subparsers.add_parser(
        "webserver", help="Start a local web server"
    )

    # arguments for mashpit build
    subparser_build.add_argument(
        "type", type=str, help="mashpit database type", choices=["taxon", "accession"]
    )
    subparser_build.add_argument("name", type=str, help="mashpit database name")
    subparser_build.add_argument("--quiet", action="store_true", help="disable logs")
    subparser_build.add_argument(
        "--number",
        type=int,
        default=1000,
        help="maximum number of hashes for sourmash, default is 1000",
    )
    subparser_build.add_argument(
        "--ksize", type=int, default=31, help="kmer size for sourmash, default is 31"
    )
    subparser_build.add_argument("--species", type=str, help="species name")
    subparser_build.add_argument("--email", type=str, help="Entrez email")
    subparser_build.add_argument("--key", type=str, help="Entrez api key")
    subparser_build.add_argument(
        "--pd_version",
        type=str,
        help="a specified Pathogen Detection version (PDG accession). Default is the latest.",
    )
    subparser_build.add_argument(
        "--list", type=str, help="Path to a list of NCBI BioSample accessions"
    )
    subparser_build.set_defaults(func=build.build)

    # "update" subcommand
    # TODO: remove the name argument? because the name is already in the database folder
    subparser_update.add_argument(
        dest="database", type=str, help="path for the database folder"
    )
    subparser_update.add_argument(dest="name", type=str, help="database name")
    subparser_update.add_argument(
        "--metadata", type=str, help="metadata file in csv format"
    )
    subparser_update.add_argument("--quiet", action="store_true", help="disable logs")
    subparser_update.set_defaults(func=update.update)

    # "query" subcommand
    subparser_query.add_argument(
        dest="sample", type=str, help="file path to the query sample"
    )
    subparser_query.add_argument(
        dest="database", type=str, help="path to the database folder"
    )
    subparser_query.add_argument(
        "--number",
        type=int,
        help="number of isolates in the query output, default is 200",
        default=200,
    )
    subparser_query.add_argument(
        "--threshold",
        type=float,
        help="minimum jaccard similarity for mashtree, default is 0.85",
        default=0.85,
    )
    subparser_query.add_argument(
        "--annotation", type=str, help="mashtree tip annoatation, default is none"
    )
    subparser_query.set_defaults(func=query.query)

    # "webserver" subcommand
    subparser_webserver.set_defaults(func=webserver.webserver)

    args = parser.parse_args(commandline)
    return args


def main():
    if len(sys.argv[1:]) < 1:
        print(
            "Subcommand is required to run mashpit. Use -h or --help to show help information.\n"
        )
        print("Subcommand options:\n")

        exit(0)
    args = commandToArgs(sys.argv[1:])
    args.func(args)
    return
