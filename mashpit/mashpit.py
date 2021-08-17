#!/usr/bin/env python3

import sys 
import argparse

from mashpit import __version__
from mashpit import create 
from mashpit import config
from mashpit import metadata
from mashpit import sketch
from mashpit import query
from mashpit import split

def commandToArgs(commandline):
    parser = argparse.ArgumentParser(description="""A sketch-based surveillance platform.""")
    parser.add_argument('-v', '--version',
                        action='version',
                        version='%(prog)s ' + __version__
                        )
    subparsers = parser.add_subparsers(dest="subparser_name", help=None, metavar="subcommand")
    subparsers.required = True
    # Add "create" subcommand
    description = 'Create new mashpit database'
    subparser = subparsers.add_parser("create", help=description,description=description)
    subparser.add_argument(dest='database',type=str,help="Name for the database.")
    subparser.set_defaults(func=create.create)
    # Add "config" subcommand
    description = 'Add Entrez email and key to environment variables'
    subparser = subparsers.add_parser("config",help=description,description=description)
    subparser.add_argument(dest="email", type=str, help="Entrez email address")
    subparser.add_argument("-k", "--key", help="Entrez api key")
    subparser.set_defaults(func=config.config)
    # Add "metadata" subcommad
    description = 'Collect metadata from NCBI based on bioproject/biosample accession or keywords'
    subparser = subparsers.add_parser("metadata", help=description,description=description)
    subparser.add_argument(dest="database", type=str, help="Name of the database")
    subparser.add_argument(dest="method", choices=["bioproject_list", "biosample_list", "keyword"], help= "Metadata collecting method. Available options: bioproject_list, biosample_list, keyword")
    subparser.add_argument("-l", "--list", help="File name of a list of bioproject or biosample")
    subparser.add_argument("-t", "--term", help="Query keyword")
    subparser.set_defaults(func=metadata.metadata)
    # Add "sketch" subcommand
    description = 'Build sketches for the records in the database'
    subparser = subparsers.add_parser("sketch", help=description,description=description)
    subparser.add_argument(dest="database", type=str,help="Name of the database")
    subparser.add_argument("-n","--number",type=int,default=1000,help="Number of genomes in a batch to be downloaded and sketched. Default is 1000.")
    subparser.set_defaults(func=sketch.sketch)
    # Add "split" subcommand
    description = 'Split large signature file to speed up the query'
    subparser = subparsers.add_parser("split", help=description,description=description)
    subparser.add_argument(dest="database", type=str,help="Name of the database")
    subparser.add_argument("-n","--number",type=int,default=16,help="Number of files to be splited into. Default is 16.")
    subparser.set_defaults(func=split.split)
    # Add "query" subcommand
    description = 'Find the most similar assemblies to the target sample'
    subparser = subparsers.add_parser("query",help=description,description=description)
    subparser.add_argument(dest="sample", type=str, help="target sample file path")
    subparser.add_argument(dest="database", type=str, help="name of the database")
    subparser.add_argument("-n","--number", type=int, help="number of separated signature file")
    subparser.add_argument("-f", "--force", help="overwrite query record if query table exists", action="store_true")
    subparser.set_defaults(func=query.query)

    args = parser.parse_args(commandline)
    return args


def main():
    if len(sys.argv[1:])<1:
        print("Subcommand is required to run mashpit. Use -h or --help to show help information.\n")
        print("Subcommand options:\n")
        print("create: create new mashpit database.\t")
        print("metadata: collect metadata from NCBI based on bioproject/biosample accession or keywords")
        print("sketch: build sketches for the records in the database")
        print("query: find the most similar assemblies to the target sample")
        exit(0)
    args = commandToArgs(sys.argv[1:])
    args.func(args)
    return 