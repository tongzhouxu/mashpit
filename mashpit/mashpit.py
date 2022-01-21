#!/usr/bin/env python3

import sys 
import argparse

from mashpit import __version__
from mashpit import build 
from mashpit import config
from mashpit import sketch
from mashpit import query
from mashpit import split

def commandToArgs(commandline):
    parser = argparse.ArgumentParser(description="A sketch-based surveillance platform.")
    parser.add_argument('-v', '--version',
                        action='version',
                        version='%(prog)s ' + __version__
                        )
    subparsers = parser.add_subparsers(dest="subparser_name", help=None, metavar="subcommand")
    subparsers.required = True

    # "build" subcommand: build a standard or custom database
    description = 'Build mashpit database'
    subparser = subparsers.add_parser("build", help=description,description=description)
    subparser.add_argument(dest="type", choices=["standard", "biosample_list", "keyword"], help= "Database type: standard or custom")
    subparser.add_argument(dest="name", type=str, help="Mashpit database name")
    subparser.add_argument("-s", "--species", default="Salmonella",help="Query keyword")
    subparser.add_argument("-l", "--list", help="File containing a list of biosample accessions")
    subparser.add_argument("-t", "--term", help="Query keyword")
    subparser.set_defaults(func=build.build)

    # "config" subcommand
    description = 'Add Entrez email and key to environment variables'
    subparser = subparsers.add_parser("config",help=description,description=description)
    subparser.add_argument(dest="email", type=str, help="Entrez email address")
    subparser.add_argument("-k", "--key", help="Entrez api key")
    subparser.set_defaults(func=config.config)
    
    # "sketch" subcommand
    description = 'Build sketches for the records in the database'
    subparser = subparsers.add_parser("sketch", help=description,description=description)
    subparser.add_argument(dest="name", type=str,help="Mashpit database name")
    subparser.set_defaults(func=sketch.sketch)
    
    # "query" subcommand
    description = 'Find the most similar assemblies to the target sample'
    subparser = subparsers.add_parser("query",help=description,description=description)
    subparser.add_argument(dest="sample", type=str, help="file name of the query sample")
    subparser.add_argument(dest="database", type=str, help="name of the database")
    subparser.add_argument("-n","--number", type=int, help="number of splited signature file")
    subparser.set_defaults(func=query.query)

    args = parser.parse_args(commandline)
    return args


def main():
    if len(sys.argv[1:])<1:
        print("Subcommand is required to run mashpit. Use -h or --help to show help information.\n")
        print("Subcommand options:\n")
        
        exit(0)
    args = commandToArgs(sys.argv[1:])
    args.func(args)
    return 