#!/usr/bin/env python3
"""Rewrite an already-built database's signature storage into (or out of)
shards that can be loaded and compared in parallel at query time - see
mashpit.build.reshard_database. Works on any Mashpit database (taxon,
accession, or custom) and touches only the .sig file/directory; the
database's .db file is never modified.
"""

from mashpit.build import reshard_database


def reshard(args):
    summary = reshard_database(args.database, getattr(args, "shards", None))
    if summary["shards"] == 1:
        print(
            "Database merged into a single signature file: %d signature(s)"
            % summary["signature_count"]
        )
    else:
        print(
            "Database resharded: %d signature(s) across %d shard(s)"
            % (summary["signature_count"], summary["shards"])
        )
