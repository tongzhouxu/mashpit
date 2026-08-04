#!/usr/bin/env python3
"""Add or update custom METADATA columns on an already-built database.

Works on any Mashpit database (taxon, accession, or custom), since the
METADATA table schema is identical regardless of build type, and on
databases built before this feature existed - the ANNOTATION_LOG table
that records change history is created on first use if it isn't already
present.
"""

from mashpit.build import (
    annotate_database,
    list_accessions,
    read_annotation_log,
    resolve_database_file,
)


def print_accessions(db_path):
    print("asm_acc\tbiosample_acc")
    for asm_acc, biosample_acc in list_accessions(db_path):
        print("%s\t%s" % (asm_acc, biosample_acc))


def print_history(db_path):
    rows = read_annotation_log(db_path)
    if not rows:
        print("No annotation history recorded for %s" % db_path)
        return

    for (
        timestamp,
        key_column,
        columns_added,
        columns_updated,
        values_file,
        values_file_sha256,
        rows_updated,
        rows_unmatched_in_file,
        rows_unmatched_in_database,
    ) in rows:
        print(
            "%s  key=%s  file=%s (sha256 %s)"
            % (timestamp, key_column, values_file, values_file_sha256[:12])
        )
        if columns_added:
            print("  columns added   : %s" % columns_added)
        print("  columns updated : %s" % columns_updated)
        print(
            "  rows updated: %d, unmatched in file: %d, unmatched in database: %d"
            % (rows_updated, rows_unmatched_in_file, rows_unmatched_in_database)
        )


def annotate(args):
    db_path = resolve_database_file(args.database)

    if args.list_accessions:
        print_accessions(db_path)
        return

    if args.history:
        print_history(db_path)
        return

    if not args.values:
        raise SystemExit(
            "--values is required unless --history or --list-accessions is given"
        )

    summary = annotate_database(db_path, args.values, getattr(args, "key", None))

    print("Database annotated: %s" % db_path)
    if summary["columns_added"]:
        print("New columns added   : %s" % ", ".join(summary["columns_added"]))
    print("Columns updated     : %s" % ", ".join(summary["columns_updated"]))
    print("Rows updated        : %d" % summary["rows_updated"])

    if summary["unmatched_in_file"]:
        print(
            "Warning: %d id(s) in --values had no matching database row: %s"
            % (
                len(summary["unmatched_in_file"]),
                ", ".join(summary["unmatched_in_file"][:10]),
            )
        )
    if summary["unmatched_in_database"]:
        print(
            "Note: %d database row(s) were not mentioned in --values and "
            "keep their existing values" % len(summary["unmatched_in_database"])
        )
