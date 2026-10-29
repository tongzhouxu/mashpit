"""Strict, offline merging of compatible database snapshots."""

import json
import math
from contextlib import closing
import os
from pathlib import Path
import shutil
import sqlite3
import tempfile

import pandas as pd
from sourmash import SourmashSignature

from mashpit.build import (
    create_database,
    load_all_signatures,
    resolve_database_file,
    resolve_sig_path,
    write_signatures,
)
from mashpit.report import clean, read_cluster_members, summarize_clusters


class MergeConflict(ValueError):
    """Conflict details can be inspected without publishing a partial database."""

    def __init__(self, conflicts):
        self.conflicts = conflicts
        super().__init__("Merge conflicts: " + json.dumps(conflicts, default=str))


def combine(frames, key, label):
    """Union columns and fill missing values, rejecting ambiguous records."""
    rows = {}
    columns = list(dict.fromkeys(c for frame in frames for c in frame.columns))
    for frame in frames:
        for row in frame.to_dict("records"):
            identity = row.get(key)
            if clean(identity) == "Unknown":
                raise ValueError("%s has a missing %s" % (label, key))
            previous = rows.setdefault(str(identity), {})
            for column, value in row.items():
                if clean(value) == "Unknown":
                    continue
                old = previous.get(column)
                if old is not None and old != value:
                    raise MergeConflict(
                        [
                            dict(
                                table=label,
                                id=identity,
                                column=column,
                                values=[old, value],
                            )
                        ]
                    )
                previous[column] = value
    return pd.DataFrame(list(rows.values()), columns=columns)


def read_table(conn, table):
    if not conn.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?", (table,)
    ).fetchone():
        return pd.DataFrame()
    return pd.read_sql_query('SELECT * FROM "%s"' % table, conn)


def merge_databases(databases, output, prefixes=None):
    """Build a fresh output directory; input SQLite files are opened read-only."""
    if len(databases) < 2:
        raise ValueError("At least two input databases are required")
    output = Path(output).absolute()
    if output.exists():
        raise FileExistsError("Output already exists: %s" % output)
    paths = [resolve_database_file(p).resolve() for p in databases]
    if len(set(paths)) != len(paths):
        raise ValueError("Input databases must be distinct")
    if any(p.parent == output or output in p.parents for p in paths):
        raise ValueError("Output must not contain an input database")
    prefixes = prefixes or [""] * len(paths)
    if len(prefixes) != len(paths) or any("/" in p or "\\" in p for p in prefixes):
        raise ValueError("Supply one prefix per input, without path separators")
    sources, signatures, tables, archives = [], {}, {}, {}
    settings = None
    sketch_settings = None
    deduplicated = 0
    mappings = []
    for path, prefix in zip(paths, prefixes):
        with closing(sqlite3.connect(path.as_uri() + "?mode=ro", uri=True)) as conn:
            desc = dict(conn.execute("SELECT name, value FROM DESC"))
            kind = desc.get("Type")
            if kind not in {"Custom", "Taxonomy"}:
                raise ValueError("Only Custom and Taxonomy databases can be merged")
            supported_tables = {"DESC", "METADATA", "REPRESENTATIVE", "ANNOTATION_LOG"}
            if kind == "Taxonomy":
                supported_tables.update({"CLUSTER_MEMBERS", "REP_CLUSTER"})
            input_tables = {
                row[0]
                for row in conn.execute(
                    "SELECT name FROM sqlite_master WHERE type='table'"
                )
                if not row[0].startswith("sqlite_")
            }
            unsupported = input_tables - supported_tables
            if unsupported:
                raise ValueError(
                    "Unsupported tables would be lost during merge: %s"
                    % ", ".join(sorted(unsupported))
                )
            keys = ["Type", "Hash_number", "Kmer_size"]
            if kind == "Taxonomy":
                if any(prefixes):
                    raise ValueError("Prefixes apply only to custom databases")
                keys += ["Species", "Version", "Representative_method", "Tree_radius"]
            current = {k: desc.get(k) for k in keys}
            if any(v is None or str(v).strip() == "" for v in current.values()):
                raise ValueError("Missing required database settings: %s" % current)
            for key in ["Hash_number", "Kmer_size"]:
                current[key] = str(int(current[key]))
                if int(current[key]) <= 0:
                    raise ValueError("Sketch settings must be positive")
            if kind == "Taxonomy":
                radius = float(current["Tree_radius"])
                if not math.isfinite(radius) or radius < 0:
                    raise ValueError("Tree radius must be finite and nonnegative")
                current["Tree_radius"] = str(radius)
            if settings is not None and current != settings:
                raise ValueError(
                    "Incompatible database settings: %s versus %s" % (settings, current)
                )
            settings = current
            metadata = read_table(conn, "METADATA")
            if metadata.empty:
                raise ValueError("Input database has no metadata")
            for key in ["asm_acc", "biosample_acc"]:
                if (
                    key not in metadata
                    or metadata[key].map(clean).eq("Unknown").any()
                    or metadata[key].duplicated().any()
                ):
                    raise ValueError(
                        "Input metadata must have unique nonempty %s" % key
                    )
            original_ids = set(metadata.asm_acc)
            original_biosamples = dict(zip(metadata.asm_acc, metadata.biosample_acc))
            for key in ["asm_acc", "biosample_acc"]:
                metadata[key] = prefix + metadata[key].astype(str)
            tables.setdefault("METADATA", []).append(metadata)
            loaded = load_all_signatures(resolve_sig_path(path))
            seen = set()
            for sig in loaded:
                mh = sig.minhash
                config = (
                    mh.ksize,
                    mh.moltype,
                    mh.seed,
                    mh.num,
                    mh.scaled,
                    mh.track_abundance,
                )
                if sketch_settings is not None and config != sketch_settings:
                    raise ValueError("Incompatible signature configurations")
                sketch_settings = config
                if (
                    mh.ksize != int(desc["Kmer_size"])
                    or mh.num != int(desc["Hash_number"])
                    or mh.scaled
                ):
                    raise ValueError(
                        "Signature settings disagree with database description"
                    )
                if sig.name not in original_ids:
                    raise ValueError(
                        "Signature ID is missing from metadata: %s" % sig.name
                    )
                if sig.name in seen:
                    raise ValueError("Duplicate signature ID in input: %s" % sig.name)
                seen.add(sig.name)
                name = prefix + sig.name
                if name in signatures:
                    if signatures[name].minhash != mh:
                        raise MergeConflict(
                            [
                                dict(
                                    table="signatures",
                                    id=name,
                                    reason="different sketches",
                                )
                            ]
                        )
                    deduplicated += 1
                else:
                    signatures[name] = SourmashSignature(mh, name=name)
                mappings.append(
                    dict(
                        source=str(path),
                        original_id=sig.name,
                        output_id=name,
                        original_biosample=original_biosamples[sig.name],
                        output_biosample=prefix + str(original_biosamples[sig.name]),
                    )
                )
            if seen != original_ids:
                raise ValueError("Signature IDs do not match metadata in %s" % path)
            history = read_table(conn, "ANNOTATION_LOG")
            if not history.empty:
                if "merge_source" not in history:
                    history["merge_source"] = str(path)
                else:
                    history["merge_source"] = history["merge_source"].fillna(str(path))
                tables.setdefault("ANNOTATION_LOG", []).append(history)
            if kind == "Taxonomy":
                if "REP_CLUSTER" in input_tables:
                    labels = read_table(conn, "REP_CLUSTER")
                    if (
                        "PDS_acc" not in labels
                        or labels.PDS_acc.map(clean).eq("Unknown").any()
                        or labels.PDS_acc.duplicated().any()
                    ):
                        raise ValueError(
                            "REP_CLUSTER requires unique nonempty PDS_acc values"
                        )
                    tables.setdefault("REP_CLUSTER", []).append(labels)
                reps = read_table(conn, "REPRESENTATIVE")
                if reps.empty or set(reps.asm_acc) != original_ids:
                    raise ValueError("Representative IDs do not match metadata")
                # Build-time file paths and attempts belong to each source, not the merged snapshot.
                tables.setdefault("REPRESENTATIVE", []).append(
                    reps.reindex(
                        columns=["asm_acc", "PDT_acc", "PDS_acc", "tree_radius"]
                    )
                )
                clusters = set(metadata.PDS_acc)
                members = read_table(conn, "CLUSTER_MEMBERS")
                if members.empty:
                    members = read_cluster_members(conn, clusters)
                if members.empty or not clusters <= set(members.PDS_acc):
                    raise ValueError(
                        "Taxon merge requires complete cluster membership for represented clusters"
                    )
                representatives_by_id = reps.set_index("asm_acc").to_dict("index")
                member_pairs = set(zip(members.target_acc, members.PDS_acc))
                for row in metadata.to_dict("records"):
                    representative = representatives_by_id[row["asm_acc"]]
                    target = row.get("PDT_acc")
                    if (
                        clean(row.get("PDS_acc")) == "Unknown"
                        or clean(target) == "Unknown"
                        or representative["PDS_acc"] != row["PDS_acc"]
                        or representative["PDT_acc"] != target
                        or float(representative["tree_radius"])
                        != float(current["Tree_radius"])
                        or (target, row["PDS_acc"]) not in member_pairs
                    ):
                        raise ValueError(
                            "Inconsistent taxon representative or membership metadata"
                        )
                tables.setdefault("CLUSTER_MEMBERS", []).append(members)
                for archive in sorted((path.parent / "membership").glob("*.tsv.gz")):
                    if not archive.name.startswith(desc["Version"] + "."):
                        continue
                    archives.setdefault(archive.name, []).append(
                        pd.read_csv(archive, sep="\t", dtype=str)
                    )
            sources.append(
                dict(
                    database=str(path),
                    description=desc,
                    prefix=prefix,
                    representatives=(
                        reps.astype(object)
                        .where(pd.notna(reps), None)
                        .to_dict("records")
                        if kind == "Taxonomy"
                        else []
                    ),
                )
            )
    merged = {}
    for table, frames in tables.items():
        if table == "ANNOTATION_LOG":
            merged[table] = pd.concat(frames, ignore_index=True)
        else:
            key = {
                "METADATA": "asm_acc",
                "REPRESENTATIVE": "asm_acc",
                "CLUSTER_MEMBERS": "target_acc",
                "REP_CLUSTER": "PDS_acc",
            }[table]
            merged[table] = combine(frames, key, table)
    if merged["METADATA"].biosample_acc.duplicated().any():
        raise MergeConflict(
            [
                dict(
                    table="METADATA",
                    column="biosample_acc",
                    reason="BioSample IDs collide across different assemblies; use custom prefixes",
                )
            ]
        )
    archive_frames = {
        name: combine(frames, "target_acc", name) for name, frames in archives.items()
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(
        tempfile.mkdtemp(prefix="." + output.name + "-merge-", dir=str(output.parent))
    )
    try:
        with sqlite3.connect(str(stage / (output.name + ".db"))) as conn:
            create_database(conn)
            for table, frame in merged.items():
                if table == "ANNOTATION_LOG":
                    frame = frame.drop(columns=["id"], errors="ignore")
                    conn.execute(
                        "ALTER TABLE ANNOTATION_LOG ADD COLUMN merge_source TEXT"
                    )
                    frame.to_sql(table, conn, index=False, if_exists="append")
                else:
                    frame.to_sql(table, conn, index=False, if_exists="replace")
            conn.execute(
                "CREATE UNIQUE INDEX metadata_biosample ON METADATA(biosample_acc)"
            )
            conn.execute("CREATE UNIQUE INDEX metadata_assembly ON METADATA(asm_acc)")
            if "CLUSTER_MEMBERS" in merged:
                conn.execute(
                    "CREATE INDEX cluster_members_pds ON CLUSTER_MEMBERS(PDS_acc)"
                )
            if "REP_CLUSTER" in merged:
                conn.execute(
                    "CREATE UNIQUE INDEX rep_cluster_pds ON REP_CLUSTER(PDS_acc)"
                )
            desc = dict(settings)
            desc[
                (
                    "Sequence_count"
                    if settings["Type"] == "Custom"
                    else "Representative_count"
                )
            ] = str(len(signatures))
            conn.execute("DELETE FROM DESC")
            conn.executemany("INSERT INTO DESC VALUES (?, ?)", sorted(desc.items()))
            if conn.execute("PRAGMA integrity_check").fetchone()[0] != "ok":
                raise ValueError("Output database failed integrity check")
        write_signatures(
            [signatures[k] for k in sorted(signatures)], stage / (output.name + ".sig")
        )
        if settings["Type"] == "Taxonomy":
            merged["REPRESENTATIVE"].to_csv(
                stage / "representatives.tsv", sep="\t", index=False
            )
            summarize_clusters(merged["CLUSTER_MEMBERS"]).to_csv(
                stage / "cluster_summary.tsv", sep="\t", index=False
            )
        for index, path in enumerate(paths, 1):
            source_dir = stage / "merge_sources" / str(index)
            source_dir.mkdir(parents=True)
            for artifact in sorted(path.parent.iterdir()):
                if artifact.is_file() and artifact.name in {
                    "representatives.tsv",
                    "cluster_summary.tsv",
                    "unavailable_assemblies.tsv",
                    "unavailable_sequences.tsv",
                    "merge_provenance.json",
                    "database.json",
                    "rep_cluster_annotations.json",
                }:
                    shutil.copy2(artifact, source_dir / artifact.name)
        for name, frame in archive_frames.items():
            (stage / "membership").mkdir(exist_ok=True)
            frame.to_csv(
                stage / "membership" / name, sep="\t", index=False, compression="gzip"
            )
        provenance = dict(
            sources=sources,
            id_mappings=mappings,
            deduplicated_signatures=deduplicated,
            signature_count=len(signatures),
            sketch_configuration=sketch_settings,
        )
        (stage / "merge_provenance.json").write_text(
            json.dumps(provenance, indent=2, allow_nan=False) + "\n"
        )
        # Reserve the destination so a concurrent merge cannot replace an existing output.
        output.mkdir()
        try:
            os.replace(stage, output)
        except BaseException:
            output.rmdir()
            raise
    finally:
        if stage.exists():
            shutil.rmtree(stage)
    return provenance


def merge(args):
    try:
        result = merge_databases(args.databases, args.output, args.prefixes)
    except MergeConflict as exc:
        if args.conflict_report:
            with Path(args.conflict_report).open("x") as handle:
                json.dump(exc.conflicts, handle, indent=2, default=str)
        raise
    print("Merged %d signatures into %s" % (result["signature_count"], args.output))
