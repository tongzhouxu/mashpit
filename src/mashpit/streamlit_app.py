#!/usr/bin/env python3

import os
import re
import sqlite3
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
import streamlit as st


APP_NAME = "Mashpit Explorer"


def safe_filename(filename):
    basename = os.path.basename(filename.replace("\\", "/"))
    cleaned = re.sub(r"[^A-Za-z0-9._-]", "_", basename)
    return cleaned or "query.fasta"


def validate_database(database_text):
    if not database_text.strip():
        return None, "Select a Mashpit database directory."

    database = Path(database_text).expanduser().resolve()
    if not database.is_dir():
        return None, "The database directory does not exist."

    database_files = list(database.glob("*.db"))
    signature_files = list(database.glob("*.sig"))
    if len(database_files) != 1 or len(signature_files) != 1:
        return None, "The directory must contain exactly one .db and one .sig file."

    if database_files[0].stem != signature_files[0].stem:
        return None, "The database and signature filenames do not match."

    return database, None


def read_optional_bytes(path):
    return path.read_bytes() if path.is_file() else None


def read_database_summary(database):
    sql_path = next(database.glob("*.db"))
    conn = sqlite3.connect(f"file:{sql_path}?mode=ro", uri=True)
    try:
        rows = dict(conn.execute("SELECT name, value FROM DESC").fetchall())
    finally:
        conn.close()
    return rows


def run_query(uploaded_assembly, database, number, threshold, annotation, tie_tolerance_hashes):
    with tempfile.TemporaryDirectory(prefix="mashpit-query-") as temporary:
        work_dir = Path(temporary)
        assembly_name = safe_filename(uploaded_assembly.name)
        assembly_path = work_dir / assembly_name
        assembly_path.write_bytes(uploaded_assembly.getvalue())

        command = [
            sys.executable,
            "-m",
            "mashpit.mashpit",
            "query",
            str(assembly_path),
            str(database),
            "--number",
            str(number),
            "--threshold",
            str(threshold),
            "--tie-tolerance-hashes",
            str(tie_tolerance_hashes),
        ]
        if annotation.strip():
            command.extend(["--annotation", annotation.strip()])

        completed = subprocess.run(
            command,
            cwd=work_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        # query_name mirrors mashpit.query's own naming: split on the first
        # "." in the uploaded filename (e.g. "GCA_1.1_genomic.fna" -> "GCA_1").
        query_name = assembly_name.split(".")[0]
        representative_path = work_dir / f"{query_name}_representative_matches.csv"

        if not representative_path.is_file():
            # The representative-matches CSV is written before mashtree
            # generation, so a missing CSV means a genuine failure (bad
            # database, sketching error, etc.), not just "no tree could be
            # built".
            details = completed.stderr.strip() or completed.stdout.strip()
            raise RuntimeError(
                details or "Mashpit query failed without producing any output."
            )

        representative_df = pd.read_csv(representative_path, index_col=0)

        # Only taxon databases have SNP clusters to group by; accession
        # databases never get a cluster_candidates file.
        cluster_path = work_dir / f"{query_name}_cluster_candidates.csv"
        cluster_df = None
        cluster_csv = None
        if cluster_path.is_file():
            cluster_df = pd.read_csv(cluster_path, index_col=0)
            cluster_csv = cluster_path.read_bytes()

        tree_png = read_optional_bytes(work_dir / f"{query_name}_tree.png")
        tree_skip_reason = None
        if tree_png is None:
            # generate_mashtree exits(1) - and mashpit query with it - when
            # the top hit is below the threshold or fewer than two hits
            # qualify. That is a legitimate outcome, not a crash, so the
            # already-written CSVs are still shown.
            tree_skip_reason = (
                completed.stderr.strip()
                or "A tree was not generated: the top hit may be below the "
                "similarity threshold, or fewer than two candidates qualified."
            )

        return {
            "representative_df": representative_df,
            "representative_csv": representative_path.read_bytes(),
            "cluster_df": cluster_df,
            "cluster_csv": cluster_csv,
            "tree_png": tree_png,
            "tree_newick": read_optional_bytes(work_dir / f"{query_name}_tree.newick"),
            "tree_skip_reason": tree_skip_reason,
            "query_name": query_name,
            "log": completed.stdout + completed.stderr,
        }


def display_results(results, db_summary):
    representative_df = results["representative_df"]
    cluster_df = results["cluster_df"]

    columns = st.columns(3)
    columns[0].metric("Candidate hits", len(representative_df))
    top_score = (
        representative_df["similarity_score"].max()
        if "similarity_score" in representative_df
        else None
    )
    columns[1].metric(
        "Top similarity", f"{top_score:.3f}" if top_score is not None else "n/a"
    )
    columns[2].metric("Database type", db_summary.get("Type", "n/a"))

    if cluster_df is not None:
        # Cluster candidates are the headline view: multiple representatives
        # of the same SNP cluster collapse into one row (hits_in_results /
        # total_representatives), ranked by best_similarity_score - never by
        # hit count alone. near_top flags clusters whose best hit is within
        # the sketch's own resolution of the single best score, since at
        # typical sketch sizes a gap like 0.999 vs 0.998 is one hash of
        # difference and not a meaningful distinction.
        near_top_count = int(cluster_df["near_top"].sum())
        st.subheader("Candidate clusters")
        st.caption(
            f"{near_top_count} cluster(s) are statistically tied with the top hit "
            "at this database's sketch resolution (see the near_top column)."
        )
        st.dataframe(cluster_df, use_container_width=True, hide_index=True)
        st.download_button(
            "Download cluster candidates",
            results["cluster_csv"],
            file_name=f"{results['query_name']}_cluster_candidates.csv",
            mime="text/csv",
        )

        with st.expander("Representative-level evidence"):
            st.dataframe(representative_df, use_container_width=True, hide_index=True)
            st.download_button(
                "Download representative matches",
                results["representative_csv"],
                file_name=f"{results['query_name']}_representative_matches.csv",
                mime="text/csv",
            )
    else:
        # Accession databases have no SNP clusters to group by.
        st.subheader("Query results")
        st.dataframe(representative_df, use_container_width=True, hide_index=True)
        st.download_button(
            "Download results",
            results["representative_csv"],
            file_name=f"{results['query_name']}_representative_matches.csv",
            mime="text/csv",
        )

    if results["tree_png"] is not None:
        st.subheader("Candidate sketch-distance tree")
        st.image(results["tree_png"], use_container_width=True)
        if results["tree_newick"] is not None:
            st.download_button(
                "Download Newick tree",
                results["tree_newick"],
                file_name=f"{results['query_name']}_tree.newick",
                mime="text/plain",
            )
    else:
        st.info(results["tree_skip_reason"])

    with st.expander("Query log"):
        st.code(results["log"] or "No console output was produced.")


def main():
    st.set_page_config(page_title=APP_NAME, page_icon="🧬", layout="wide")
    st.title(APP_NAME)
    st.write(
        "Screen an assembled genome against a local Mashpit database. "
        "Your genome and results remain on this computer."
    )

    with st.sidebar:
        st.header("Query settings")
        number = st.number_input(
            "Maximum representative hits",
            min_value=1,
            max_value=1000,
            value=50,
            step=1,
            help="How many representative genome hits to consider, before "
            "they get grouped into cluster candidates.",
        )
        threshold = st.slider(
            "Tree similarity threshold",
            min_value=0.0,
            max_value=1.0,
            value=0.85,
            step=0.01,
        )
        tie_tolerance_hashes = st.number_input(
            "Tie tolerance (sketch hashes)",
            min_value=0,
            max_value=100,
            value=2,
            step=1,
            help="A cluster is flagged near_top when its best hit is within "
            "this many sketch hashes of the single best hit overall - e.g. "
            "at a 1000-hash sketch, 2 hashes is a 0.002 similarity gap, "
            "which is within normal MinHash sampling noise.",
        )
        annotation = st.text_input(
            "Optional tree annotation",
            help="Metadata column such as isolation_source or geo_loc_name.",
        )

    with st.form("mashpit-query-form"):
        database_text = st.text_input(
            "Mashpit database directory",
            placeholder="/path/to/database",
            help="Choose the directory containing one matching .db and .sig file.",
        )
        uploaded_assembly = st.file_uploader(
            "Query assembly",
            type=["fa", "fasta", "fna", "fas", "gz"],
            help="Upload an assembled genome in FASTA format.",
        )
        submitted = st.form_submit_button("Run Mashpit", type="primary")

    if submitted:
        database, database_error = validate_database(database_text)
        if database_error:
            st.error(database_error)
        elif uploaded_assembly is None:
            st.error("Upload a query assembly before running Mashpit.")
        else:
            try:
                with st.spinner("Comparing the query with database representatives..."):
                    st.session_state["mashpit_results"] = run_query(
                        uploaded_assembly,
                        database,
                        int(number),
                        float(threshold),
                        annotation,
                        int(tie_tolerance_hashes),
                    )
                    st.session_state["mashpit_db_summary"] = read_database_summary(
                        database
                    )
            except Exception as error:
                st.session_state.pop("mashpit_results", None)
                st.session_state.pop("mashpit_db_summary", None)
                st.error(str(error))

    results = st.session_state.get("mashpit_results")
    if results is not None:
        display_results(results, st.session_state.get("mashpit_db_summary", {}))


if __name__ == "__main__":
    main()
