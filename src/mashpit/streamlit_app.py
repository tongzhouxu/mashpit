#!/usr/bin/env python3

import json
import os
import re
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


def run_query(uploaded_assembly, database, number, threshold, annotation):
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
        if completed.returncode != 0:
            details = completed.stderr.strip() or completed.stdout.strip()
            raise RuntimeError(
                details or "Mashpit query failed without an error message."
            )

        query_name = assembly_name.split(".")[0]
        candidate_path = work_dir / f"{query_name}_cluster_candidates.csv"
        representative_path = work_dir / f"{query_name}_representative_matches.csv"
        summary_path = work_dir / f"{query_name}_run_summary.json"
        if not candidate_path.is_file() or not representative_path.is_file():
            raise RuntimeError("Mashpit did not produce the expected result tables.")
        if not summary_path.is_file():
            raise RuntimeError("Mashpit did not produce a run summary.")

        candidate_bytes = candidate_path.read_bytes()
        representative_bytes = representative_path.read_bytes()
        summary_bytes = summary_path.read_bytes()
        return {
            "candidates": pd.read_csv(candidate_path),
            "representatives": pd.read_csv(representative_path),
            "summary": json.loads(summary_bytes),
            "candidate_csv": candidate_bytes,
            "representative_csv": representative_bytes,
            "summary_json": summary_bytes,
            "tree_png": read_optional_bytes(work_dir / f"{query_name}_tree.png"),
            "tree_newick": read_optional_bytes(
                work_dir / f"{query_name}_tree.newick"
            ),
            "query_name": query_name,
            "log": completed.stdout + completed.stderr,
        }


def display_results(results):
    summary = results["summary"]
    st.success(summary["conclusion"])
    st.caption(summary["interpretation_warning"])

    leading = summary.get("leading_candidates", [])
    columns = st.columns(3)
    columns[0].metric("Leading candidates", len(leading))
    columns[1].metric("Sketch size", summary["sketch_size"])
    columns[2].metric("K-mer size", summary["kmer_size"])

    st.subheader("Cluster candidates")
    st.dataframe(results["candidates"], width="stretch", hide_index=True)
    st.download_button(
        "Download cluster candidates",
        results["candidate_csv"],
        file_name=f"{results['query_name']}_cluster_candidates.csv",
        mime="text/csv",
    )

    with st.expander("Representative evidence"):
        st.dataframe(results["representatives"], width="stretch", hide_index=True)
        st.download_button(
            "Download representative evidence",
            results["representative_csv"],
            file_name=f"{results['query_name']}_representative_matches.csv",
            mime="text/csv",
        )

    if results["tree_png"] is not None:
        st.subheader("Candidate sketch-distance tree")
        st.image(results["tree_png"], width="stretch")
        if results["tree_newick"] is not None:
            st.download_button(
                "Download Newick tree",
                results["tree_newick"],
                file_name=f"{results['query_name']}_tree.newick",
                mime="text/plain",
            )
    else:
        st.info("A tree was not generated because fewer than two candidates qualified.")

    st.download_button(
        "Download run summary",
        results["summary_json"],
        file_name=f"{results['query_name']}_run_summary.json",
        mime="application/json",
    )
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
            "Maximum candidate clusters",
            min_value=1,
            max_value=1000,
            value=50,
            step=1,
        )
        threshold = st.slider(
            "Tree similarity threshold",
            min_value=0.0,
            max_value=1.0,
            value=0.85,
            step=0.01,
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
                    )
            except Exception as error:
                st.session_state.pop("mashpit_results", None)
                st.error(str(error))

    results = st.session_state.get("mashpit_results")
    if results is not None:
        display_results(results)


if __name__ == "__main__":
    main()
