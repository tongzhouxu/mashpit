#!/usr/bin/env python3

import base64
import os
import re
import sqlite3
import subprocess
import sys
import tempfile
from pathlib import Path

import altair as alt
import pandas as pd

from mashpit.report import prepare_members, read_cluster_members
from mashpit.tree_plot import render_tree
import streamlit as st
import streamlit.components.v1 as components


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

        members = pd.DataFrame()
        if cluster_df is not None:
            sql_path = next(database.glob("*.db"))
            with sqlite3.connect(f"file:{sql_path}?mode=ro", uri=True) as conn:
                members = read_cluster_members(conn, cluster_df["PDS_acc"])

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
            "members": members,
            "cluster_csv": cluster_csv,
            "tree_png": tree_png,
            "tree_newick": read_optional_bytes(work_dir / f"{query_name}_tree.newick"),
            "tree_skip_reason": tree_skip_reason,
            "query_name": query_name,
            "log": completed.stdout + completed.stderr,
        }


def metadata_chart(frame, column, title):
    counts = frame[column].fillna("Unknown").value_counts().rename_axis("Value").reset_index(name="Isolates")
    # Keep totals honest while limiting long legends / category lists.
    if len(counts) > 10:
        remaining = int(counts.iloc[10:]["Isolates"].sum())
        counts = pd.concat([counts.head(10), pd.DataFrame([{"Value": "Other categories", "Isolates": remaining}])])
    chart = alt.Chart(counts).mark_bar(color="#238b8d", cornerRadiusEnd=3).encode(
        x=alt.X("Isolates:Q", title="Isolates", axis=alt.Axis(tickMinStep=1)),
        y=alt.Y("Value:N", sort="-x", title=None, axis=alt.Axis(labelLimit=260)),
        tooltip=["Value:N", "Isolates:Q"],
    ).properties(title=title, height=max(260, len(counts) * 32 + 140))
    st.altair_chart(chart, use_container_width=True)


def display_cluster_context(cluster_df, members):
    st.subheader("Cluster landscape")
    st.caption(
        "Ranked by best representative similarity. Cluster sizes describe all isolates "
        "in the database release; metadata do not establish the source of the query isolate."
    )
    if "cluster_size" in cluster_df and cluster_df["cluster_size"].notna().any():
        plot = cluster_df.head(20).copy()
        counts = plot.melt(
            id_vars=["PDS_acc", "best_similarity_score"],
            value_vars=["environmental_count", "clinical_count", "unknown_source_count"],
            var_name="Source", value_name="Isolates",
        )
        counts["Source"] = counts["Source"].map({
            "environmental_count": "Environmental / other",
            "clinical_count": "Clinical", "unknown_source_count": "Unknown / other",
        })
        st.altair_chart(alt.Chart(counts).mark_bar().encode(
            y=alt.Y("PDS_acc:N", sort=plot["PDS_acc"].tolist(), title=None),
            x=alt.X("Isolates:Q", title="Full cluster size (isolates)"),
            color=alt.Color("Source:N", scale=alt.Scale(
                domain=["Environmental / other", "Clinical", "Unknown / other"],
                range=["#238b8d", "#e5a44e", "#b9c2cb"],
            ), legend=alt.Legend(orient="bottom")),
            tooltip=["PDS_acc:N", "Source:N", "Isolates:Q", alt.Tooltip("best_similarity_score:Q", format=".4f")],
        ).properties(height=max(300, len(plot) * 36 + 150)), use_container_width=True)
        if len(cluster_df) > 20:
            st.caption("Chart shows the first 20 ranked clusters; all candidates are in the table and selector below.")
    else:
        st.info("Full cluster metadata are unavailable in this database. Rebuild the taxon database to retain all member isolates, source categories, and computed types. Representative counts are not cluster sizes.")

    preferred = ["PDS_acc", "best_similarity_score", "near_top", "cluster_size", "env_cli_ratio",
                 "sampling_year_start", "sampling_year_end", "dated_isolates", "computed_serotypes",
                 "hits_in_results", "total_representatives", "SNP_tree_link"]
    st.dataframe(cluster_df[[c for c in preferred if c in cluster_df]], hide_index=True,
                 use_container_width=True, column_config={
                     "PDS_acc": "SNP cluster", "cluster_size": "Isolates in cluster",
                     "near_top": "Near top",
                     "sampling_year_start": st.column_config.NumberColumn("First year", format="%d"),
                     "sampling_year_end": st.column_config.NumberColumn("Last year", format="%d"),
                     "computed_serotypes": "Computed serotypes (counts)",
                     "hits_in_results": "Representative hits",
                     "total_representatives": "Representatives in database",
                     "env_cli_ratio": "Env : clinical", "dated_isolates": "Isolates with usable dates",
                     "best_similarity_score": st.column_config.NumberColumn("Best similarity", format="%.4f"),
                     "SNP_tree_link": st.column_config.LinkColumn("NCBI cluster", display_text="Open ↗"),
                 })
    if members.empty or cluster_df.empty:
        return
    st.subheader("Inside a cluster")
    selected = st.selectbox("Explore candidate cluster", cluster_df["PDS_acc"].tolist())
    frame = prepare_members(members[members["PDS_acc"].eq(selected)])
    if frame.empty:
        st.info("No full-member metadata are stored for this cluster.")
        return
    years = frame["Collection year"].dropna().astype(int)
    metrics = st.columns(4)
    metrics[0].metric("Isolates", f"{len(frame):,}")
    env = int(frame["Source group"].eq("Environmental / other").sum())
    clinical = int(frame["Source group"].eq("Clinical").sum())
    metrics[1].metric("Env : clinical", f"{env} : {clinical}")
    metrics[1].caption(f"{len(frame) - env - clinical:,} unknown / other")
    metrics[2].metric("Sampling years", f"{years.min()}–{years.max()}" if len(years) else "Unknown")
    metrics[2].caption(f"{len(years):,} / {len(frame):,} isolates with usable dates")
    known = int(frame["Computed serotype"].ne("Unknown").sum())
    metrics[3].metric("Computed serotype coverage", f"{known / len(frame):.0%}")
    metrics[3].caption(f"{known:,} / {len(frame):,} isolates")
    left, right = st.columns(2, gap="large")
    with left:
        metadata_chart(frame, "Computed serotype", "Computed serotype · all predictions retained")
        metadata_chart(frame, "Reported serovar", "Reported serovar · submitted metadata")
    with right:
        if len(years):
            annual = years.value_counts().rename_axis("Year").reset_index(name="Isolates")
            st.altair_chart(alt.Chart(annual).mark_bar(color="#e5a44e").encode(
                x=alt.X("Year:O", title="Collection year"), y="Isolates:Q",
                tooltip=["Year:O", "Isolates:Q"],
            ).properties(title="Sampling through time", height=280), use_container_width=True)
        else:
            st.info("No usable collection years are available.")
        metadata_chart(frame, "Country / region", "Geographic context")
    st.caption("Sampling years use valid ISO year, month, or date values. Date ranges and free text remain in the isolate table. Computed serotype and reported serovar are shown separately; conflicting predictions are retained.")
    with st.expander("Isolation sources and hosts"):
        left, right = st.columns(2)
        with left:
            metadata_chart(frame, "isolation_source", "Isolation source")
        with right:
            metadata_chart(frame, "host", "Host")
    with st.expander("All isolates in this cluster", expanded=False):
        fields = ["target_acc", "biosample_acc", "asm_acc", "strain", "Computed serotype",
                  "computed_types", "serovar", "epi_type", "collection_date", "geo_loc_name",
                  "isolation_source", "host"]
        st.dataframe(frame[fields], hide_index=True, use_container_width=True)
        st.download_button("Download cluster isolates", frame[fields].to_csv(index=False).encode(),
                           file_name=f"{selected}_isolates.csv", mime="text/csv")


@st.cache_data(show_spinner=False, max_entries=8)
def render_report_tree(newick, query_name, font_size, spacing):
    return render_tree(newick, query_name, font_size, spacing)


def display_results(results, db_summary):
    representative_df = results["representative_df"]
    cluster_df = results["cluster_df"]
    st.divider()
    st.subheader(f"Results · {results['query_name']}")
    st.caption(f"{db_summary.get('Type', 'Unknown')} database · local genome screening")
    columns = st.columns(3)
    columns[0].metric("Representative hits", len(representative_df))
    top_score = representative_df["similarity_score"].max() if "similarity_score" in representative_df else None
    columns[1].metric("Top similarity", f"{top_score:.4f}" if pd.notna(top_score) else "n/a")
    columns[2].metric("Candidate clusters", len(cluster_df) if cluster_df is not None else "n/a")
    overview, tree_tab, evidence = st.tabs(["Cluster context", "Sketch-distance tree", "Isolate evidence & downloads"])
    with overview:
        if cluster_df is not None:
            near_top_count = int(cluster_df["near_top"].sum())
            st.caption(f"{near_top_count} cluster(s) within the configured sketch-hash tolerance of the top hit. This is a screening heuristic, not a statistical confidence interval or outbreak confirmation.")
            display_cluster_context(cluster_df, results.get("members", pd.DataFrame()))
        else:
            st.caption("This database has no SNP cluster membership. Metadata below describe the returned isolates.")
            frame = prepare_members(representative_df)
            left, right = st.columns(2)
            with left:
                metadata_chart(frame, "Computed serotype", "Computed serotype")
                metadata_chart(frame, "isolation_source", "Isolation source")
            with right:
                metadata_chart(frame, "Country / region", "Geographic context")
                metadata_chart(frame, "Reported serovar", "Reported serovar")
    with tree_tab:
        st.caption("Neighbor-joining tree of sketch distances among the query and qualifying representatives. Branch distances are not SNP counts. The query is highlighted in red.")
        if results.get("tree_newick") is not None:
            controls = st.columns(2)
            font = controls[0].slider("Tip font size (pt)", 7, 18, 10)
            spacing = controls[1].slider("Tip spacing", 0.8, 2.0, 1.0, 0.1)
            with st.spinner("Rendering tree…"):
                png, svg, tips = render_report_tree(results["tree_newick"].decode(), results["query_name"], font, spacing)
            st.caption(f"{tips:,} tips · height follows tip count and font size. Large trees scroll vertically; SVG preserves detail at any zoom.")
            encoded = base64.b64encode(svg).decode("ascii")
            svg_height = re.search(r'<svg[^>]*height="([\d.]+)pt"', svg.decode())
            preview_height = min(650, max(220, int(float(svg_height.group(1)) * 4 / 3) + 32)) if svg_height else 650
            components.html(
                '<div style="background:white;padding:12px;width:max-content">'
                '<img alt="Candidate sketch-distance tree" style="max-width:none" '
                f'src="data:image/svg+xml;base64,{encoded}"></div>',
                height=preview_height, scrolling=True,
            )
            downloads = st.columns(3)
            downloads[0].download_button("Download PNG", png, file_name=f"{results['query_name']}_tree.png", mime="image/png")
            downloads[1].download_button("Download SVG", svg, file_name=f"{results['query_name']}_tree.svg", mime="image/svg+xml")
            downloads[2].download_button("Download Newick", results["tree_newick"], file_name=f"{results['query_name']}_tree.newick", mime="text/plain")
        elif results.get("tree_png") is not None:
            st.image(results["tree_png"], use_container_width=False)
        else:
            st.info(results.get("tree_skip_reason") or "No tree was generated.")
    with evidence:
        st.subheader("Representative-level matches")
        st.caption("Similarity scores belong to these representatives. Other members of their clusters were not individually compared with the query.")
        st.dataframe(representative_df, use_container_width=True, hide_index=True)
        st.download_button("Download representative matches", results["representative_csv"],
                           file_name=f"{results['query_name']}_representative_matches.csv", mime="text/csv")
        if cluster_df is not None:
            st.download_button("Download cluster summary", results["cluster_csv"],
                               file_name=f"{results['query_name']}_cluster_candidates.csv", mime="text/csv")
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
