"""Standalone query reports using the same cluster summaries as the GUI."""
import base64
from datetime import datetime, timezone
from html import escape
from pathlib import Path

import pandas as pd


STYLE = """:root { color-scheme: light; font-family: system-ui, sans-serif; color: #203344; background: #f3f6f8; }
body { margin: 0; padding: 32px 20px; }
main { max-width: 1280px; margin: auto; }
h1 { font-size: 2.2rem; margin: 8px 0; overflow-wrap: anywhere; }
h2 { font-size: 1.35rem; margin-top: 0; }
p { line-height: 1.6; }
a { color: #126480; }
header { margin-bottom: 24px; }
.eyebrow { color: #126480; font-weight: 700; letter-spacing: .1em; text-transform: uppercase; }
.muted { color: #536878; }
.metrics { display: flex; flex-wrap: wrap; gap: 12px; margin: 24px 0; }
.metric, section { background: white; border: 1px solid #dce5ea; border-radius: 10px; padding: 22px; }
.metric { flex: 1; min-width: 170px; }
.metric strong { display: block; font-size: 1.6rem; margin-top: 8px; }
section { margin-bottom: 20px; }
.notice { border-left: 4px solid #c88c27; padding: 12px 16px; background: #fff7e7; }
.scroll { overflow: auto; max-height: 640px; }
table { border-collapse: collapse; width: 100%; font-size: .9rem; }
th, td { padding: 12px; border-bottom: 1px solid #e6edf1; text-align: left; vertical-align: top; }
th { background: #edf4f7; position: sticky; top: 0; white-space: nowrap; }
tr:nth-child(even) { background: #fafcfd; }
.downloads { display: flex; flex-wrap: wrap; gap: 10px; }
.downloads a { padding: 10px 15px; border: 1px solid #bdd2dd; border-radius: 6px; text-decoration: none; }
summary { cursor: pointer; font-weight: 600; padding: 8px 0; }
img { max-width: 100%; height: auto; }
code { overflow-wrap: anywhere; }
@media print { body { padding: 0; background: white; } .scroll { max-height: none; overflow: visible; }
  section { break-inside: avoid; } th { position: static; } .downloads { display: none; } }
"""


def table(frame):
    return '<div class="scroll">' + frame.to_html(
        index=False, border=0, escape=True, na_rep="Unavailable",
        float_format=lambda value: f"{value:.6f}",
    ) + "</div>"


def count(value):
    return "Unavailable" if pd.isna(value) else f"{int(value):,}"


def composition(row, column):
    total, value = row.get("cluster_size"), row.get(column)
    if pd.isna(total) or pd.isna(value) or total <= 0:
        return "Unavailable"
    return f"{int(value):,} ({100 * value / total:.1f}%)"


def download(path, label, mime="text/csv"):
    path = Path(path)
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return (f'<a download="{escape(path.name, quote=True)}" '
            f'href="data:{mime};base64,{encoded}">{escape(label)}</a>')


def render_query_report(query_name, database_name, settings, representatives,
                        clusters, tree_status, number_results, tie_tolerance_hashes):
    """Return self-contained HTML; CSVs and an available tree are embedded."""
    title = escape(query_name)
    release = escape(str(settings.get("Version", "Not recorded")))
    created = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    best = (f"{representatives['similarity_score'].max():.6f}"
            if not representatives.empty else "Unavailable")
    parts = [f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Mashpit report — {title}</title><style>{STYLE}</style></head><body><main>
<header><div class="eyebrow">Mashpit · query report</div><h1>{title}</h1>
<p class="muted">Database: <strong>{escape(database_name)}</strong> · Release: <strong>{release}</strong><br>
Generated {created}</p></header>
<div class="metrics"><div class="metric">Returned representatives<strong>{len(representatives):,}</strong></div>
<div class="metric">Best sketch similarity<strong>{best}</strong></div>''']
    if clusters is not None:
        parts.append(f'<div class="metric">Candidate SNP clusters<strong>{len(clusters):,}</strong></div>')
    parts.append('</div>')
    if clusters is not None:
        parts.append('<section><h2>SNP cluster context</h2>')
        near_top = int(clusters["near_top"].sum()) if "near_top" in clusters else 0
        if near_top > 1:
            parts.append(f'<p class="notice">{near_top} clusters are near the top within sketch resolution. '
                         'The screening result does not distinguish a unique cluster.</p>')
        parts.append('<p>Clusters are ranked by their best representative similarity. '
                     'Source percentages use all members of each cluster, including unknown sources. '
                     'These describe the database cluster, not the source of the query isolate.</p>')
        if clusters.get("cluster_size", pd.Series(dtype=float)).isna().any():
            parts.append('<p class="notice">Full cluster metadata are unavailable for some candidates. '
                         'Representative counts are not used as cluster sizes.</p>')
        overview = pd.DataFrame({
            "SNP cluster": clusters["PDS_acc"],
            "Best similarity": clusters["best_similarity_score"],
            "Near top": clusters["near_top"].map({True: "Yes", False: "No"}),
            "Cluster isolates": clusters["cluster_size"].map(count),
            "Clinical": clusters.apply(lambda row: composition(row, "clinical_count"), axis=1),
            "Environmental / other": clusters.apply(lambda row: composition(row, "environmental_count"), axis=1),
            "Unknown / other": clusters.apply(lambda row: composition(row, "unknown_source_count"), axis=1),
            "Representative hits": clusters["hits_in_results"],
        })
        parts.append(table(overview))
        context_columns = {
            "PDS_acc": "SNP cluster", "env_cli_ratio": "Environmental:clinical",
            "sampling_year_start": "First sampling year", "sampling_year_end": "Last sampling year",
            "dated_isolates": "Isolates with usable dates", "computed_serotypes": "Computed serotypes",
            "computed_serotype_known": "Isolates with computed serotype",
        }
        context = clusters[[column for column in context_columns if column in clusters]].copy()
        for column in ("sampling_year_start", "sampling_year_end", "dated_isolates", "computed_serotype_known"):
            if column in context:
                context[column] = context[column].map(lambda value: "Unavailable" if pd.isna(value) else str(int(value)))
        parts.append('<details><summary>Sampling dates and computed serotypes</summary>')
        parts.append(table(context.rename(columns=context_columns)))
        parts.append('<p class="muted">Date coverage counts valid ISO dates or years. Computed serotypes '
                     'are separate from submitted serovar metadata.</p></details></section>')
    else:
        parts.append('<section><h2>Isolate context</h2><p>This database has no SNP cluster membership. '
                     'The metadata below describe the returned isolates.</p></section>')

    parts.append('<section><h2>Local similarity tree</h2>')
    if tree_status.get("status") == "generated" and tree_status.get("png"):
        encoded = base64.b64encode(Path(tree_status["png"]).read_bytes()).decode("ascii")
        parts.append(f'<img alt="Similarity tree of the query and returned representatives" src="data:image/png;base64,{encoded}">')
    else:
        parts.append(f'<p>{escape(str(tree_status.get("reason") or "Tree unavailable."))}</p>')
    parts.append('<p class="muted">This tree uses sketch distances; it is not an SNP phylogeny.</p></section>')
    parts.append('<section><h2>Representative evidence</h2><details><summary>View all returned matches and metadata</summary>')
    parts.append(table(representatives))
    parts.append('</details></section><section><h2>Download results</h2><div class="downloads">')
    parts.append(download(query_name + "_representative_matches.csv", "Representative matches (CSV)"))
    if clusters is not None:
        parts.append(download(query_name + "_cluster_candidates.csv", "Cluster candidates (CSV)"))
    if tree_status.get("newick"):
        parts.append(download(tree_status["newick"], "Tree (Newick)", "text/plain"))
    parts.append('</div><p class="muted">Downloads and the displayed image are embedded. '
                 'This HTML file can be opened and shared on its own without a server or internet connection.</p></section>')
    parameters = escape(
        f"Database type: {settings.get('Type', 'Not recorded')}; "
        f"hashes: {settings.get('Hash_number', 'Not recorded')}; "
        f"k-mer size: {settings.get('Kmer_size', 'Not recorded')}; "
        f"requested representatives: {number_results}; near-top tolerance: {tie_tolerance_hashes} hashes."
    )
    parts.append(f'<footer><p class="muted">{parameters}</p><p>Similarity is screening evidence. '
                 'It does not establish outbreak membership or a confirmed SNP cluster assignment.</p></footer></main></body></html>')
    return "\n".join(parts)
