"""Cluster metadata summaries, independent of the GUI and sketch ranking."""
import re

import pandas as pd


MEMBER_COLUMNS = [
    "target_acc", "PDS_acc", "asm_acc", "biosample_acc", "strain",
    "collection_date", "epi_type", "geo_loc_name", "isolation_source",
    "host", "serovar", "computed_types",
]
MISSING = {"", "missing", "nan", "none", "null", "na", "n/a", "not available",
           "not collected", "not provided", "unknown"}


def clean(value):
    if pd.isna(value):
        return "Unknown"
    value = str(value).strip()
    return "Unknown" if value.lower() in MISSING else value


def source_group(value):
    value = clean(value).lower()
    if value in {"clinical", "clinical/human", "pathogen: clinical"}:
        return "Clinical"
    if value in {"environmental/other", "environmental", "environmental/food/other",
                 "pathogen: environmental/food/other"}:
        return "Environmental / other"
    return "Unknown / other"


def computed_serotype(value):
    """Keep distinct predictions and method labels; never silently pick a winner."""
    value = clean(value)
    if value == "Unknown":
        return value
    # Commas can be part of Salmonella antigenic formulae (4,[5],12:i:-).
    # Split on a delimiter only when it introduces another key/value pair.
    fields = re.split(r"[;,|]\s*(?=[A-Za-z][A-Za-z0-9_ ()/.-]*\s*[=:])", value)
    matches = []
    for field in fields:
        match = re.match(
            r"\s*((?:[^=;|]*[ _])?sero(?:type|var)(?:\s*\([^)]*\))?)\s*[=:]\s*(.+)",
            field, flags=re.I,
        )
        if match and clean(match.group(2)) != "Unknown":
            matches.append(f"{match.group(1).strip()}: {match.group(2).strip()}")
    return "; ".join(matches) or "Unknown"


def collection_year(value):
    # Only a single ISO year/month/date is suitable for an annual histogram.
    # Ranges and free text remain available in the raw metadata download.
    value = clean(value)
    if not re.fullmatch(r"\d{4}(?:-\d{2}(?:-\d{2})?)?", value):
        return None
    try:
        pd.Timestamp(value)
        return int(value[:4])
    except (ValueError, OverflowError):
        return None


def prepare_members(frame):
    frame = frame.copy()
    for column in MEMBER_COLUMNS:
        if column not in frame:
            frame[column] = "Unknown"
        frame[column] = frame[column].map(clean)
    frame["Source group"] = frame["epi_type"].map(source_group)
    frame["Computed serotype"] = frame["computed_types"].map(computed_serotype)
    frame["Reported serovar"] = frame["serovar"]
    frame["Collection year"] = frame["collection_date"].map(collection_year)
    frame["Country / region"] = frame["geo_loc_name"].map(lambda x: x.split(":", 1)[0])
    return frame


def store_cluster_members(conn, metadata):
    """Persist all source-release members, including isolates without assemblies."""
    frame = metadata.reindex(columns=MEMBER_COLUMNS).copy()
    frame = frame[frame["PDS_acc"].map(clean).ne("Unknown")]
    frame = frame[frame["target_acc"].map(clean).ne("Unknown")]
    # A target is an isolate; assembly and BioSample are not reliable counting keys.
    frame = frame.drop_duplicates(["PDS_acc", "target_acc"])
    frame.to_sql("CLUSTER_MEMBERS", conn, if_exists="replace", index=False, chunksize=10000)
    conn.execute("CREATE INDEX cluster_members_pds ON CLUSTER_MEMBERS(PDS_acc)")
    conn.commit()


def read_cluster_members(conn, clusters):
    exists = conn.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name='CLUSTER_MEMBERS'"
    ).fetchone()
    if not exists:
        return pd.DataFrame()
    clusters = list(dict.fromkeys(clusters))
    frames = []
    for start in range(0, len(clusters), 400):
        batch = clusters[start:start + 400]
        placeholders = ",".join("?" for _ in batch)
        frames.append(pd.read_sql_query(
            f"SELECT * FROM CLUSTER_MEMBERS WHERE PDS_acc IN ({placeholders})",
            conn, params=batch,
        ))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def summarize_clusters(members):
    rows = []
    if members.empty:
        return pd.DataFrame(columns=["PDS_acc", "cluster_size"])
    for cluster, frame in prepare_members(members).groupby("PDS_acc", sort=False):
        counts = frame["Source group"].value_counts()
        env = int(counts.get("Environmental / other", 0))
        clinical = int(counts.get("Clinical", 0))
        years = frame["Collection year"].dropna()
        types = frame["Computed serotype"]
        known_types = types[types.ne("Unknown")]
        rows.append({
            "PDS_acc": cluster, "cluster_size": len(frame),
            "environmental_count": env, "clinical_count": clinical,
            "unknown_source_count": len(frame) - env - clinical,
            "env_cli_ratio": f"{env}:{clinical}",
            "sampling_year_start": int(years.min()) if len(years) else None,
            "sampling_year_end": int(years.max()) if len(years) else None,
            "dated_isolates": len(years),
            "computed_serotypes": "; ".join(
                f"{name} ({count})" for name, count in known_types.value_counts().items()
            ) or "Unknown",
            "computed_serotype_known": len(known_types),
        })
    return pd.DataFrame(rows)
