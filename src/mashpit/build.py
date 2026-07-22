#!/usr/bin/env python3

import csv
import glob
import logging
import os
import re
import shutil
import sqlite3
import subprocess
import tarfile
import time
import urllib.request
import zipfile
from collections import defaultdict
from html.parser import HTMLParser
from pathlib import Path

import pandas as pd
import requests
import screed
from Bio import Entrez, Phylo
from sourmash import MinHash, SourmashSignature, save_signatures
from tqdm import tqdm


PDT_RE = re.compile(r"(PDT\d+)")
PDS_ARCHIVE_RE = re.compile(r"(PDS\d+\.\d+)\.tar\.gz$")


class LinkParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []

    def handle_starttag(self, tag, attrs):
        if tag != "a":
            return
        href = dict(attrs).get("href")
        if href:
            self.links.append(href)


class DownloadProgressBar(tqdm):
    def update_to(self, blocks=1, block_size=1, total_size=None):
        if total_size is not None:
            self.total = total_size
        self.update(blocks * block_size - self.n)


def create_connection(sql_path):
    return sqlite3.connect(sql_path)


def create_database(conn):
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS METADATA (
            biosample_acc TEXT PRIMARY KEY,
            taxid INTEGER,
            strain TEXT,
            collected_by TEXT,
            collection_date TEXT,
            geo_loc_name TEXT,
            isolation_source TEXT,
            lat_lon TEXT,
            serovar TEXT,
            sub_species TEXT,
            species TEXT,
            genus TEXT,
            host TEXT,
            host_disease TEXT,
            outbreak TEXT,
            srr TEXT,
            PDT_acc TEXT,
            PDS_acc TEXT,
            asm_acc TEXT
        )
        """
    )
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS DESC (
            name TEXT PRIMARY KEY,
            value TEXT
        )
        """
    )
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS REPRESENTATIVE (
            asm_acc TEXT PRIMARY KEY,
            PDT_acc TEXT,
            PDS_acc TEXT,
            tree_radius REAL,
            selection_round INTEGER,
            fasta_path TEXT,
            download_attempts INTEGER
        )
        """
    )
    conn.commit()


def download_url(url, output_path):
    output_path = str(output_path)
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    with DownloadProgressBar(
        unit="B",
        unit_scale=True,
        miniters=1,
        unit_divisor=1024,
        desc=Path(output_path).name,
        leave=False,
    ) as progress:
        urllib.request.urlretrieve(
            url,
            filename=output_path,
            reporthook=progress.update_to,
        )


def list_links(url):
    response = requests.get(url, timeout=120)
    response.raise_for_status()
    parser = LinkParser()
    parser.feed(response.text)
    return parser.links


def normalize_value(value):
    if pd.isna(value):
        return ""
    return str(value).strip()


def valid_assembly_accession(value):
    value = normalize_value(value)
    return value.startswith("GCA_") or value.startswith("GCF_")


def run_command(command, check=True):
    logging.debug("Running: %s", " ".join(str(item) for item in command))
    result = subprocess.run(
        [str(item) for item in command],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    if result.stdout:
        logging.debug(result.stdout)
    if result.stderr:
        logging.debug(result.stderr)
    if check and result.returncode != 0:
        raise RuntimeError(
            "Command failed with exit code %d: %s\n%s"
            % (
                result.returncode,
                " ".join(str(item) for item in command),
                result.stderr.strip(),
            )
        )
    return result


def prepare(args):
    if shutil.which("datasets") is None:
        raise RuntimeError("NCBI datasets was not found in PATH")

    timestamp = time.strftime("%Y%m%d%H%M%S")
    print_handler = logging.StreamHandler()
    file_handler = logging.FileHandler("mashpit-%s.log" % timestamp)
    print_handler.setLevel(logging.ERROR if args.quiet else logging.INFO)
    file_handler.setLevel(logging.DEBUG)

    logging.basicConfig(
        format="%(asctime)s [%(levelname)s] %(message)s",
        level=logging.DEBUG,
        datefmt="%d-%b-%y %H:%M:%S",
        handlers=[print_handler, file_handler],
    )

    db_folder = Path.cwd() / args.name
    if db_folder.exists():
        raise FileExistsError("Folder already exists: %s" % db_folder)
    db_folder.mkdir()

    tmp_folder = Path.cwd() / ("tmp_%s" % args.name)
    if tmp_folder.exists():
        shutil.rmtree(str(tmp_folder))
    tmp_folder.mkdir()

    conn = create_connection(db_folder / ("%s.db" % args.name))
    create_database(conn)
    return db_folder, tmp_folder, conn


def validate_pathogen_name(species):
    requested = species.strip().replace(" ", "_")
    base_url = "https://ftp.ncbi.nlm.nih.gov/pathogen/Results/"
    available = [item.strip("/") for item in list_links(base_url) if item.endswith("/")]

    for item in available:
        if item.lower() == requested.lower():
            return item

    raise ValueError(
        "Taxon %s is not available from NCBI Pathogen Detection" % requested
    )


def resolve_release(pathogen_name, pd_version):
    if pd_version:
        release_url = (
            "https://ftp.ncbi.nlm.nih.gov/pathogen/Results/%s/%s/"
            % (pathogen_name, pd_version)
        )
        response = requests.get(release_url, timeout=120)
        if response.status_code == 404:
            raise ValueError("Invalid Pathogen Detection release: %s" % pd_version)
        response.raise_for_status()
        return release_url, pd_version

    latest_url = (
        "https://ftp.ncbi.nlm.nih.gov/pathogen/Results/%s/latest_snps/"
        % pathogen_name
    )
    metadata_links = list_links(latest_url + "Metadata/")
    metadata_names = [
        Path(link).name
        for link in metadata_links
        if Path(link).name.endswith(".metadata.tsv")
    ]
    if not metadata_names:
        raise RuntimeError("No metadata file found under %s" % latest_url)

    pdg_acc = metadata_names[0].replace(".metadata.tsv", "")
    return latest_url, pdg_acc


def download_release_files(release_url, pdg_acc, tmp_folder):
    metadata_dir = tmp_folder / "metadata"
    tree_dir = tmp_folder / "trees"
    archive_dir = tmp_folder / "tree_archives"
    metadata_dir.mkdir()
    tree_dir.mkdir()
    archive_dir.mkdir()

    metadata_url = release_url + "Metadata/"
    cluster_url = release_url + "Clusters/"
    tree_url = release_url + "SNP_trees/"

    metadata_names = [
        Path(link).name
        for link in list_links(metadata_url)
        if Path(link).name.endswith(".metadata.tsv")
    ]
    if not metadata_names:
        raise RuntimeError("No metadata TSV found")

    metadata_name = metadata_names[0]
    isolate_name = "%s.reference_target.all_isolates.tsv" % pdg_acc

    metadata_path = metadata_dir / metadata_name
    isolate_path = metadata_dir / isolate_name

    download_url(metadata_url + metadata_name, metadata_path)
    download_url(cluster_url + isolate_name, isolate_path)

    archives = []
    for link in list_links(tree_url):
        name = Path(link).name
        match = PDS_ARCHIVE_RE.match(name)
        if match:
            archives.append((match.group(1), tree_url + name))

    if not archives:
        raise RuntimeError("No SNP-tree archives found under %s" % tree_url)

    logging.info("Downloading %d SNP-tree archives", len(archives))

    tree_paths = {}
    for index, (pds_acc, url) in enumerate(archives, start=1):
        archive_path = archive_dir / ("%s.tar.gz" % pds_acc)
        cluster_tree_dir = tree_dir / pds_acc
        cluster_tree_dir.mkdir()

        try:
            download_url(url, archive_path)
            with tarfile.open(str(archive_path), "r:gz") as archive:
                safe_extract_tar(archive, cluster_tree_dir)
            tree_path = find_newick_file(cluster_tree_dir)
            if tree_path:
                tree_paths[pds_acc] = tree_path
        except Exception as error:
            logging.error("Tree download failed for %s: %s", pds_acc, error)

        if index % 100 == 0 or index == len(archives):
            logging.info("Prepared %d/%d tree archives", index, len(archives))

    return metadata_path, isolate_path, tree_paths


def safe_extract_tar(archive, destination):
    destination = Path(destination).resolve()
    for member in archive.getmembers():
        target = (destination / member.name).resolve()
        if target != destination and not str(target).startswith(
            str(destination) + os.sep
        ):
            raise RuntimeError("Unsafe path in archive: %s" % member.name)
    archive.extractall(str(destination))


def find_newick_file(folder):
    candidates = []
    for pattern in ("*.newick", "*.nwk", "*.tree", "*.tre"):
        candidates.extend(Path(folder).rglob(pattern))
    if not candidates:
        return None
    return sorted(candidates, key=lambda path: path.stat().st_size, reverse=True)[0]


def load_metadata(metadata_path, isolate_path):
    metadata = pd.read_csv(metadata_path, sep="\t", dtype=str, low_memory=False)
    isolates = pd.read_csv(isolate_path, sep="\t", dtype=str, low_memory=False)

    if "target_acc" not in metadata.columns:
        raise ValueError("Metadata file is missing target_acc")
    if "asm_acc" not in metadata.columns:
        raise ValueError("Metadata file is missing asm_acc")
    if "target_acc" not in isolates.columns or "PDS_acc" not in isolates.columns:
        raise ValueError("Cluster file is missing target_acc or PDS_acc")

    metadata = metadata.merge(
        isolates[["target_acc", "PDS_acc"]].drop_duplicates("target_acc"),
        on="target_acc",
        how="left",
    )

    metadata["target_key"] = metadata["target_acc"].map(target_key)
    metadata["asm_acc"] = metadata["asm_acc"].map(normalize_value)
    metadata["PDS_acc"] = metadata["PDS_acc"].map(normalize_value)

    eligible = metadata[
        metadata["asm_acc"].map(valid_assembly_accession)
        & metadata["PDS_acc"].ne("")
        & metadata["target_key"].ne("")
    ].copy()

    eligible = eligible.drop_duplicates(
        subset=["PDS_acc", "target_key", "asm_acc"],
        keep="first",
    )
    return metadata, eligible


def target_key(value):
    match = PDT_RE.search(normalize_value(value))
    if match:
        return match.group(1)
    return normalize_value(value).split(".")[0]


def load_tree_graph(tree_path):
    tree = Phylo.read(str(tree_path), "newick")
    adjacency = defaultdict(list)
    terminal_by_key = {}

    for parent in tree.find_clades(order="level"):
        for child in parent.clades:
            length = child.branch_length
            if length is None:
                length = 0.0
            length = float(length)
            if length < 0:
                logging.warning(
                    "Negative branch length %.6g in %s; clamping to 0.0",
                    length,
                    tree_path,
                )
                length = 0.0
            adjacency[parent].append((child, length))
            adjacency[child].append((parent, length))

    for terminal in tree.get_terminals():
        key = target_key(terminal.name)
        if key:
            terminal_by_key[key] = terminal

    return tree, adjacency, terminal_by_key


def tree_distances(adjacency, start):
    distances = {start: 0.0}
    stack = [(start, None)]

    while stack:
        node, parent = stack.pop()
        base = distances[node]
        for neighbor, length in adjacency[node]:
            if neighbor is parent:
                continue
            distances[neighbor] = base + length
            stack.append((neighbor, node))

    return distances


def select_tree_representatives(tree_path, cluster_rows, radius, excluded):
    tree, adjacency, terminal_by_key = load_tree_graph(tree_path)

    row_by_key = {}
    node_by_key = {}

    for _, row in cluster_rows.sort_values("asm_acc").iterrows():
        asm_acc = row["asm_acc"]
        key = row["target_key"]
        if asm_acc in excluded or key not in terminal_by_key:
            continue
        if key in row_by_key:
            logging.debug(
                "Tree tip %s matched by multiple assemblies; keeping %s, dropping %s",
                key,
                row_by_key[key]["asm_acc"],
                asm_acc,
            )
            continue
        row_by_key[key] = row
        node_by_key[key] = terminal_by_key[key]

    keys = sorted(row_by_key)
    if not keys:
        return []

    if len(keys) == 1:
        return [row_by_key[keys[0]]]

    first_key = keys[0]
    first_distances = tree_distances(adjacency, node_by_key[first_key])
    endpoint_a = max(keys, key=lambda key: first_distances[node_by_key[key]])

    distances_a = tree_distances(adjacency, node_by_key[endpoint_a])
    endpoint_b = max(keys, key=lambda key: distances_a[node_by_key[key]])
    distances_b = tree_distances(adjacency, node_by_key[endpoint_b])

    center_key = min(
        keys,
        key=lambda key: (
            max(
                distances_a[node_by_key[key]],
                distances_b[node_by_key[key]],
            ),
            key,
        ),
    )

    selected = [center_key]
    center_distances = tree_distances(adjacency, node_by_key[center_key])
    nearest = {
        key: center_distances[node_by_key[key]]
        for key in keys
    }

    while True:
        farthest_key = max(keys, key=lambda key: (nearest[key], key))
        if nearest[farthest_key] <= radius:
            break

        selected.append(farthest_key)
        distances = tree_distances(adjacency, node_by_key[farthest_key])

        for key in keys:
            distance = distances[node_by_key[key]]
            if distance < nearest[key]:
                nearest[key] = distance

    return [row_by_key[key] for key in selected]


def select_all_representatives(
    eligible,
    tree_paths,
    radius,
    exclusions,
    round_number,
    pds_filter=None,
):
    representatives = []
    cluster_summary = []

    if pds_filter is not None:
        pds_filter = set(pds_filter)
        eligible = eligible[eligible["PDS_acc"].isin(pds_filter)].copy()

    grouped = eligible.groupby("PDS_acc", sort=True)
    total_clusters = len(grouped)

    for index, (pds_acc, cluster_rows) in enumerate(grouped, start=1):
        tree_path = tree_paths.get(pds_acc)
        if tree_path is None:
            cluster_summary.append(
                {
                    "PDS_acc": pds_acc,
                    "eligible_isolates": len(cluster_rows),
                    "representatives": 0,
                    "status": "tree_missing",
                }
            )
            continue

        try:
            selected = select_tree_representatives(
                tree_path,
                cluster_rows,
                radius,
                exclusions,
            )
            for row in selected:
                item = row.to_dict()
                item["selection_round"] = round_number
                representatives.append(item)

            cluster_summary.append(
                {
                    "PDS_acc": pds_acc,
                    "eligible_isolates": len(cluster_rows),
                    "representatives": len(selected),
                    "status": "complete" if selected else "no_eligible_tips",
                }
            )
        except Exception as error:
            logging.error("Representative selection failed for %s: %s", pds_acc, error)
            cluster_summary.append(
                {
                    "PDS_acc": pds_acc,
                    "eligible_isolates": len(cluster_rows),
                    "representatives": 0,
                    "status": "selection_error",
                }
            )

        if index % 100 == 0 or index == total_clusters:
            logging.info("Selected representatives for %d/%d clusters", index, total_clusters)

    return pd.DataFrame(representatives), pd.DataFrame(cluster_summary)


def chunks(values, size):
    for start in range(0, len(values), size):
        yield values[start : start + size]


def verify_fasta(assembly_root, accession):
    candidates = list(
        (assembly_root / "ncbi_dataset" / "data" / accession).glob("*_genomic.fna")
    )
    candidates = [path for path in candidates if path.is_file() and path.stat().st_size > 0]
    if len(candidates) == 1:
        return candidates[0]
    if len(candidates) > 1:
        return max(candidates, key=lambda path: path.stat().st_size)
    return None


def download_batch(accessions, output_dir, api_key=None):
    output_dir.mkdir(parents=True, exist_ok=True)
    accession_file = output_dir / "accessions.txt"
    package_file = output_dir / "package.zip"

    with accession_file.open("w", encoding="utf-8") as handle:
        for accession in accessions:
            handle.write(accession + "\n")

    command = [
        "datasets",
        "download",
        "genome",
        "accession",
        "--inputfile",
        str(accession_file),
        "--include",
        "genome",
        "--dehydrated",
        "--filename",
        str(package_file),
    ]
    if api_key:
        command.extend(["--api-key", api_key])

    run_command(command)

    with zipfile.ZipFile(str(package_file), "r") as archive:
        archive.extractall(str(output_dir))

    command = ["datasets", "rehydrate", "--directory", str(output_dir)]
    if api_key:
        command.extend(["--api-key", api_key])
    run_command(command)


def download_representatives(
    accessions,
    assembly_root,
    attempts,
    batch_size,
    retry_delay,
    api_key=None,
):
    assembly_root.mkdir(parents=True, exist_ok=True)
    verified = {}
    attempt_counts = defaultdict(int)
    errors = {}

    pending = []
    for accession in sorted(set(accessions)):
        existing = verify_fasta(assembly_root, accession)
        if existing:
            verified[accession] = existing
        else:
            pending.append(accession)

    for attempt in range(1, attempts + 1):
        if not pending:
            break

        logging.info(
            "Assembly download attempt %d/%d for %d accessions",
            attempt,
            attempts,
            len(pending),
        )

        next_pending = []

        for batch_number, batch in enumerate(chunks(pending, batch_size), start=1):
            batch_dir = assembly_root / (
                "attempt_%02d_batch_%05d" % (attempt, batch_number)
            )

            for accession in batch:
                attempt_counts[accession] += 1

            try:
                download_batch(batch, batch_dir, api_key)
            except Exception as error:
                logging.error("Batch download failed: %s", error)
                for accession in batch:
                    errors[accession] = str(error)

            for accession in batch:
                found = verify_fasta(batch_dir, accession)
                if found:
                    verified[accession] = found
                else:
                    next_pending.append(accession)
                    errors.setdefault(accession, "genomic FASTA not found")

        pending = sorted(set(next_pending))
        if pending and attempt < attempts:
            time.sleep(retry_delay)

    return verified, set(pending), dict(attempt_counts), errors


def sketch_assemblies(representatives, verified, signature_dir, hash_number, kmer_size):
    signature_dir.mkdir(parents=True, exist_ok=True)
    signature_paths = {}

    for index, accession in enumerate(sorted(verified), start=1):
        fasta_path = verified[accession]
        signature_path = signature_dir / ("%s.sig" % accession)

        if signature_path.exists() and signature_path.stat().st_size > 0:
            signature_paths[accession] = signature_path
            continue

        minhash = MinHash(n=hash_number, ksize=kmer_size)
        with screed.open(str(fasta_path)) as records:
            for record in records:
                minhash.add_sequence(record.sequence, force=True)

        signature = SourmashSignature(
            minhash,
            name=accession,
            filename=str(fasta_path),
        )
        with signature_path.open("wt", encoding="utf-8") as handle:
            save_signatures([signature], fp=handle)

        signature_paths[accession] = signature_path

        if index % 1000 == 0:
            logging.info("Sketched %d assemblies", index)

    return signature_paths


def merge_signatures(args, db_folder, signature_paths):
    signatures = []
    from sourmash import load_one_signature

    for accession in sorted(signature_paths):
        signatures.append(load_one_signature(str(signature_paths[accession])))

    output_path = db_folder / ("%s.sig" % args.name)
    with output_path.open("wt", encoding="utf-8") as handle:
        save_signatures(signatures, fp=handle)


def insert_metadata(conn, metadata, representatives, radius, verified, attempts):
    metadata_by_asm = {
        normalize_value(row["asm_acc"]): row
        for _, row in metadata.iterrows()
        if valid_assembly_accession(row.get("asm_acc", ""))
    }

    columns = [
        "biosample_acc",
        "taxid",
        "strain",
        "collected_by",
        "collection_date",
        "geo_loc_name",
        "isolation_source",
        "lat_lon",
        "serovar",
        "sub_species",
        "species",
        "genus",
        "host",
        "host_disease",
        "outbreak",
        "Run",
        "target_acc",
        "PDS_acc",
        "asm_acc",
    ]

    for _, representative in representatives.iterrows():
        accession = representative["asm_acc"]
        row = metadata_by_asm.get(accession)
        if row is None:
            continue

        values = []
        for column in columns:
            value = normalize_value(row.get(column, ""))
            values.append(value if value else "missing")

        conn.execute(
            """
            INSERT OR REPLACE INTO METADATA (
                biosample_acc, taxid, strain, collected_by, collection_date,
                geo_loc_name, isolation_source, lat_lon, serovar, sub_species,
                species, genus, host, host_disease, outbreak, srr, PDT_acc,
                PDS_acc, asm_acc
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            values,
        )

        conn.execute(
            """
            INSERT OR REPLACE INTO REPRESENTATIVE (
                asm_acc, PDT_acc, PDS_acc, tree_radius, selection_round,
                fasta_path, download_attempts
            ) VALUES (?, ?, ?, ?, ?, ?, ?)
            """,
            (
                accession,
                normalize_value(representative.get("target_acc", "")),
                normalize_value(representative.get("PDS_acc", "")),
                radius,
                int(representative.get("selection_round", 1)),
                str(verified[accession]),
                int(attempts.get(accession, 0)),
            ),
        )

    conn.commit()


def insert_accession_metadata(conn, accession_to_biosample, verified):
    columns = [
        "biosample_acc",
        "taxid",
        "strain",
        "collected_by",
        "collection_date",
        "geo_loc_name",
        "isolation_source",
        "lat_lon",
        "serovar",
        "sub_species",
        "species",
        "genus",
        "host",
        "host_disease",
        "outbreak",
        "srr",
        "PDT_acc",
        "PDS_acc",
        "asm_acc",
    ]

    for accession in sorted(verified):
        values = {column: "missing" for column in columns}
        values["asm_acc"] = accession
        values["biosample_acc"] = accession_to_biosample.get(accession, "missing")

        conn.execute(
            """
            INSERT OR REPLACE INTO METADATA (
                biosample_acc, taxid, strain, collected_by, collection_date,
                geo_loc_name, isolation_source, lat_lon, serovar, sub_species,
                species, genus, host, host_disease, outbreak, srr, PDT_acc,
                PDS_acc, asm_acc
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            [values[column] for column in columns],
        )

    conn.commit()


def write_build_tables(tmp_folder, representatives, cluster_summary, unavailable, errors):
    representatives.to_csv(
        tmp_folder / "representatives.tsv",
        sep="\t",
        index=False,
    )
    cluster_summary.to_csv(
        tmp_folder / "cluster_summary.tsv",
        sep="\t",
        index=False,
    )

    with (tmp_folder / "unavailable_assemblies.tsv").open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["assembly_accession", "error"])
        for accession in sorted(unavailable):
            writer.writerow([accession, errors.get(accession, "")])


def replace_cluster_rows(table, replacement, affected_clusters):
    if table is None or table.empty:
        kept = pd.DataFrame()
    else:
        kept = table[~table["PDS_acc"].isin(affected_clusters)].copy()

    if replacement is None or replacement.empty:
        return kept.reset_index(drop=True)

    return pd.concat([kept, replacement], ignore_index=True)


def build_taxon(args):
    db_folder, tmp_folder, conn = prepare(args)

    radius = float(getattr(args, "radius", 20.0))
    download_attempts = int(getattr(args, "download_attempts", 3))
    batch_size = int(getattr(args, "download_batch_size", 500))
    retry_delay = float(getattr(args, "retry_delay", 5.0))
    max_reselection_rounds = int(getattr(args, "max_reselection_rounds", 5))

    pathogen_name = validate_pathogen_name(args.species)
    release_url, pdg_acc = resolve_release(pathogen_name, args.pd_version)

    logging.info("Pathogen Detection release: %s", pdg_acc)
    logging.info("Tree coverage radius: %s", radius)

    metadata_path, isolate_path, tree_paths = download_release_files(
        release_url,
        pdg_acc,
        tmp_folder,
    )
    metadata, eligible = load_metadata(metadata_path, isolate_path)

    logging.info("Metadata rows: %d", len(metadata))
    logging.info("Assembly-backed rows: %d", len(eligible))
    logging.info("Trees available: %d", len(tree_paths))

    exclusions = set()
    verified = {}
    attempt_counts = defaultdict(int)
    all_errors = {}

    assembly_root = tmp_folder / "assemblies"

    representatives, cluster_summary = select_all_representatives(
        eligible,
        tree_paths,
        radius,
        exclusions,
        round_number=1,
    )

    if representatives.empty:
        raise RuntimeError("No representatives could be selected")

    reselection_round = 0
    while True:
        needed = sorted(set(representatives["asm_acc"]) - set(verified))

        newly_verified, failed, attempts, errors = download_representatives(
            needed,
            assembly_root,
            download_attempts,
            batch_size,
            retry_delay,
            getattr(args, "key", None),
        )

        verified.update(newly_verified)
        all_errors.update(errors)
        for accession, count in attempts.items():
            attempt_counts[accession] += count

        if not failed:
            final_representatives = representatives[
                representatives["asm_acc"].isin(verified)
            ].copy()
            final_summary = cluster_summary
            break

        affected_clusters = set(
            representatives.loc[
                representatives["asm_acc"].isin(failed),
                "PDS_acc",
            ]
        )

        logging.warning(
            "%d selected assemblies remained unavailable after retries; "
            "reselecting representatives for %d affected clusters",
            len(failed),
            len(affected_clusters),
        )

        exclusions.update(failed)

        if reselection_round == max_reselection_rounds:
            raise RuntimeError(
                "Representative selection did not stabilize after %d reselection rounds"
                % max_reselection_rounds
            )

        reselection_round += 1

        reselected, reselected_summary = select_all_representatives(
            eligible,
            tree_paths,
            radius,
            exclusions,
            round_number=reselection_round + 1,
            pds_filter=affected_clusters,
        )

        representatives = replace_cluster_rows(
            representatives,
            reselected,
            affected_clusters,
        )
        cluster_summary = replace_cluster_rows(
            cluster_summary,
            reselected_summary,
            affected_clusters,
        )

        if representatives.empty:
            raise RuntimeError("No representatives remain after reselection")

    missing_final = set(final_representatives["asm_acc"]) - set(verified)
    if missing_final:
        raise RuntimeError(
            "Final representatives lack FASTA files: %s"
            % ", ".join(sorted(missing_final)[:10])
        )

    signature_paths = sketch_assemblies(
        final_representatives,
        verified,
        tmp_folder / "signatures",
        args.number,
        args.ksize,
    )
    merge_signatures(args, db_folder, signature_paths)

    insert_metadata(
        conn,
        metadata,
        final_representatives,
        radius,
        verified,
        attempt_counts,
    )

    description = {
        "Type": "Taxonomy",
        "Species": pathogen_name,
        "Version": pdg_acc,
        "Hash_number": str(args.number),
        "Kmer_size": str(args.ksize),
        "Representative_method": "tree_radius",
        "Tree_radius": str(radius),
        "Representative_count": str(len(final_representatives)),
        "Unavailable_assembly_count": str(len(exclusions)),
    }
    conn.executemany(
        "INSERT OR REPLACE INTO DESC(name, value) VALUES (?, ?)",
        sorted(description.items()),
    )
    conn.commit()

    write_build_tables(
        tmp_folder,
        final_representatives,
        final_summary,
        exclusions,
        all_errors,
    )

    shutil.copy2(
        tmp_folder / "representatives.tsv",
        db_folder / "representatives.tsv",
    )
    shutil.copy2(
        tmp_folder / "cluster_summary.tsv",
        db_folder / "cluster_summary.tsv",
    )
    shutil.copy2(
        tmp_folder / "unavailable_assemblies.tsv",
        db_folder / "unavailable_assemblies.tsv",
    )

    conn.close()
    shutil.rmtree(str(tmp_folder))

    logging.info("Database complete: %s", db_folder)
    logging.info("Final representatives: %d", len(final_representatives))
    logging.info("Unavailable assemblies excluded: %d", len(exclusions))

def fetch_accession_metadata(biosample, email, api_key):
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    search = Entrez.esearch(db="assembly", term=biosample, retmax="1")
    record = Entrez.read(search)
    if not record["IdList"]:
        return None

    summary = Entrez.esummary(
        db="assembly",
        id=record["IdList"][0],
        report="full",
    )
    result = Entrez.read(summary)
    accession = result["DocumentSummarySet"]["DocumentSummary"][0][
        "AssemblyAccession"
    ]
    return str(accession)


def build_accession(args):
    db_folder, tmp_folder, conn = prepare(args)

    if not args.list or not Path(args.list).is_file():
        raise FileNotFoundError("Accession list not found: %s" % args.list)
    if not args.email:
        raise ValueError("Entrez email is required for accession builds")

    biosamples = [
        line.strip()
        for line in Path(args.list).read_text().splitlines()
        if line.strip()
    ]

    accessions = []
    accession_to_biosample = {}
    for biosample in biosamples:
        accession = fetch_accession_metadata(
            biosample,
            args.email,
            getattr(args, "key", None),
        )
        if accession:
            accessions.append(accession)
            accession_to_biosample[accession] = biosample

    verified, failed, attempts, errors = download_representatives(
        accessions,
        tmp_folder / "assemblies",
        int(getattr(args, "download_attempts", 3)),
        int(getattr(args, "download_batch_size", 500)),
        float(getattr(args, "retry_delay", 5.0)),
        getattr(args, "key", None),
    )

    if not verified:
        raise RuntimeError("No assemblies were downloaded")

    signature_paths = sketch_assemblies(
        pd.DataFrame({"asm_acc": sorted(verified)}),
        verified,
        tmp_folder / "signatures",
        args.number,
        args.ksize,
    )
    merge_signatures(args, db_folder, signature_paths)

    insert_accession_metadata(conn, accession_to_biosample, verified)

    conn.executemany(
        "INSERT OR REPLACE INTO DESC(name, value) VALUES (?, ?)",
        [
            ("Type", "Accession"),
            ("Hash_number", str(args.number)),
            ("Kmer_size", str(args.ksize)),
            ("Assembly_count", str(len(verified))),
            ("Unavailable_assembly_count", str(len(failed))),
        ],
    )
    conn.commit()
    conn.close()

    with (db_folder / "unavailable_assemblies.tsv").open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["assembly_accession", "error"])
        for accession in sorted(failed):
            writer.writerow([accession, errors.get(accession, "")])

    shutil.rmtree(str(tmp_folder))
    logging.info("Database complete: %s", db_folder)


def build(args):
    if args.type == "taxon":
        build_taxon(args)
    elif args.type == "accession":
        build_accession(args)
    else:
        raise ValueError("Unsupported database type: %s" % args.type)
