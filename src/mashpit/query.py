#!/usr/bin/env python3
import os
import glob
import logging
import ntpath
import screed
import sourmash
import heapq
import tempfile
import time
import pandas as pd

from skbio import DistanceMatrix
from skbio.tree import nj
from phytreeviz import TreeViz
from mashpit.build import create_connection
from operator import itemgetter
from sourmash import (
    SourmashSignature,
    save_signatures,
    load_one_signature,
    load_file_as_signatures,
)


class MashtreeSkipped(Exception):
    pass


def get_query_sig(query_path, query_name, hash_number, kmer_size, output_dir):
    mh = sourmash.MinHash(n=hash_number, ksize=kmer_size)
    for record in screed.open(query_path):
        mh.add_sequence(record.sequence, True)
    sig = SourmashSignature(mh, name=query_path)
    query_sig_path = os.path.join(output_dir, f"{query_name}.sig")
    with open(query_sig_path, "wt") as fp:
        save_signatures([sig], fp)
    return query_sig_path


def generate_query_table(conn, sorted_asm_similarity_dict):
    c = conn.cursor()
    c.execute("SELECT * FROM METADATA")
    names = [description[0] for description in c.description]

    # List to collect DataFrames
    df_list = []

    for i in sorted_asm_similarity_dict:
        sql_query = pd.read_sql_query(
            "SELECT * FROM METADATA WHERE asm_acc = ?", conn, params=(i,)
        )
        df_query = pd.DataFrame(sql_query, columns=names)
        df_query["similarity_score"] = sorted_asm_similarity_dict[i]
        df_list.append(df_query)

    # Concatenate all DataFrames at once
    output_df = pd.concat(df_list, ignore_index=True)

    # Additional processing (if it is a taxon database, add the link of the SNP cluster to the output)
    c.execute("SELECT value FROM DESC WHERE name = 'Type';")
    db_type = c.fetchone()[0]
    if db_type == "Taxonomy":
        pds_list = output_df["PDS_acc"].to_list()
        cluster_link = [
            f"https://www.ncbi.nlm.nih.gov/pathogens/isolates/#{pds}"
            for pds in pds_list
        ]
        output_df["SNP_tree_link"] = cluster_link

    return output_df


def generate_cluster_table(conn, representative_df, hash_number, tie_tolerance_hashes=2):
    # Groups the per-representative hits by SNP cluster (PDS_acc) so that
    # multiple representatives supporting the same cluster show up as one
    # row instead of scattered duplicates. hits_in_results/
    # total_representatives let a reader judge how much of that evidence is
    # real (similarity-backed) versus how large the cluster's representative
    # set happens to be; ranking itself always stays on best_similarity_score
    # - representative count is context, never a substitute ranking key.
    #
    # near_top is defined relative to the sketch's own resolution
    # (1 / hash_number) rather than a fixed similarity cutoff: at the default
    # hash_number=1000, a 0.001 gap is a single hash of difference and within
    # normal MinHash sampling noise, so a fixed threshold like 0.85 can't
    # tell "genuinely behind" apart from "statistically tied for best".
    grouped = (
        representative_df.groupby("PDS_acc")["similarity_score"]
        .agg(
            hits_in_results="count",
            best_similarity_score="max",
            mean_similarity_score="mean",
            min_similarity_score="min",
        )
        .reset_index()
    )

    total_counts = pd.read_sql_query(
        "SELECT PDS_acc, COUNT(*) AS total_representatives "
        "FROM REPRESENTATIVE GROUP BY PDS_acc",
        conn,
    )
    cluster_df = grouped.merge(total_counts, on="PDS_acc", how="left")

    overall_best_score = representative_df["similarity_score"].max()
    tolerance = tie_tolerance_hashes / hash_number
    cluster_df["similarity_gap_from_best"] = (
        cluster_df["best_similarity_score"] - overall_best_score
    )
    cluster_df["near_top"] = cluster_df["similarity_gap_from_best"] >= -tolerance

    cluster_df["SNP_tree_link"] = cluster_df["PDS_acc"].map(
        lambda pds: f"https://www.ncbi.nlm.nih.gov/pathogens/isolates/#{pds}"
    )

    cluster_df = cluster_df.sort_values(
        by=["best_similarity_score", "hits_in_results", "PDS_acc"],
        ascending=[False, False, True],
    ).reset_index(drop=True)

    return cluster_df


def generate_mashtree(
    output_df, min_similarity, query_name, sig_path, added_annotation, database_sig
):
    # check if the top query similarity is smaller than the threshold
    if float(output_df["similarity_score"].iloc[0]) < min_similarity:
        raise MashtreeSkipped(
            "Top query similarity is smaller than the threshold. Mashtree can not be generated. "
        )
    # select top results that are above the threshold
    acc_list = output_df[output_df["similarity_score"] >= min_similarity][
        "asm_acc"
    ].to_list()
    # if the number of top results is smaller than 2, mashtree can not be generated
    if len(acc_list) < 2:
        raise MashtreeSkipped(
            "Number of top results is smaller than 2. Mashtree can not be generated. "
        )
    leaves = [query_name] + acc_list
    # select database signatures that are present in acc_list
    sigs = []
    for sig in database_sig:
        if sig.name in acc_list:
            sigs.append(sig)
    # sort sigs based on the order of acc_list
    sigs = sorted(sigs, key=lambda x: acc_list.index(x.name))
    matrix = []
    # first row
    first_row = [0]
    for similarity_score in output_df[output_df["similarity_score"] >= min_similarity][
        "similarity_score"
    ].to_list():
        first_row.append(1 - similarity_score)
    matrix.append(first_row)
    for acc in acc_list:
        distance_list = []
        distance_list.append(
            1 - (output_df[output_df["asm_acc"] == acc]["similarity_score"].iloc[0])
        )
        for sig in sigs:
            if sig.name == acc:
                query_sig = sig
        for sig in sigs:
            distance_list.append(1 - (query_sig.jaccard(sig)))
        matrix.append(distance_list)

    dm = DistanceMatrix(matrix, leaves)
    newick_str = nj(dm, result_constructor=str)
    with open(f"{query_name}_tree.newick", "w") as f:
        f.write(newick_str)

    # add annotation
    if added_annotation is not None:
        annotated_leaves = []
        for leaf in leaves:
            if leaf in acc_list:
                annotation_value = str(
                    output_df[output_df["asm_acc"] == leaf][added_annotation].iloc[0]
                )
                leaf = leaf + " " + annotation_value
            annotated_leaves.append(leaf)
        dm = DistanceMatrix(matrix, annotated_leaves)
        newick_str = nj(dm, result_constructor=str)
        with open(f"{query_name}_tree.newick", "w") as f:
            f.write(newick_str)
    tv = TreeViz(f"{query_name}_tree.newick")
    tv.set_node_label_props(query_name, color="red")
    tv.savefig(f"{query_name}_tree.png", dpi=300)


def query(args):
    t = time.localtime()
    current_time = time.strftime("%Y%m%d%H%M%S", t)
    print_handler = logging.StreamHandler()
    file_handler = logging.FileHandler(f"mashpit-{current_time}.log")
    print_handler.setLevel(logging.INFO)
    file_handler.setLevel(logging.DEBUG)
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s] %(message)s",
        level=logging.DEBUG,
        datefmt="%d-%b-%y %H:%M:%S",
        handlers=[print_handler, file_handler],
    )
    logging.info("Start querying")
    time_start = time.time()
    number_results = args.number
    min_similarity = args.threshold
    added_annotation = args.annotation
    query_path = os.path.abspath(args.sample)
    query_name = ntpath.basename(query_path)

    # remove the extension if any
    if "." in query_name:
        query_name = query_name.split(".")[0]
    db_folder = os.path.abspath(args.database)
    if not os.path.exists(db_folder):
        logging.error("Database path not found.")
        exit(1)
    # find files ending with .db
    sql_path = glob.glob(os.path.join(db_folder, "*.db"))
    if len(sql_path) == 0:
        logging.error("Database path not found.")
        exit(1)
    if len(sql_path) > 1:
        logging.error("Multiple database files found.")
        exit(1)
    sql_path = sql_path[0]
    # find files ending with .sig
    sig_path = glob.glob(os.path.join(db_folder, "*.sig"))
    if len(sig_path) == 0:
        logging.error("Signature file not found.")
        exit(1)
    if len(sig_path) > 1:
        logging.error("Multiple signature files found.")
        exit(1)
    sig_path = sig_path[0]
    # check if the basename of .db and .sig are the same
    if (
        os.path.basename(sql_path).split(".")[0]
        != os.path.basename(sig_path).split(".")[0]
    ):
        logging.error("Database sql and signature files do not match.")
        exit(1)

    # check the hash number and kmer size in the database
    conn = create_connection(sql_path)
    c = conn.cursor()
    c.execute("SELECT value FROM DESC where name = 'Hash_number';")
    hash_number = int(c.fetchone()[0])
    c.execute("SELECT value FROM DESC where name = 'Kmer_size';")
    kmer_size = int(c.fetchone()[0])

    # Sketch the query sample and load the signature. Written to a private
    # temp directory rather than beside the input assembly, so a read-only
    # or shared input location doesn't break the query and no stray file of
    # the same name there is ever touched.
    with tempfile.TemporaryDirectory(prefix="mashpit-query-sketch-") as sketch_dir:
        query_sig_path = get_query_sig(
            query_path, query_name, hash_number, kmer_size, sketch_dir
        )
        query_sig = load_one_signature(query_sig_path)

    time_finish_sketch = time.time()
    logging.info(
        f"Query sample sketched in {time_finish_sketch-time_start:.2f} seconds"
    )

    asm_similarity_dict = {}
    database_sig = list(load_file_as_signatures(sig_path))

    time_load_database_sig = time.time()
    logging.info(
        f"Database signatures loaded in {time_load_database_sig-time_finish_sketch:.2f} seconds"
    )

    for sig in database_sig:
        similarity = query_sig.jaccard(sig)
        asm_similarity_dict[str(sig)] = similarity

    time_calculate_similarity = time.time()
    logging.info(
        f"Jaccard similarity calculated in {time_calculate_similarity-time_load_database_sig:.2f} seconds"
    )

    # get the top results
    top_items = heapq.nlargest(
        number_results, asm_similarity_dict.items(), key=itemgetter(1)
    )
    sorted_asm_similarity_dict = dict(top_items)

    time_sort = time.time()
    logging.info(
        f"Top {number_results} results sorted in {time_sort-time_calculate_similarity:.2f} seconds"
    )

    output_df = generate_query_table(conn, sorted_asm_similarity_dict)
    output_df.to_csv(query_name + "_representative_matches.csv", index=True)

    c.execute("SELECT value FROM DESC WHERE name = 'Type';")
    db_type = c.fetchone()[0]
    if db_type == "Taxonomy":
        cluster_df = generate_cluster_table(
            conn, output_df, hash_number, args.tie_tolerance_hashes
        )
        cluster_df.to_csv(query_name + "_cluster_candidates.csv", index=True)

    # A skipped tree (top hit below threshold, or fewer than two qualifying
    # hits) is a legitimate outcome, not a failure: the CSVs above are
    # already valid and complete, so it must not turn into a non-zero exit
    # for the whole query.
    try:
        generate_mashtree(
            output_df,
            min_similarity,
            query_name,
            sig_path,
            added_annotation,
            database_sig,
        )
        time_mashtree = time.time()
        logging.info(f"Mashtree generated in {time_mashtree-time_sort:.2f} seconds")
    except MashtreeSkipped as error:
        logging.warning(str(error))
