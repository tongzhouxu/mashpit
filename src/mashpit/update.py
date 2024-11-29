import os
import time
import pandas as pd
import glob
import urllib.request
import logging
import shutil

from sourmash import load_file_as_signatures, save_signatures, load_one_signature
from mashpit.build import (
    create_connection,
    PDGParser,
    download_metadata,
    calculate_centroid,
    import_metadata,
    download_and_sketch_assembly,
)


def check_database_version(db_pathogen_name, db_version):
    # check if there is a new version
    url_metadata_folder = (
        "https://ftp.ncbi.nlm.nih.gov/pathogen/Results/"
        + db_pathogen_name
        + "/latest_snps/Metadata/"
    )
    response = urllib.request.urlopen(url_metadata_folder)
    metadata_html = response.read().decode("utf-8")
    parser = PDGParser()
    parser.feed(metadata_html)
    metadata_file_name = parser.list[0]
    latest_pdg_acc = metadata_file_name.strip(".metadata.tsv")
    if db_version == latest_pdg_acc:
        logging.error("%s is already the latest version" % db_version)
        exit(0)
    else:
        logging.info(
            "Updating %s from %s to %s" % (db_pathogen_name, db_version, latest_pdg_acc)
        )

    return latest_pdg_acc


def compare_tables(database_cluster_center, latest_cluster_center_metadata):
    logging.info("Comparing two databases")
    pds_to_remove = []
    asm_acc_remove = []
    pds_to_add = []
    asm_acc_add = []
    # 1. If PDS is no longer in latest metadata, remove old metadata record and signature
    pds_only_in_old = list(
        set(database_cluster_center["PDS_acc_no_suffix"].to_list())
        - set(latest_cluster_center_metadata["PDS_acc_no_suffix"].to_list())
    )
    if len(pds_only_in_old) > 0:
        for pds in pds_only_in_old:
            pds_to_remove.append(pds)
        for pds in database_cluster_center[
            database_cluster_center["PDS_acc_no_suffix"].isin(pds_to_remove)
        ]["asm_acc"].to_list():
            asm_acc_remove.append(pds)
    # 2. If PDS is present in both tables, check if the centroid that has changed
    pds_common = list(
        set(database_cluster_center["PDS_acc_no_suffix"].to_list()).intersection(
            latest_cluster_center_metadata["PDS_acc_no_suffix"].to_list()
        )
    )
    if len(pds_common) > 0:
        pds_changed = []
        for pds in pds_common:
            database_center = database_cluster_center[
                database_cluster_center["PDS_acc_no_suffix"] == pds
            ]["PDS_acc"].iloc[0]
            latest_center = latest_cluster_center_metadata[
                latest_cluster_center_metadata["PDS_acc_no_suffix"] == pds
            ]["PDS_acc"].iloc[0]
            if database_center != latest_center:
                pds_changed.append(pds)
        if len(pds_changed) > 0:
            for pds in pds_changed:
                pds_to_add.append(pds)
                pds_to_remove.append(pds)
            for pds in database_cluster_center[
                database_cluster_center["PDS_acc_no_suffix"].isin(pds_changed)
            ]["asm_acc"].to_list():
                asm_acc_remove.append(pds)
                asm_acc_add.append(pds)
    # 3. If PDS is only in the new table, add metadata and download signature files
    pds_new = list(
        set(latest_cluster_center_metadata["PDS_acc_no_suffix"].to_list())
        - set(database_cluster_center["PDS_acc_no_suffix"].to_list())
    )
    if len(pds_new) > 0:
        for pds in pds_new:
            pds_to_add.append(pds)
        for pds in latest_cluster_center_metadata[
            latest_cluster_center_metadata["PDS_acc_no_suffix"].isin(pds_new)
        ]["asm_acc"].to_list():
            asm_acc_add.append(pds)

    return pds_to_remove, asm_acc_remove, pds_to_add, asm_acc_add


def delete_metadata(conn, asm_acc):
    sql = "DELETE FROM METADATA WHERE asm_acc=?"
    cur = conn.cursor()
    try:
        cur.execute(sql, (asm_acc,))
    except:
        logging.error("Unable to delete %s" % asm_acc)
    cur.execute(sql, (asm_acc,))
    conn.commit()


def update_taxon(args, conn, tmp_folder, database_sig_path):
    c = conn.cursor()
    c.execute("SELECT value FROM DESC where name = 'Species';")
    db_pathogen_name = c.fetchone()[0]
    logging.info("Pathogen name: %s" % db_pathogen_name)
    c.execute("SELECT value FROM DESC where name = 'Version';")
    db_version = c.fetchone()[0]
    logging.info("Pathogen Detection version: %s" % db_version)
    c.execute("SELECT value FROM DESC where name = 'Hash_number';")
    hash_number = int(c.fetchone()[0])
    logging.info("Hash number: %s" % hash_number)
    c.execute("SELECT value FROM DESC where name = 'Kmer_size';")
    kmer_size = int(c.fetchone()[0])
    logging.info("Kmer size: %s" % kmer_size)

    latest_pdg_acc = check_database_version(db_pathogen_name, db_version)

    # get the cluster-centroid table from the original database
    database_cluster_center = pd.read_sql_query("select * from METADATA ", conn)
    # calculate the new cluster-centroid table for the new database
    metadata_file_name, isolate_pds_file = download_metadata(
        db_pathogen_name, latest_pdg_acc, tmp_folder
    )
    latest_metadata = pd.read_csv(
        os.path.join(tmp_folder, metadata_file_name), header=0, sep="\t"
    )
    df_isolate_pds = pd.read_csv(
        os.path.join(tmp_folder, isolate_pds_file), header=0, sep="\t"
    )
    # add cluster accession to metadata
    latest_metadata = latest_metadata.join(
        df_isolate_pds[["target_acc", "PDS_acc"]].set_index("target_acc"),
        on="target_acc",
    )
    latest_metadata_asm = latest_metadata[~latest_metadata["asm_acc"].isnull()]
    calculate_centroid(latest_metadata_asm, latest_pdg_acc, tmp_folder)
    latest_cluster_center = pd.read_csv(
        os.path.join(tmp_folder, f"{latest_pdg_acc}_cluster_center.tsv"), sep="\t"
    )
    latest_cluster_center_metadata = latest_cluster_center.join(
        latest_metadata_asm.drop("PDS_acc", axis=1).set_index("target_acc"),
        on="target_acc",
    )
    # add a new PDS accession column without suffix/PDS version
    database_cluster_center["PDS_acc_no_suffix"] = (
        database_cluster_center["PDS_acc"].str.split(".").str[0]
    )
    latest_cluster_center_metadata["PDS_acc_no_suffix"] = (
        latest_cluster_center_metadata["PDS_acc"].str.split(".").str[0]
    )
    pds_to_remove, asm_acc_remove, pds_to_add, asm_acc_add = compare_tables(
        database_cluster_center, latest_cluster_center_metadata
    )
    logging.info("Found %s PDS to remove" % len(pds_to_remove))
    logging.info("Found %s PDS to add" % len(pds_to_add))

    if len(pds_to_remove) > 0:
        for asm_acc in asm_acc_remove:
            delete_metadata(conn, asm_acc)
    if len(pds_to_add) > 0:
        import_metadata(
            latest_metadata,
            latest_cluster_center_metadata[
                latest_cluster_center_metadata["PDS_acc_no_suffix"].isin(pds_to_add)
            ],
            conn,
        )
        download_and_sketch_assembly(asm_acc_add, hash_number, kmer_size, tmp_folder)

    database_sig = load_file_as_signatures(database_sig_path)
    siglist = []
    for sig in database_sig:
        if str(sig) not in asm_acc_remove:
            siglist.append(sig)
    new_sig_path = os.path.join(tmp_folder, "signature", "*.sig")
    sig_path_list = glob.glob(new_sig_path)
    for sig in sig_path_list:
        siglist.append(load_one_signature(sig))
    with open(database_sig_path, "w") as f:
        save_signatures(siglist, fp=f)
    shutil.rmtree(tmp_folder)

    # update the PDS version info in DESC table in the db
    c.execute(f"Update DESC set value = '{latest_pdg_acc}' where name = 'Version';")

    conn.commit()
    conn.close()
    logging.info(
        f"Database {args.name} updated successfully to version {latest_pdg_acc}"
    )


def update_accession(args, conn, tmp_folder):
    # update the metadata or add new columns
    metadata_file = args.metadata
    if not os.path.exists(metadata_file):
        logging.error("Metadata file not found.")
        exit(1)
    # check if the metadata file has the biosample_acc columns
    metadata = pd.read_csv(metadata_file, header=0)
    if "biosample_acc" not in metadata.columns:
        logging.error("Metadata file does not have biosample_acc column.")
        exit(1)
    # check if other columns are present in the database
    c = conn.cursor()
    c.execute("SELECT * FROM METADATA")
    columns = [description[0] for description in c.description]
    for col in metadata.columns:
        if col not in columns:
            # add the column to the database
            c.execute(f"ALTER TABLE METADATA ADD COLUMN '{col}'")
    conn.commit()
    # import the metadata
    import_metadata(metadata, conn)


def update(args):
    t = time.localtime()
    current_time = time.strftime("%Y%m%d%H%M%S", t)
    print_handler = logging.StreamHandler()
    file_handler = logging.FileHandler("mashpit-" + current_time + ".log")
    print_handler.setLevel(logging.INFO)
    file_handler.setLevel(logging.DEBUG)
    if args.quiet:
        print_handler.setLevel(logging.ERROR)
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s] %(message)s",
        level=logging.DEBUG,
        datefmt="%d-%b-%y %H:%M:%S",
        handlers=[print_handler, file_handler],
    )

    db_folder = os.path.abspath(args.database)
    # create a tmp folder
    cwd = os.getcwd()
    tmp_folder = os.path.join(cwd, "tmp")
    if os.path.isdir(tmp_folder):
        shutil.rmtree(tmp_folder)
    os.mkdir(tmp_folder)
    if not os.path.exists(db_folder):
        logging.error("Path not found.")
        exit(1)
    sql_path = os.path.join(db_folder, args.name + ".db")
    sig_path = os.path.join(db_folder, args.name + ".sig")
    if not (os.path.exists(sql_path)):
        logging.error(f"{sql_path} not found.")
        exit(1)
    if not (os.path.exists(sig_path)):
        logging.error(f"{sig_path} not found.")
        exit(1)

    conn = create_connection(sql_path)
    logging.info("Updating database %s" % args.name)
    c = conn.cursor()
    c.execute("SELECT value FROM DESC where name = 'Type';")
    db_type = c.fetchone()[0]
    logging.info("Database type: %s" % db_type)
    if db_type == "Taxonomy":
        update_taxon(args, conn, tmp_folder, sig_path)
    elif db_type == "Accession":
        update_accession(args, conn, tmp_folder)
