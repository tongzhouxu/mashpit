#!/usr/bin/env python3
import os
import time
import glob
import csv
import itertools
import zipfile
import screed
import shutil
import logging
import sqlite3
import requests
import psutil
import subprocess

import urllib.request
import dask.dataframe as dd
import pandas as pd
import xml.etree.ElementTree as ET

from tqdm import tqdm
from Bio import Entrez
from html.parser import HTMLParser
from sourmash import SourmashSignature, save_signatures, load_one_signature, MinHash


def log_memory_usage(stage):
    process = psutil.Process()
    mem_info = process.memory_info()
    logging.info(
        f"{stage} - RSS: {mem_info.rss / 1024 ** 2:.2f} MB, VMS: {mem_info.vms / 1024 ** 2:.2f} MB"
    )


# create connection to sqlite db file
def create_connection(sql_path):
    conn = sqlite3.connect(sql_path)
    return conn


# create sqlite3 database with the tables
def create_database(conn):
    # metadata and accessions of a sample
    metadata_table = """CREATE TABLE IF NOT EXISTS METADATA (
                              biosample_acc    TEXT PRIMARY KEY, 
                              taxid            INTEGER,
                              strain           TEXT, 
                              collected_by     TEXT,
                              collection_date  TEXT,
                              geo_loc_name     TEXT,
                              isolation_source TEXT,
                              lat_lon          TEXT,
                              serovar          TEXT,
                              sub_species      TEXT,
                              species          TEXT,
                              genus            TEXT,
                              host             TEXT,
                              host_disease     TEXT,
                              outbreak         TEXT,
                              srr              TEXT,
                              PDT_acc          TEXT,
                              PDS_acc          TEXT,
                              asm_acc          TEXT
                        );"""
    # a description table to keep a record of the dabtabase type, (species and PD version if it is a standard database)
    description_table = """CREATE TABLE IF NOT EXISTS DESC (
                              name    TEXT PRIMARY KEY, 
                              value            TEXT
                        );"""
    # create tables
    c = conn.cursor()
    c.execute(metadata_table)
    c.execute(description_table)
    conn.commit()


class PDGParser(HTMLParser):
    list = []

    def handle_data(self, data):
        if str(data).endswith("tsv"):
            self.list.append(data)


class Species_Name_Parser(HTMLParser):
    list = []

    def handle_data(self, data):
        if str(data).endswith("/"):
            self.list.append(data.strip("/"))


class DownloadProgressBar(tqdm):
    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)


def download_url(url, output_path):
    with DownloadProgressBar(
        unit="B",
        unit_scale=True,
        miniters=1,
        unit_divisor=1024,
        desc=url.split("/")[-1],
        leave=False,
    ) as t:
        urllib.request.urlretrieve(url, filename=output_path, reporthook=t.update_to)


def download_metadata(pathogen_name, pd_version, tmp_folder):
    logging.info("Downloading metadata files...")
    time_start_download = time.time()
    # get the metadata and SNP distance files
    if pd_version is not None:
        r = requests.get(
            f"https://ftp.ncbi.nlm.nih.gov/pathogen/Results/{pathogen_name}/{pd_version}/Metadata/"
        )
        if r.status_code == 404:
            logging.error("invalid PDG accession")
            exit(0)
        elif r.status_code == 200:
            url_metadata_folder = f"https://ftp.ncbi.nlm.nih.gov/pathogen/Results/{pathogen_name}/{pd_version}/Metadata/"
    else:
        url_metadata_folder = f"https://ftp.ncbi.nlm.nih.gov/pathogen/Results/{pathogen_name}/latest_snps/Metadata/"
    response = urllib.request.urlopen(url_metadata_folder)
    metadata_html = response.read().decode("utf-8")
    parser = PDGParser()
    parser.feed(metadata_html)
    metadata_file_name = parser.list[0]
    pdg_acc = metadata_file_name.strip(".metadata.tsv")
    distance_file = f"{pdg_acc}.reference_target.SNP_distances.tsv"
    isolate_pds_file = f"{pdg_acc}.reference_target.all_isolates.tsv"
    parser.close()
    url_metadata = url_metadata_folder + metadata_file_name
    if pd_version is not None:
        url_snp_distances = f"https://ftp.ncbi.nlm.nih.gov/pathogen/Results/{pathogen_name}/{pd_version}/Clusters/{distance_file}"
        url_isolate_pds = f"https://ftp.ncbi.nlm.nih.gov/pathogen/Results/{pathogen_name}/{pd_version}/Clusters/{isolate_pds_file}"
    else:
        url_snp_distances = f"https://ftp.ncbi.nlm.nih.gov/pathogen/Results/{pathogen_name}/latest_snps/Clusters/{distance_file}"
        url_isolate_pds = f"https://ftp.ncbi.nlm.nih.gov/pathogen/Results/{pathogen_name}/latest_snps/Clusters/{isolate_pds_file}"
    # check if the metadata file exists
    if not os.path.isfile(os.path.join(tmp_folder, metadata_file_name)):
        logging.info(f"Downloading metadata file from {url_metadata_folder}")
        download_url(url_metadata, os.path.join(tmp_folder, metadata_file_name))
    else:
        logging.info(f"Using existing metadata file {metadata_file_name}")
    # check if the SNP distance file exists
    if not os.path.isfile(os.path.join(tmp_folder, distance_file)):
        logging.info(f"Downloading pairwise SNP distance file from {url_snp_distances}")
        download_url(url_snp_distances, os.path.join(tmp_folder, distance_file))
    else:
        logging.info(f"Using existing pairwise SNP distance file {distance_file}")
    # check if the isolate PDS acc file exists
    if not os.path.isfile(os.path.join(tmp_folder, isolate_pds_file)):
        logging.info(f"Downloading isolates_PDS_acc file from {url_isolate_pds}")
        download_url(url_isolate_pds, os.path.join(tmp_folder, isolate_pds_file))
    else:
        logging.info(f"Using existing isolates_PDS_acc file {isolate_pds_file}")

    time_end_download = time.time()
    logging.info(
        f"Downloaded metadata and distance files in {round(time_end_download-time_start_download,2)} seconds."
    )
    # return the metadata file name
    return metadata_file_name, isolate_pds_file


def calculate_centroid(df_metadata, pdg_acc, tmp_folder):
    logging.info("Calculating centroid...")
    time_start_calculate_centroid = time.time()
    # drop columns to reduce memory usage
    input_file_name = f"{pdg_acc}.reference_target.SNP_distances.tsv"
    output_file_name = f"{pdg_acc}_distance.tsv"
    input_file = os.path.join(tmp_folder, input_file_name)
    output_file = os.path.join(tmp_folder, output_file_name)
    columns_to_select = [0, 4, 8, 11]

    with open(input_file, "r") as f_in, open(output_file, "w") as f_out:
        reader = csv.reader(f_in, delimiter="\t")
        writer = csv.writer(f_out, delimiter="\t")
        for row in reader:
            selected_cols = [row[i] for i in columns_to_select]
            writer.writerow(selected_cols)

    # filter out isolates without asm_acc
    target_columns = ["target_acc_1", "target_acc_2"]
    distance_filtered_name = f"{pdg_acc}_filtered_distance.tsv"
    distance_filtered = os.path.join(tmp_folder, distance_filtered_name)
    # Use the csv module to open the TSV file and create a reader object
    with open(output_file, "r") as distance_file, open(
        distance_filtered, "w"
    ) as filtered_distance_file:

        reader = csv.reader(distance_file, delimiter="\t")
        writer = csv.writer(filtered_distance_file, delimiter="\t")

        # Extract the header row from the reader object and write it to the new file
        header = next(reader)
        writer.writerow(header)

        # Find the indices of the target columns in the header row
        target_column_indices = [header.index(column) for column in target_columns]

        # Iterate over each row in the reader object and write it to the new file if it contains the target columns
        target_values = set(df_metadata["target_acc"])
        for row in reader:
            if all(row[index] in target_values for index in target_column_indices):
                writer.writerow(row)

    # Define the column to group by
    group_by_column = "PDS_acc"
    # two lists to be added into the database
    cluster_list_db = []
    centroid_list_db = []
    # Open the TSV file using the csv module and read it line by line
    with open(distance_filtered, "r") as file:
        reader = csv.reader(file, delimiter="\t")
        header = next(reader)  # read the header row

        # Define a generator expression that yields each row from the file
        rows = (row for row in reader)

        # Use itertools.groupby() to group the rows by the specified column
        for group_value, group_rows in itertools.groupby(
            rows, lambda x: x[header.index(group_by_column)]
        ):
            # Process the rows for the current group
            cluster_list_db.append(group_value)
            all_rows = list(group_rows)
            col1 = [row[0] for row in all_rows]
            col2 = [row[1] for row in all_rows]
            col3 = [row[3] for row in all_rows]
            df_cluster = pd.DataFrame(
                {
                    "target_acc_1": col1,
                    "target_acc_2": col2,
                    "delta_positions_unambiguous": col3,
                }
            )
            df_cluster = pd.concat(
                [
                    df_cluster,
                    df_cluster.rename(
                        columns={
                            "target_acc_1": "target_acc_2",
                            "target_acc_2": "target_acc_1",
                        }
                    ),
                ]
            )
            df_cluster["delta_positions_unambiguous"] = df_cluster[
                "delta_positions_unambiguous"
            ].astype("float")
            df_dask = dd.from_pandas(
                df_cluster, npartitions=4
            )  # Convert the pandas dataframe to a Dask dataframe with 4 partitions
            centroid_target_acc = (
                df_dask.groupby("target_acc_1")["delta_positions_unambiguous"]
                .sum()
                .compute()
                .idxmin()
            )
            centroid_list_db.append(centroid_target_acc)

    df_cluster_center = pd.DataFrame(
        {"PDS_acc": cluster_list_db, "target_acc": centroid_list_db}
    )
    df_cluster_center.to_csv(
        os.path.join(tmp_folder, f"{pdg_acc}_cluster_center.tsv"), sep="\t", index=False
    )
    time_end_calculate_centroid = time.time()
    logging.info(
        f"Calculated centroid in {round(time_end_calculate_centroid-time_start_calculate_centroid,2)} seconds."
    )
    logging.info(f"Number of clusters: {len(cluster_list_db)}")


def download_and_sketch_assembly(gca_acc_list, hash_number, kmer_size, tmp_folder):
    logging.info("Downloading and sketching assemblies...")
    time_start_download_and_sketch_assembly = time.time()
    gca_acc_file = os.path.join(tmp_folder, "gca_list.txt")
    with open(gca_acc_file, "w") as f:
        for item in gca_acc_list:
            f.write("%s\n" % item)
    subprocess.run(
        [
            "datasets",
            "download",
            "genome",
            "accession",
            "--inputfile",
            gca_acc_file,
            "--filename",
            "assembly.zip",
            "--dehydrated",
        ]
    )
    # move the zip file to the tmp folder
    zip_file = os.path.join(tmp_folder, "assembly.zip")
    shutil.move("assembly.zip", zip_file)
    with zipfile.ZipFile(zip_file, "r") as zip_ref:
        zip_ref.extractall(os.path.join(tmp_folder, "assembly"))
    # try to rehydrate the dataset, if failed, try again after 5 seconds
    while True:
        try:
            subprocess.run(
                [
                    "datasets",
                    "rehydrate",
                    "--directory",
                    os.path.join(tmp_folder, "assembly"),
                ],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            break
        except subprocess.CalledProcessError:
            time.sleep(5)
            logging.error(f"Failed to rehydrate {tmp_folder}/assembly. Trying again...")
    time_end_download = time.time()
    logging.info(
        f"Downloaded {len(gca_acc_list)} assemblies in {round(time_end_download-time_start_download_and_sketch_assembly,2)} seconds."
    )
    os.mkdir(os.path.join(tmp_folder, "signature"))
    for acc in gca_acc_list:
        # TODO: other file name patterns?
        fna_file_list = glob.glob(
            os.path.join(
                tmp_folder, "assembly", "ncbi_dataset", "data", acc, "*_genomic.fna"
            )
        )
        # skip if there are more than one fna files or no fasta file
        if len(fna_file_list) != 1:
            continue
        mh = MinHash(n=hash_number, ksize=kmer_size)
        with screed.open(fna_file_list[0], "r") as f:
            for record in f:
                mh.add_sequence(record.sequence, True)
        signame = acc
        sig = [SourmashSignature(mh, name=signame)]
        with open(os.path.join(tmp_folder, "signature", acc + ".sig"), "w") as f:
            save_signatures(sig, fp=f)
    time_end_sketch = time.time()
    logging.info(
        f"Sketched {len(gca_acc_list)} assemblies in {round(time_end_sketch-time_end_download,2)} seconds."
    )


def merge_sig(args, db_folder, tmp_folder):
    logging.info("Merging all signature files...")
    time_start_merge_sig = time.time()
    all_sig_path = os.path.join(tmp_folder, "signature", "*.sig")
    sig_path_list = glob.glob(all_sig_path)
    siglist = []
    for sig_path in sig_path_list:
        siglist.append(load_one_signature(sig_path))
    sig_file_name = os.path.join(db_folder, args.name + ".sig")
    with open(sig_file_name, "w") as f:
        save_signatures(siglist, fp=f)
    time_end_merge_sig = time.time()
    logging.info(
        f"Merged all signature files in {round(time_end_merge_sig-time_start_merge_sig,2)} seconds."
    )


# import metadata to database directly from pathogen detection metadata dataframe
def import_metadata(df_metadata, df_cluster_center_metadata, conn):
    info = {
        "biosample_acc": None,
        "taxid": None,
        "strain": None,
        "collected_by": None,
        "collection_date": None,
        "geo_loc_name": None,
        "isolation_source": None,
        "lat_lon": None,
        "serovar": None,
        "sub_species": None,
        "species": None,
        "genus": None,
        "host": None,
        "host_disease": None,
        "outbreak": None,
        "srr": None,
        "PDT_acc": None,
        "PDS_acc": None,
        "asm_acc": None,
    }

    for index, row in df_cluster_center_metadata.iterrows():
        info["biosample_acc"] = row["biosample_acc"]
        info["taxid"] = row["taxid"]
        info["strain"] = row["strain"]
        info["collected_by"] = row["collected_by"]
        info["collection_date"] = row["collection_date"]
        info["geo_loc_name"] = row["geo_loc_name"]
        info["isolation_source"] = row["isolation_source"]
        info["lat_lon"] = row["lat_lon"]
        info["serovar"] = row["serovar"]
        info["host"] = row["host"]
        info["host_disease"] = row["host_disease"]
        info["outbreak"] = ",".join(
            list(
                df_metadata[df_metadata["PDS_acc"] == row["PDS_acc"]]
                .dropna(subset=["outbreak"])
                .outbreak.unique()
            )
        )
        info["srr"] = row["Run"]
        info["PDT_acc"] = row["target_acc"]
        info["PDS_acc"] = row["PDS_acc"]
        info["asm_acc"] = row["asm_acc"]
        # check whether the information is missing
        for key in info:
            if info[key] is None:
                info[key] = "missing"

        info_list = list(info.values())

        with conn:
            insert_metadata(conn, info_list)
    return


def insert_metadata(conn, info):
    sql = """INSERT OR IGNORE INTO METADATA(
             biosample_acc,
             taxid,
             strain,
             collected_by,
             collection_date,
             geo_loc_name,
             isolation_source,
             lat_lon,
             serovar,
             sub_species,
             species,
             genus,
             host,
             host_disease,
             outbreak,
             srr,
             PDT_acc,
             PDS_acc,
             asm_acc) 
             VALUES(:biosample_acc, :taxid, :strain, :collected_by, :collection_date,
                    :geo_loc_name, :isolation_source, :lat_lon, :serovar, :sub_species,
                    :species, :genus, :host, :host_disease, :outbreak, :srr,
                    :PDT_acc, :PDS_acc, :asm_acc);"""
    c = conn.cursor()
    c.execute(sql, info)


def prepare(args):
    # check if datasets is installed
    if not shutil.which("datasets"):
        logging.error(
            "NCBI datasets is not installed. Please install it from https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/"
        )
        exit(1)
    t = time.localtime()
    current_time = time.strftime("%Y%m%d%H%M%S", t)
    print_handler = logging.StreamHandler()
    file_handler = logging.FileHandler(f"mashpit-{current_time}.log")
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

    # print out database infomation
    logging.info("-" * 50)
    if args.type == "taxon":
        logging.info("Database type: Taxon")
        if args.species is None:
            logging.error("Taxon name not provided")
            exit(1)
        logging.info(f"Species: {args.species}")
    if args.type == "accession":
        logging.info("Database type: Accession")
        if args.list is None:
            logging.error("Accession list not provided")
            exit(1)
        logging.info(f"Biosample accession list: {args.list}")
        if args.email == None:
            logging.error("Entrez email not provided.")
            exit(1)
    if args.name is None:
        logging.error("Database name not provided")
        exit(1)
    logging.info(f"Database name: {args.name}")
    logging.info(f"Maximum number of hashes: {args.number}")
    logging.info(f"Kmer size: {args.ksize}")
    logging.info("-" * 50)

    # create a database folder
    cwd = os.getcwd()
    db_folder = os.path.join(cwd, args.name)
    if os.path.isdir(db_folder):
        logging.error(f"Folder {db_folder} already exists.")
        exit(1)
    os.mkdir(db_folder)
    logging.info(f"Database folder created at {db_folder}")

    # create a tmp folder
    tmp_folder = os.path.join(cwd, "tmp")
    if os.path.isdir(tmp_folder):
        shutil.rmtree(tmp_folder)
    os.mkdir(tmp_folder)

    # create a sqlite database
    sql_path = os.path.join(db_folder, f"{args.name}.db")
    conn = create_connection(sql_path)
    create_database(conn)

    return db_folder, tmp_folder, conn


def build_taxon(args):
    # prepare the database folder and sqlite database
    db_folder, tmp_folder, conn = prepare(args)
    log_memory_usage("start")
    # format the pathogen name
    pathogen_name = args.species.strip().replace(" ", "_")
    # validate the species name
    results_response = urllib.request.urlopen(
        "https://ftp.ncbi.nlm.nih.gov/pathogen/Results/"
    )
    results_html = results_response.read().decode("utf-8")
    results_parser = Species_Name_Parser()
    results_parser.feed(results_html)
    pathogen_name_lower = pathogen_name.lower()
    for item in results_parser.list:
        # convert the list item to lowercase
        item_lower = item.lower()
        # compare the lowercase versions of the string and the list item
        if pathogen_name_lower == item_lower:
            pathogen_name = item
            break
    else:
        logging.error(
            f"Taxon {pathogen_name} not available on Pathogen Detection. Check for available names at https://ftp.ncbi.nlm.nih.gov/pathogen/Results/"
        )
        exit(1)
    logging.info(f"Taxon name validated. Using {pathogen_name} as the taxon name.")
    # download the latest PD metadata files and get the PDG accession
    metadata_file_name, isolate_pds_file = download_metadata(
        pathogen_name, args.pd_version, tmp_folder
    )
    log_memory_usage("After downloading files")
    df_metadata = pd.read_csv(
        os.path.join(tmp_folder, metadata_file_name), header=0, sep="\t"
    )
    df_isolate_pds = pd.read_csv(
        os.path.join(tmp_folder, isolate_pds_file), header=0, sep="\t"
    )
    # add cluster accession to metadata
    df_metadata = df_metadata.join(
        df_isolate_pds[["target_acc", "PDS_acc"]].set_index("target_acc"),
        on="target_acc",
    )
    pdg_acc = metadata_file_name.replace(".metadata.tsv", "")
    # filter out isolates without genome assembly
    df_metadata_asm = df_metadata[~df_metadata["asm_acc"].isnull()]
    # calculate the centroid for each cluster
    calculate_centroid(df_metadata_asm, pdg_acc, tmp_folder)
    df_cluster_center = pd.read_csv(
        os.path.join(tmp_folder, f"{pdg_acc}_cluster_center.tsv"), sep="\t"
    )
    df_cluster_center_metadata = df_cluster_center.join(
        df_metadata_asm.drop("PDS_acc", axis=1).set_index("target_acc"), on="target_acc"
    )
    log_memory_usage("After centroid calculation")
    # download the centroid assembly using NCBI datasets
    gca_acc_list = df_cluster_center_metadata["asm_acc"].to_list()
    hash_number = args.number
    kmer_size = args.ksize
    download_and_sketch_assembly(gca_acc_list, hash_number, kmer_size, tmp_folder)
    log_memory_usage("After downloading assemblies")
    # merge all signature files
    merge_sig(args, db_folder, tmp_folder)
    log_memory_usage("After merging signature files")
    # import metadata to database
    import_metadata(df_metadata, df_cluster_center_metadata, conn)
    log_memory_usage("After importing metadata")
    # add description to database
    c = conn.cursor()
    c.execute("INSERT INTO DESC VALUES ('Type','Taxonomy');")
    c.execute(f"INSERT INTO DESC VALUES ('Species','{pathogen_name}');")
    c.execute(f"INSERT INTO DESC VALUES ('Version','{pdg_acc}');")
    c.execute(f"INSERT INTO DESC VALUES ('Hash_number','{str(hash_number)}');")
    c.execute(f"INSERT INTO DESC VALUES ('Kmer_size','{str(kmer_size)}');")
    conn.commit()

    # delete temp files and folders
    logging.info("Deleting temp file and folders")
    shutil.rmtree(tmp_folder)
    logging.info("Complete.")


def build_accession(args):
    db_folder, tmp_folder, conn = prepare(args)
    Entrez.email = args.email
    if args.key != None:
        Entrez.api_key = args.key
    # check if the list file exists
    if not os.path.isfile(args.list):
        logging.error(f"File {args.list} not found.")
        exit(1)
    # download assembly and get metadata
    gca_acc_list = []
    f = open(args.list, "r")
    biosample_list = f.readlines()
    for biosample in biosample_list:
        # if api key is not provided, wait for 1 second
        if args.key == None:
            time.sleep(1)
        # check if the assembly is available on NCBI
        esearch_handle = Entrez.esearch(db="assembly", term=biosample, retmax="3")
        esearch_record = Entrez.read(esearch_handle)
        # skip if no assembly found
        if len(esearch_record["IdList"]) == 0:
            continue
        # get assembly accession
        if args.key == None:
            time.sleep(1)
        esummary_handle = Entrez.esummary(
            db="assembly", id=esearch_record["IdList"][0], report="full"
        )
        esummary_record = Entrez.read(esummary_handle)
        asm_acc = esummary_record["DocumentSummarySet"]["DocumentSummary"][0][
            "AssemblyAccession"
        ]
        gca_acc_list.append(asm_acc)
        # get metadata
                # if api key is not provided, wait for 1 second
        if args.key == None:
            time.sleep(1)
        esearch_handle = Entrez.esearch(db="biosample", term=biosample, retmax="3")
        esearch_record = Entrez.read(esearch_handle)
        if args.key == None:
            time.sleep(1)
        efetch_handle = Entrez.efetch(db="biosample", id=esearch_record["IdList"][0])
        efetch_record = efetch_handle.read()
        root = ET.fromstring(efetch_record)
        info = {
            "biosample_acc": None,
            "taxid": None,
            "strain": None,
            "collected_by": None,
            "collection_date": None,
            "geo_loc_name": None,
            "isolation_source": None,
            "lat_lon": None,
            "serovar": None,
            "sub_species": None,
            "species": None,
            "genus": None,
            "host": None,
            "host_disease": None,
            "outbreak": None,
            "srr": None,
            "PDT_acc": None,
            "PDS_acc": None,
            "asm_acc": None,
        }
        info["biosample_acc"] = biosample
        info["asm_acc"] = asm_acc
        for domain in root[0]:
            if domain.tag == "Attributes":
                for sub_attrib in domain:
                    if sub_attrib.attrib["attribute_name"] in info:
                        info[sub_attrib.attrib["attribute_name"]] = sub_attrib.text
        insert_metadata(conn, info)
    f.close()
    hash_number = args.number
    kmer_size = args.ksize
    download_and_sketch_assembly(gca_acc_list, hash_number, kmer_size, tmp_folder)
    # merge all generated signature files
    merge_sig(args, db_folder, tmp_folder)
    c = conn.cursor()
    c.execute("INSERT INTO DESC VALUES ('Type','Accession');")
    c.execute("INSERT INTO DESC VALUES ('Hash_number','" + str(hash_number) + "');")
    c.execute("INSERT INTO DESC VALUES ('Kmer_size','" + str(kmer_size) + "');")
    conn.commit()
    logging.info("Deleting tmp folders")
    shutil.rmtree(tmp_folder)
    logging.info("Complete.")


def build(args):
    if args.type == "taxon":
        build_taxon(args)
    elif args.type == "accession":
        build_accession(args)
