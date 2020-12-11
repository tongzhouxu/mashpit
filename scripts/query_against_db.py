#!/usr/bin/env python3

import argparse
import csv
import os
import ntpath
import pandas as pd
from scripts.create_db import create_connection
from scripts.create_db import create_table
import screed
import sourmash
from sourmash import SourmashSignature, save_signatures, load_one_signature, load_signatures


def parse_args():
    parser = argparse.ArgumentParser(usage='query_against_db.py <sample name> <database name> [--force, -f]')
    parser.add_argument("sample", help="<string>: path to the sample file")
    parser.add_argument("database", help="<string>: name of the database")
    parser.add_argument("-f", "--force", help="overwrite if query table exists", action="store_true")
    return parser.parse_args()


def check_output_existence(conn, args, sample_name):
    c = conn.cursor()
    c.execute("""SELECT count(name) FROM sqlite_master WHERE type='table' AND name='""" + sample_name + """_output' """)
    if c.fetchone()[0] == 1:
        if args.force:
            c.execute("""DROP TABLE """ + sample_name + """_output""")
            c.execute("""DROP TABLE """ + sample_name + """_distance""")
        else:
            print("Query result exists. Use option -f to overwrite the result")
            exit()
    c.close()
    conn.commit()


def get_target_sig(sample_name):
    genome = sample_name
    mh = sourmash.MinHash(n=1000, ksize=21)
    for record in screed.open(genome):
        mh.add_sequence(record.sequence, True)
    sig = SourmashSignature(mh, name=genome)
    with open(sample_name + '.sig', 'wt') as fp:
        save_signatures([sig], fp)


def select_by_srr(conn, srr):
    c = conn.cursor()
    cursor = c.execute("SELECT biosample_acc FROM SRA WHERE srr=?", (srr,))
    return cursor.fetchone()[0]


def insert_distance(sample_name, conn, info):
    sql = ''' INSERT INTO ''' + sample_name + '''_distance(biosample_acc,jaccard_similarity)
              VALUES(?,?) '''
    cur = conn.cursor()
    cur.execute(sql, info)
    return cur.lastrowid


def main():
    args = parse_args()
    sample_path = args.sample
    sample_name = ntpath.basename(sample_path)

    cwd = os.getcwd()
    db_path = os.path.join(cwd, args.database + '.db')
    database_sig_path = os.path.join(cwd, args.database + '.sig')
    target_sig_path = os.path.join(cwd, sample_path + '.sig')

    # check for the existence of the database and tables
    if os.path.exists(db_path):
        pass
    else:
        print("Database does not exist. Please make sure the name is correct or run create_db.py and "
              "metadata_sra_db.py first")
        exit(0)
    conn = create_connection(db_path)
    c = conn.cursor()
    c.execute(''' SELECT count(name) FROM sqlite_master WHERE type='table' AND name='BIOSAMPLE' ''')
    if c.fetchone()[0] == 0:
        print("No BIOSAMPLE table found in the database. Please make sure the name is correct or run create_db.py and "
              "metadata_sra_db.py first")
        exit(0)

    sql_create_distance = """CREATE TABLE IF NOT EXISTS """ + sample_name + """_distance (
                              biosample_acc         TEXT PRIMARY KEY, 
                              jaccard_similarity    REAL
                      );"""
    sql_create_output = """CREATE TABLE IF NOT EXISTS """ + sample_name + """_output (
                                  biosample_acc    TEXT PRIMARY KEY, 
                                  taxid            INTEGER,
                                  strain           TEXT, 
                                  collected_by     TEXT,
                                  collection_date  TEXT,
                                  geo_loc_name     TEXT,
                                  isolation_source TEXT,
                                  lat_lon          TEXT,
                                  genotype         TEXT,
                                  host             TEXT,
                                  host_disease     TEXT,
                                  outbreak         TEXT,
                                  jaccard_similarity    REAL
                          );"""

    check_output_existence(conn, args, sample_name)
    create_table(conn, sql_create_distance)
    create_table(conn, sql_create_output)
    conn.commit()
    get_target_sig(sample_path)

    database_sig = load_signatures(database_sig_path)
    target_sig = load_one_signature(target_sig_path)
    c = conn.cursor()
    for sig in database_sig:
        biosample_acc = select_by_srr(conn,sig.name())
        distance = target_sig.jaccard(sig)
        insert_distance(sample_name, conn, [biosample_acc, distance])

    # combine the tables
    c.execute("""INSERT INTO """ + sample_name + """_output SELECT 
                 biosample.biosample_acc,biosample.taxid,
                 biosample.strain,
                 biosample.collected_by,
                 biosample.collection_date,
                 biosample.geo_loc_name,
                 biosample.isolation_source,
                 biosample.lat_lon,
                 biosample.serovar,
                 biosample.host,
                 biosample.host_disease,
                 biosample.outbreak,
                 """ + sample_name + """_distance.jaccard_similarity 
                 From BIOSAMPLE
                 INNER JOIN """ + sample_name + """_distance
                 ON biosample.biosample_acc=""" + sample_name + """_distance.biosample_acc""")
    print("Printing out the top 50 results.")
    print(pd.read_sql_query("SELECT * FROM " + sample_name + "_output ORDER BY jaccard_similarity DESC LIMIT 50", conn))
    print("Output file has been stored in the database table " + sample_name + "_output and exported in CSV format.")
    c.execute("SELECT * FROM " + sample_name + "_output ORDER BY jaccard_similarity DESC")
    with open(sample_name + "_output.csv", "w") as csv_file:
        csv_writer = csv.writer(csv_file, delimiter="\t")
        csv_writer.writerow([i[0] for i in c.description])
        csv_writer.writerows(c)

    conn.commit()


if __name__ == '__main__':
    main()
