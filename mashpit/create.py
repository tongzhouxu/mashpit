#!/usr/bin/env python3

import os
import sqlite3
from sqlite3 import Error

# Define the method to create the database and connect to the database
def create_connection(db_file):
    conn = None
    try:
        conn = sqlite3.connect(db_file, check_same_thread=False)
        return conn
    except Error as e:
        print(e)

    return conn


# Define the method to create the tables
def create_table(conn, create_table_sql):
    try:
        c = conn.cursor()
        c.execute(create_table_sql)
    except Error as e:
        print(e)


def create(args):
    cwd = os.getcwd()
    db_path = os.path.join(cwd, args.database + '.db')
    conn = create_connection(db_path)

    # define and create the tables in the database
    sql_create_biosample = """CREATE TABLE IF NOT EXISTS BIOSAMPLE (
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
                              outbreak         TEXT
                        );"""

    sql_create_sra = """CREATE TABLE IF NOT EXISTS SRA (
                          srr           TEXT PRIMARY KEY,
                          biosample_acc TEXT,
                          FOREIGN KEY (biosample_acc) REFERENCES BIOSAMPLE(biosample_acc)
                          ON DELETE CASCADE
                          ON UPDATE CASCADE
                        );"""

    if conn is not None:
        create_table(conn, sql_create_biosample)
        create_table(conn, sql_create_sra)
        conn.commit()
        conn.close()
    else:
        print("Cannot create the database connection.")
