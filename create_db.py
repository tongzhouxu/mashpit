#!/usr/bin/env python3

import sqlite3
from sqlite3 import Error


# Define the methods to create the database
def create_connection(db_file):
    conn = None
    try:
        conn = sqlite3.connect(db_file)
        return conn
    except Error as e:
        print(e)

    return conn


# Define the methods to create the tables
def create_table(conn, create_table_sql):
    try:
        c = conn.cursor()
        c.execute(create_table_sql)
    except Error as e:
        print(e)


def main():
    conn = create_connection('mashpit.db')

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
                              genotype         TEXT,
                              host             TEXT,
                              host_disease     TEXT
                        );"""

    sql_create_sra = """CREATE TABLE IF NOT EXISTS SRA (
                          srr           TEXT PRIMARY KEY,
                          biosample_acc TEXT,
                          FOREIGN KEY (biosample_acc) REFERENCES BIOSAMPLE(biosample_acc)
                          ON DELETE CASCADE
                          ON UPDATE CASCADE
                        );"""

    sql_create_sketch = """CREATE TABLE IF NOT EXISTS SKETCH (
                           sketchid      INTEGER PRIMARY KEY AUTOINCREMENT,
                           biosample_acc TEXT NOT NULL,
                           srr           TEXT NOT NULL,
                           path          TEXT NOT NULL,
                           source        TEXT NOT NULL, 
                           software      TEXT NOT NULL, 
                           seed          INTEGER NOT NULL,
                           UNIQUE    (path),
                           FOREIGN KEY (srr) REFERENCES SRA(srr)
                           ON DELETE CASCADE
                           ON UPDATE CASCADE,
                           FOREIGN KEY (biosample_acc) REFERENCES BIOSAMPLE(biosample_acc)
                           ON DELETE CASCADE
                           ON UPDATE CASCADE
                        );"""

    if conn is not None:
        create_table(conn, sql_create_biosample)
        create_table(conn, sql_create_sra)
        create_table(conn, sql_create_sketch)
        conn.commit()
        conn.close()
    else:
        print("Cannot create the database connection.")


if __name__ == '__main__':
    main()
