#!/usr/bin/env python3

import sqlite3
import subprocess
import xml.etree.ElementTree as ET
from sqlite3 import Error
from Bio import Entrez

## Define the methods to create the database and the tables 
def create_connection(db_file):
    conn = None
    try:
        conn = sqlite3.connect(db_file)
        return conn
    except Error as e:
        print(e)
 
    return conn

def create_table(conn, create_table_sql):
    try:
        c = conn.cursor()
        c.execute(create_table_sql)
    except Error as e:
        print(e)
    

def main():
    ## define and create the tables in the database
    conn = create_connection('mashpit.db')
    
                     
    sql_create_biosample = """CREATE TABLE IF NOT EXISTS biosample (
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
    
    sql_create_sra = """CREATE TABLE IF NOT EXISTS sra (
                          srr           TEXT PRIMARY KEY,
                          biosample_acc TEXT
                        );"""
    
    sql_create_sketch = """CREATE TABLE IF NOT EXISTS sketch (
                           sketchid      INTEGER PRIMARY KEY AUTOINCREMENT,
                           biosample_acc TEXT NOT NULL,
                           srr           TEXT NOT NULL,
                           path          TEXT NOT NULL,
                           source        TEXT NOT NULL, 
                           software      TEXT NOT NULL, 
                           seed          INTEGER NOT NULL
                        );"""    

    if conn is not None:
        create_table(conn, sql_create_biosample)
        create_table(conn, sql_create_sra)
        create_table(conn, sql_create_sketch)
    else:
        print("Cannot create the database connection.")


    


if __name__ == '__main__':
  main()
