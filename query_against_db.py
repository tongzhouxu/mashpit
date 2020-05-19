#!/usr/bin/env python3

import argparse
import csv
import subprocess
from create_db import create_connection
from create_db import create_table


def parse_args():
    parser = argparse.ArgumentParser(usage='query_against_db.py -n <sample name>')
    parser.add_argument("-n", help="<string>: sample SRR name")
    return parser.parse_args()


def insert_distance(sample_name, conn, info):
    sql = ''' INSERT INTO ''' + sample_name + '''_distance(biosample_acc,mash_distance)
              VALUES(?,?) '''
    cur = conn.cursor()
    cur.execute(sql, info)
    return cur.lastrowid


def main():
    conn = create_connection('mashpit.db')
    args = parse_args()
    sample_name = args.n

    sql_create_distance = """CREATE TABLE IF NOT EXISTS """ + sample_name + """_distance (
                              biosample_acc    TEXT PRIMARY KEY, 
                              mash_distance    REAL
                      );"""
    sql_create_output = ("""CREATE TABLE IF NOT EXISTS """ + sample_name + """_output (
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
                                  mash_distance    REAL
                          )""")
    if conn is not None:
        create_table(conn, sql_create_distance)
        create_table(conn, sql_create_output)
    else:
        print("Cannot create the database connection.")
    conn.commit()

    c = conn.cursor()
    cursor = c.execute("SELECT path from SKETCH")
    msh = []
    distance = []
    for row in cursor:
        msh.append(row[0])
    # get the list of the distances
    for i in msh:
        res = subprocess.check_output(['mash', 'dist', i, sample_name + '_skesa.fa.msh'])
        res = res.decode('utf-8')
        distance.append(res.split()[2])
    cursor = c.execute("SELECT biosample_acc from SKETCH")
    biosample_acc_list = []
    for row in cursor:
        biosample_acc_list.append(row[0])
    for i in range(len(msh)):
        insert_distance(sample_name, conn, [biosample_acc_list[i], distance[i]])
    # combine the tables
    c.execute("""INSERT INTO """ + sample_name + """_output SELECT 
                 biosample.biosample_acc,biosample.taxid,
                 biosample.strain,
                 biosample.collected_by,
                 biosample.collection_date,
                 biosample.geo_loc_name,
                 biosample.isolation_source,
                 biosample.lat_lon,
                 biosample.genotype,
                 biosample.host,
                 biosample.host_disease,
                 """ + sample_name + """_distance.mash_distance 
                 From biosample
                 INNER JOIN """ + sample_name + """_distance
                 ON biosample.biosample_acc=""" + sample_name + """_distance.biosample_acc""")

    print("Exporting data into CSV……")
    c.execute("SELECT * FROM " + sample_name + "_output ORDER BY mash_distance ASC")
    with open(sample_name + "output.csv", "w") as csv_file:
        csv_writer = csv.writer(csv_file, delimiter="\t")
        csv_writer.writerow([i[0] for i in c.description])
        csv_writer.writerows(c)

    conn.commit()


if __name__ == '__main__':
    main()
