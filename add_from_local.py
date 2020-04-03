#!/usr/bin/env python3

import sqlite3
import subprocess
import argparse
from create_db import create_connection
from metadata_sra_db import insert_biosample
from metadata_sra_db import insert_sra

## local paired reads needed to be named as SRR_1.fastq and SRR_2.fastq 
## all the metadata information should be in a separate list

def parse_args():
    parser = argparse.ArgumentParser(usage='add_from_local.py -metadata <metadata information> -srr <srr information>')
    parser.add_argument("-metadata",help="<string>: full name of the metadata info file")
    parser.add_argument("-srr",help="<string>: full name of the srr info file")
    return parser.parse_args()

def main():
    args = parse_args()
    conn = create_connection('mashpit.db')
    metadata_file = open(args.metadata,'r')
    srr_file = open(args.srr,'r')
    metadata_list = metadata_file.readlines()
    srr_list = srr_file.readlines()
    for metadata in metadata_list:
        metadata = metadata.strip('\n').split
        insert_biosample(conn,metadata)
        conn.commit()
    for srr in srr_list:
        srr = srr.strip('\n').split
        insert_sra(conn,srr)
        conn.commit()
    
    metadata_file.close()
    srr_file.close()
    conn.close()




if __name__ == '__main__':
  main()