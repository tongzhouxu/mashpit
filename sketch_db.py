#!/usr/bin/env python3

import sqlite3
import subprocess
from sqlite3 import Error
from create_db import create_connection
from create_db import create_table

## Download the assembly according to the sra database
def skesa_assembly_download(SRR):
    try:
        subprocess.check_call("dump-ref-fasta http://sra-download.ncbi.nlm.nih.gov/srapub_files/"+SRR+"_"+SRR+".realign > "+SRR+"_skesa.fa",shell=True)
        print("Downloaded assembly for "+SRR)
    except subprocess.CalledProcessError:
        print("can't download skesa assembly for "+SRR)
        try:
            subprocess.check_call("prefetch "+SRR+" -O ./", shell=True)
            subprocess.check_call("fastq-dump --split-files "+SRR+".sra",shell=True)
            subprocess.check_call("skesa --reads "+SRR+"_1.fastq,"+SRR+"_2.fastq --cores 4 --memory 8 > "+SRR+"_skesa.fa",shell=True)
        except subprocess.CalledProcessError:          
            subprocess.check_call("skesa --reads "+SRR+"_1.fastq,"+SRR+"_2.fastq --cores 4 --memory 8 > "+SRR+"_skesa.fa",shell=True)

    
def skech_assembly(SRR):
    try:
        subprocess.check_call("mash sketch "+SRR+"_skesa.fa",shell=True)
        info = {"sketchid":None,"biosample_acc":None,"srr":None,"path":None,"source":None,"software":None,"seed":None}
        info["srr"] = SRR
        info["path"] = SRR+"_skesa.fa.msh"
        info["source"] = "NCBI"
        info["software"] = "mash"
        info["seed"] = 42
    except subprocess.CalledProcessError:
        print("no assembly found")
    return info

def select_by_srr(conn,srr):
    c = conn.cursor()
    cursor = c.execute("SELECT biosample_acc FROM sra WHERE srr=?", (srr,))
    return cursor.fetchone()[0]


def insert_sketch(conn, info):
    
    sql = ''' INSERT INTO sketch(sketchid,biosample_acc,srr,path,source,software,seed)
              VALUES(?,?,?,?,?,?,?) '''
    cur = conn.cursor()
    cur.execute(sql, info)
    return cur.lastrowid 

def main():
    conn = create_connection('mashpit.db')
    c = conn.cursor()
    cursor = c.execute("SELECT srr from sra")
    for row in cursor:
        skesa_assembly_download(row[0])
        info=skech_assembly(row[0])
        info["biosample_acc"] = select_by_srr(conn,row[0])
        info_list=list(info.values())
        print(info_list)
        insert_sketch(conn,info_list)
        conn.commit()
    
    conn.close()

    

if __name__ == '__main__':
  main()