#!/usr/bin/env python3

import sqlite3
import subprocess
from sqlite3 import Error
from create_db import create_connection
from create_db import create_table

## Download the assembly according to the sra database
def skesa_assembly_download(SRR):
    subprocess.check_call("cd ./skesa_assem",shell=True)
    try:
        res = subprocess.check_call("dump-ref-fasta http://sra-download.ncbi.nlm.nih.gov/srapub_files/"+SRR+"_"+SRR+".realign > "+SRR+"_skesa.fa",shell=True)
        print("res: ",res)
    except subprocess.CalledProcessError:
        print("can't download skesa assembly for "+SRR)
        try:
            subprocess.check_call("fastq-dump --split-lanes "+SRR,shell=True)
            subprocess.check_call("skesa --reads "+SRR+"_1.fq,"+SRR+"_2.fq --cores 4 --memory 16 > "+SRR+".skesa.fa",shell=True)
        except subprocess.CalledProcessError:
            subprocess.check_call("skesa --reads "+SRR+"_1.fq,"+SRR+"_2.fq --cores 4 --memory 16 > "+SRR+".skesa.fa",shell=True)
    subprocess.check_call("cd ..",shell=True)
    ## if the skesa assembly is not avalaible, download the fastq file and make the assembly using skesa


def insert_sketch(conn, info):
    
    sql = ''' INSERT INTO sketch (biosample_acc,srr,path,source,software,seed)
              VALUES(?,?,?,?,?,?) '''
    cur = conn.cursor()
    cur.execute(sql, info)
    return cur.lastrowid
    
def skech_assembly(SRR):
    try:
        subprocess.check_call("mash sketch ./skesa_assem/"+SRR+"_skesa.fa",shell=True)
        info = {"biosample_acc":None,"srr":None,"path":None,"source":None,"software":None,"seed":None}
        info["srr"] = SRR
        info["path"] = "./skesa_assem/"+SRR+"_skesa.fa.msh"
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
 
def main():
    conn = create_connection('mashpit.db')
    c = conn.cursor()
    cursor = c.execute("SELECT srr from sra")
    for row in cursor:
        """
        skesa_assembly_download(row[0])
        """
        info=skech_assembly(row[0])
        info["biosample_acc"] = select_by_srr(conn,row[0])
        info_list=list(info.values())
        insert_sketch(conn,info_list)

    

if __name__ == '__main__':
  main()