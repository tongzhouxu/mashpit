#!/usr/bin/env python3
# Author: Tongzhou Xu

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

## Define a method to insert information into the biosample database
def insert_biosample(conn, info):
    
    sql = ''' INSERT INTO biosample(biosample_acc,srr,taxid,strain,collected_by,collection_date,geo_loc_name,isolation_source,lat_lon,genotype,host,host_disease,BioSampleModel)
              VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?) '''
    cur = conn.cursor()
    cur.execute(sql, info)
    return cur.lastrowid

## Define a method to insert information into the sketch database
def insert_sketch(conn, info):
    sql = ''' INSERT INTO sra(sketchid,biosample_acc,srr,path,source,software,seed)
              VALUES(?,?,?,?,?,?,?) '''
    cur = conn.cursor()
    cur.execute(sql,info)
    return cur.lastrowid

## Using a list of bioproject in PDP
## First, find the names of the biosamples in a bioproject and insert the information into the sql database
def biosample_sra_database(bioproject_acc,conn):
  ## find the ncbi-id for the specific bioproject
  handle_search = Entrez.esearch(db="bioproject",term=bioproject_acc)
  record_search = Entrez.read(handle_search)
  id_list = record_search['IdList']
  ## find the related biosample id list
  handle_link = Entrez.elink(db='biosample',dbfrom='bioproject',id=id_list[0])
  record_link = Entrez.read(handle_link)
  ## parse all the ids
  for record in record_link[0]['LinkSetDb'][0]['Link']:
      handle_link_sra = Entrez.elink(db='sra',dbfrom='biosample',id=record['Id'])
      record_link_sra = Entrez.read(handle_link_sra)
      handle_fetch_sra = Entrez.efetch(db='sra',id=record_link_sra[0]['LinkSetDb'][0]['Link'][0]['Id'])
      xml_result_sra = handle_fetch_sra.read()
      root_sra = ET.fromstring(xml_result_sra)
      info = []
      biosample_acc = root_sra.find('EXPERIMENT_PACKAGE').find('EXPERIMENT').find('DESIGN').find('SAMPLE_DESCRIPTOR').find('IDENTIFIERS').find('EXTERNAL_ID')
      info.append(biosample_acc.text)
      srr = root_sra.find('EXPERIMENT_PACKAGE').find('RUN_SET').find('RUN')
      info.append(srr.attrib['accession'])
      taxid = root_sra.find('EXPERIMENT_PACKAGE').find('Pool').find('Member')
      info.append(taxid.attrib['tax_id'])
      attributes = root_sra.find('EXPERIMENT_PACKAGE').find('SAMPLE').find('SAMPLE_ATTRIBUTES')
      for item in attributes:
          info.append(item.find('VALUE').text)
      insert_biosample(conn, info)  
      
## Download the assembly according to the sra database
def skesa_assembly_download(SRR):
    ## TO DO: check whether the directory and the assembly already exist
    subprocess.check_call("dump-ref-fasta http://sra-download.ncbi.nlm.nih.gov/srapub_files/"+SRR+"_"+SRR+".realign > "+SRR+"_skesa.fa",shell=True) 
    ## if the skesa assembly is not avalaible, download the fastq file and make the assembly using skesa
    """
    subprocess.check_call("fastq-dump --split-lanes "+SRR,shell=True)
    subprocess.check_call("skesa --reads "+SRR+"_1.fq,"+SRR+"_2.fq --cores 4 --memory 16 > "+SRR+".skesa.fa",shell=True)
    """

## sketch the assembly and upload it into the sketch database
def sketch_assembly(SRR):
    subprocess.check_call("mash sketch "+SRR+"_skesa.fa",shell=True)

## Cluster detection using mash
def mash_distance(sketch1,sketch2):
    subprocess.check_call("mash dist "+sketch1+" "+sketch2,shell=True)

def main():
    Entrez.email = "tongzhouxu97@gmail.com"
    Entrez.api_key = "44c9216e8c0c24c97fa7871093da74808908"
    
    ## define and create the tables in the database
    conn = create_connection('mashpit.db')
    
                     
    sql_create_biosample = """CREATE TABLE IF NOT EXISTS biosample (
                              biosample_acc    TEXT PRIMARY KEY, 
                              srr              TEXT,
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
                              BioSampleModel   TEXT
                        );"""

    sql_create_taxonomy = """CREATE TABLE IF NOT EXISTS taxonomy (
                             taxid     INTEGER PRIMARY KEY, 
                             genus      TEXT, 
                             species    TEXT, 
                             subspecies TEXT
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
        create_table(conn, sql_create_taxonomy)
        create_table(conn, sql_create_biosample)
        create_table(conn, sql_create_sketch)
    else:
        print("Cannot create the database connection.")
    
    f = open("project_list.txt",'r')           
    project_list = f.readlines()
    for project in project_list:
        project = project.strip('\n')
        biosample_sra_database(project,conn)

    c = conn.cursor()
    cursor = c.execute("SELECT srr  from biosample")
    SRR = []
    for row in cursor:
        SRR.append()
        skesa_assembly_download(row[0])
        sketch_assembly(row[0])
    
    SRR2 = SRR
    for i in SRR:
        for j in SRR2:
            mash_distance(i+"_skesa.fa",j+"_skesa.fa")
 

if __name__ == '__main__':
  main()
