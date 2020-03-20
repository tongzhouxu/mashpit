#!/usr/bin/env python3

import sqlite3
import subprocess
import xml.etree.ElementTree as ET
from sqlite3 import Error
from Bio import Entrez
from create_db import create_connection

## Define a method to insert information into the biosample database
def insert_biosample(conn, info):
    
    sql = ''' INSERT INTO biosample(biosample_acc,taxid,strain,collected_by,collection_date,geo_loc_name,isolation_source,lat_lon,genotype,host,host_disease)
              VALUES(?,?,?,?,?,?,?,?,?,?,?) '''
    cur = conn.cursor()
    cur.execute(sql, info)
    return cur.lastrowid

def insert_sra(conn, info):
    
    sql = ''' INSERT INTO sra(srr,biosample_acc)
              VALUES(?,?) '''
    cur = conn.cursor()
    cur.execute(sql, info)
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
      print("processing biosample "+record['Id'])  

      handle_link_sra = Entrez.elink(db='sra',dbfrom='biosample',id=record['Id'])
      record_link_sra = Entrez.read(handle_link_sra)
      if len(record_link_sra[0]['LinkSetDb'])>=1:
          handle_fetch_sra = Entrez.efetch(db='sra',id=record_link_sra[0]['LinkSetDb'][0]['Link'][0]['Id'])
      else:
          print('No SRA record, skipping')
          continue
      xml_result_sra = handle_fetch_sra.read()
      root_sra = ET.fromstring(xml_result_sra)
      info = {'biosample_acc':None,'taxid':None,'strain':None,'collected_by':None,'collection_date':None,'geo_loc_name':None,'isolation_source':None,'lat_lon':None,'genotype':None,'host':None,'host_disease':None}
      info_sra = {'srr':None,'biosample_acc':None}
      taxid = root_sra.find('EXPERIMENT_PACKAGE').find('Pool').find('Member')
      if taxid.attrib['tax_id'] == '28901':
          info['taxid']=taxid.attrib['tax_id']
      else:
          print("Not Salmonella enterica, skipping")
          continue   

      srr = root_sra.find('EXPERIMENT_PACKAGE').find('RUN_SET').find('RUN')
      info_sra['srr']=srr.attrib['accession']

      attributes = root_sra.find('EXPERIMENT_PACKAGE').find('SAMPLE').find('SAMPLE_ATTRIBUTES')
      for item in attributes:
          if item.find('TAG').text == 'strain':
              info['strain']=item.find('VALUE').text
          elif item.find('TAG').text == 'collected_by':
              info['collected_by']=item.find('VALUE').text
          elif item.find('TAG').text == 'collection_date':
              info['collection_date']=item.find('VALUE').text
          elif item.find('TAG').text == 'geo_loc_name':
              info['geo_loc_name']=item.find('VALUE').text
          elif item.find('TAG').text == 'isolation_source':
              info['isolation_source']=item.find('VALUE').text
          elif item.find('TAG').text == 'lat_lon':
              info['lat_lon']=item.find('VALUE').text
          elif item.find('TAG').text == 'genotype':
              info['genotype']=item.find('VALUE').text
          elif item.find('TAG').text == 'host':
              info['host']=item.find('VALUE').text
          elif item.find('TAG').text == 'host_disease':
              info['host_disease']=item.find('VALUE').text

      handle_fetch_biosample = Entrez.efetch(db='biosample',id=record['Id'])
      xml_result_biosample = handle_fetch_biosample.read()
      root_biosample = ET.fromstring(xml_result_biosample)
      acc = root_biosample.find('BioSample')
      info['biosample_acc'] = acc.attrib['accession']
      info_sra['biosample_acc'] = acc.attrib['accession']
      for key in info:
          if info[key]==None:
              info[key]=='missing'
      conn = create_connection('mashpit.db')
      info_list=[]
      info_list = list(info.values())
      info_sra_list=[]
      info_sra_list = list(info_sra.values())
      with conn:
          insert_biosample(conn, info_list)
          insert_sra(conn,info_sra_list)  

def main():
    Entrez.email = "tongzhouxu97@gmail.com"
    Entrez.api_key = "44c9216e8c0c24c97fa7871093da74808908"
    

    conn = create_connection('mashpit.db')

    f = open("project_list.txt",'r')           
    project_list = f.readlines()
    for project in project_list:
        project = project.strip('\n')
        biosample_sra_database(project,conn)




if __name__ == '__main__':
  main()