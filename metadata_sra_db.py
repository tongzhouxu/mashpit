#!/usr/bin/env python3

import argparse
import xml.etree.ElementTree as ET
from Bio import Entrez
from create_db import create_connection


def parse_args():
    parser = argparse.ArgumentParser(usage='metadata_sra_db.py -source <source> -list <name> -email <Entrez email '
                                           'address> -key <Entrez API key>')
    parser.add_argument("-source", help="<int>: source of data. 0: bioproject list, 1: biosample list, 2: name")
    parser.add_argument("-list", help="<string>: full name of the list file")
    parser.add_argument("-email", help="<string>: entrez email address")
    parser.add_argument("-key", help="<string>: entrez api key")
    parser.add_argument("-term", help="<string>: specie/serovar name")
    return parser.parse_args()


# Define methods to insert information into the BIOSAMPLE and SRA database
def insert_biosample(conn, info):
    sql = '''INSERT INTO BIOSAMPLE(biosample_acc,taxid,strain,collected_by,collection_date,geo_loc_name,
    isolation_source,lat_lon,genotype,host,host_disease) VALUES(?,?,?,?,?,?,?,?,?,?,?) '''
    cur = conn.cursor()
    cur.execute(sql, info)
    return cur.lastrowid


def insert_sra(conn, info):
    sql = ''' INSERT INTO SRA(srr,biosample_acc)
              VALUES(?,?) '''
    cur = conn.cursor()
    cur.execute(sql, info)
    return cur.lastrowid


def metadata_sra_by_biosample_id(id, conn):
    handle_link_sra = Entrez.elink(db='sra', dbfrom='biosample', id=id)
    record_link_sra = Entrez.read(handle_link_sra)
    # check whether the sra exists
    if len(record_link_sra[0]['LinkSetDb']) >= 1:
        handle_fetch_sra = Entrez.efetch(db='sra', id=record_link_sra[0]['LinkSetDb'][0]['Link'][0]['Id'])
    else:
        print('No SRA record, skipping it')
        return
    # get the xml formatted sra result
    xml_result_sra = handle_fetch_sra.read()
    root_sra = ET.fromstring(xml_result_sra)
    info = {'biosample_acc': None, 'taxid': None, 'strain': None, 'collected_by': None, 'collection_date': None,
            'geo_loc_name': None, 'isolation_source': None, 'lat_lon': None, 'genotype': None, 'host': None,
            'host_disease': None}
    info_sra = {'srr': None, 'biosample_acc': None}

    taxid = root_sra.find('EXPERIMENT_PACKAGE').find('SAMPLE').find('SAMPLE_NAME').find('TAXON_ID').text
    info['taxid'] = taxid
    srr = root_sra.find('EXPERIMENT_PACKAGE').find('RUN_SET').find('RUN')
    info_sra['srr'] = srr.attrib['accession']

    attributes = root_sra.find('EXPERIMENT_PACKAGE').find('SAMPLE').find('SAMPLE_ATTRIBUTES')
    for item in attributes:
        if item.find('TAG').text == 'strain':
            info['strain'] = item.find('VALUE').text
        elif item.find('TAG').text == 'collected_by':
            info['collected_by'] = item.find('VALUE').text
        elif item.find('TAG').text == 'collection_date':
            info['collection_date'] = item.find('VALUE').text
        elif item.find('TAG').text == 'geo_loc_name':
            info['geo_loc_name'] = item.find('VALUE').text
        elif item.find('TAG').text == 'isolation_source':
            info['isolation_source'] = item.find('VALUE').text
        elif item.find('TAG').text == 'lat_lon':
            info['lat_lon'] = item.find('VALUE').text
        elif item.find('TAG').text == 'genotype':
            info['genotype'] = item.find('VALUE').text
        elif item.find('TAG').text == 'host':
            info['host'] = item.find('VALUE').text
        elif item.find('TAG').text == 'host_disease':
            info['host_disease'] = item.find('VALUE').text
    # get the xml formatted information for the biosample
    handle_fetch_biosample = Entrez.efetch(db='biosample', id=id)
    xml_result_biosample = handle_fetch_biosample.read()
    root_biosample = ET.fromstring(xml_result_biosample)
    acc = root_biosample.find('BioSample')
    info['biosample_acc'] = acc.attrib['accession']
    info_sra['biosample_acc'] = acc.attrib['accession']
    # check whether the information is missing
    for key in info:
        if info[key] is None:
            info[key] = 'Missing'

    info_list = list(info.values())
    info_sra_list = list(info_sra.values())
    print(info_list)
    print(info_sra_list)
    with conn:
        insert_biosample(conn, info_list)
        insert_sra(conn, info_sra_list)


def main():
    args = parse_args()
    Entrez.email = args.email
    Entrez.api_key = args.key

    conn = create_connection('mashpit.db')

    if int(args.source) == 0:
        f = open(args.list, 'r')
        project_list = f.readlines()
        for project in project_list:
            project = project.strip('\n')
            # find the ncbi-id for the specific bioproject
            handle_search = Entrez.esearch(db="bioproject", term=project)
            record_search = Entrez.read(handle_search)
            id_list = record_search['IdList']
            # find the related biosample ids list
            handle_link = Entrez.elink(db='biosample', dbfrom='bioproject', id=id_list[0])
            record_link = Entrez.read(handle_link)
            # parse all the biosample ids
            for record in record_link[0]['LinkSetDb'][0]['Link']:
                metadata_sra_by_biosample_id(record['Id'], conn)
            conn.commit()
        f.close()
    elif int(args.source) == 1:
        f = open(args.list, 'r')
        biosample_list = f.readlines()
        for biosample in biosample_list:
            biosample = biosample.strip('\n')
            # find the ncbi-id for the specific biosample
            handle_search = Entrez.esearch(db="biosample", term=biosample)
            record_search = Entrez.read(handle_search)
            id_list = record_search['IdList']
            metadata_sra_by_biosample_id(id_list[0], conn)
            conn.commit()
        f.close()
    elif int(args.source) == 2:
        handle_search = Entrez.esearch(db="biosample", term=args.term, retmax=100000)
        record_search = Entrez.read(handle_search)
        id_list = record_search['IdList']
        for id in id_list:
            metadata_sra_by_biosample_id(id, conn)
        conn.commit()


if __name__ == '__main__':
    main()
