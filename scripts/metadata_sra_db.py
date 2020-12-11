#!/usr/bin/env python3

import argparse
import os
import xml.etree.ElementTree as ET
from Bio import Entrez
from scripts.create_db import create_connection


def parse_args():
    parser = argparse.ArgumentParser(usage='metadata_sra_db.py <database name> <method> <Entrez email address> '
                                           '[--key <Entrez API key>] [--list <accession list>] [--term <keyword>] ')
    parser.add_argument("database", help="<string>: name of the database")
    parser.add_argument("method", choices=["bioproject_list", "biosample_list", "keyword"], help= "<string>: data collecting method. Available options: bioproject_list, biosample_list, keyword")
    parser.add_argument("email", help="<string>: Entrez email address")
    parser.add_argument("-k", "--key", help="<string>: Entrez api key")
    parser.add_argument("-l", "--list", help="<string>: list file name of bioproject or biosample")
    parser.add_argument("-t", "--term", help="<string>: query keyword")
    return parser.parse_args()


# define methods to insert information into the BIOSAMPLE and SRA database
def insert_biosample(conn, info):
    sql = '''INSERT OR IGNORE INTO BIOSAMPLE(
             biosample_acc,
             taxid,
             strain,
             collected_by,
             collection_date,
             geo_loc_name,
             isolation_source,
             lat_lon,
             serovar,
             sub_species,
             Species,
             Genus,
             host,
             host_disease,
             outbreak) 
             VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?) '''
    c = conn.cursor()
    c.execute(sql, info)
    return c.lastrowid


def insert_sra(conn, info):
    sql = ''' INSERT OR IGNORE INTO SRA(srr,biosample_acc)
              VALUES(?,?) '''
    c = conn.cursor()
    c.execute(sql, info)
    return c.lastrowid


# define the method to get metadata from NCBI according to the biosample record id
def metadata_sra_by_biosample_id(id, conn):
    # get the xml formatted information for the biosample
    handle_link_sra = Entrez.elink(db='sra', dbfrom='biosample', id=id)
    record_link_sra = Entrez.read(handle_link_sra)
    try:
        handle_fetch_sra = Entrez.efetch(db='sra', id=record_link_sra[0]['LinkSetDb'][0]['Link'][0]['Id'])
    except IndexError:
        print("No SRA record for this BIOSAMPLE! Entrez id: " + id)
        # keep a record of the biosample ids that failed to get metadata
        f = open("biosample_error_id", 'a')
        f.write(id + '\n')
        f.close()
        return

    # get the xml formatted sra result
    xml_result_sra = handle_fetch_sra.read()
    root_sra = ET.fromstring(xml_result_sra)
    # check the sequencing library layout and source
    lib_layout = root_sra.find('EXPERIMENT_PACKAGE').find('EXPERIMENT').find('DESIGN').find('LIBRARY_DESCRIPTOR').find(
        'LIBRARY_LAYOUT')[0].tag
    lib_source = root_sra.find('EXPERIMENT_PACKAGE').find('EXPERIMENT').find('DESIGN').find('LIBRARY_DESCRIPTOR').find(
        'LIBRARY_SOURCE').text
    if lib_source != 'GENOMIC':
        print("ERROR! Library source is " + lib_source)
        return
    elif lib_layout != 'PAIRED':
        print("ERROR! Library layout is " + lib_layout)
        return
    info = {'biosample_acc': None,
            'taxid': None,
            'strain': None,
            'collected_by': None,
            'collection_date': None,
            'geo_loc_name': None,
            'isolation_source': None,
            'lat_lon': None,
            'serovar': None,
            'sub_species': None,
            'Species': None,
            'Genus': None,
            'host': None,
            'host_disease': None,
            'outbreak': None}
    info_sra = {'srr': None, 'biosample_acc': None}

    taxid = root_sra.find('EXPERIMENT_PACKAGE').find('SAMPLE').find('SAMPLE_NAME').find('TAXON_ID').text
    info['taxid'] = taxid
    srr = root_sra.find('EXPERIMENT_PACKAGE').find('RUN_SET').find('RUN')
    info_sra['srr'] = srr.attrib['accession']

    attributes = root_sra.find('EXPERIMENT_PACKAGE').find('SAMPLE').find('SAMPLE_ATTRIBUTES')
    try:
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
            elif item.find('TAG').text == 'serovar':
                info['serovar'] = item.find('VALUE').text
            elif item.find('TAG').text == 'serovar':
                info['sub_species'] = item.find('VALUE').text
            elif item.find('TAG').text == 'serovar':
                info['Species'] = item.find('VALUE').text
            elif item.find('TAG').text == 'serovar':
                info['Genus'] = item.find('VALUE').text
            elif item.find('TAG').text == 'host':
                info['host'] = item.find('VALUE').text
            elif item.find('TAG').text == 'host_disease':
                info['host_disease'] = item.find('VALUE').text
            elif item.find('TAG').text == 'outbreak':
                info['outbreak'] = item.find('VALUE').text
    except AttributeError:
        print("Error! No metadata information.")
        return

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
    with conn:
        insert_biosample(conn, info_list)
        insert_sra(conn, info_sra_list)
    print("Record stored in the database!")


def main():
    args = parse_args()
    Entrez.email = args.email
    if args.key is not None:
        Entrez.api_key = args.key
    
    if not os.path.exists("biosample_error_id"):
        os.system("touch biosample_error_id")
    cwd = os.getcwd()
    db_path = os.path.join(cwd, args.database + '.db')

    # check for the existence of the database and tables
    if os.path.exists(db_path):
        pass
    else:
        print("Database does not exist. Please make sure the name is correct or run create_db.py and "
              "metadata_sra_db.py first")
        exit(0)
    conn = create_connection(db_path)
    c = conn.cursor()
    c.execute(''' SELECT count(name) FROM sqlite_master WHERE type='table' AND name='BIOSAMPLE' ''')
    if c.fetchone()[0] == 0:
        print("No BIOSAMPLE table found in the database. Run create_db.py to create the tables first")
        exit(0)
    
    # start searching on NCBI according to the selected method
    if args.method == "bioproject_list":
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
    elif args.method == "biosample_list":
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
    elif args.method == "keyword":
        handle_search = Entrez.esearch(db="biosample", term=args.term, retmax=100000)
        record_search = Entrez.read(handle_search)
        id_list = record_search['IdList']
        for id in id_list:
            metadata_sra_by_biosample_id(id, conn)
        conn.commit()


if __name__ == '__main__':
    main()
