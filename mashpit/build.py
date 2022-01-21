#!/usr/bin/env python3

import os
import sqlite3
import subprocess
import urllib.request

import xml.etree.ElementTree as ET
import pandas as pd

from sqlite3 import Error
from Bio import Entrez
from dotenv import load_dotenv
from html.parser import HTMLParser

# connection to sqlite db file
def create_connection(db_file):
    conn = None
    try:
        conn = sqlite3.connect(db_file, check_same_thread=False)
        return conn
    except Error as e:
        print(e)

    return conn

# create table in sqlite db
def create_table(conn, create_table_sql):
    try:
        c = conn.cursor()
        c.execute(create_table_sql)
    except Error as e:
        print(e)

# create sqlite3 database with the tables
def create_database(conn):
    # metadata and accessions of a sample
    metadata_table = """CREATE TABLE IF NOT EXISTS METADATA (
                              biosample_acc    TEXT PRIMARY KEY, 
                              taxid            INTEGER,
                              strain           TEXT, 
                              collected_by     TEXT,
                              collection_date  TEXT,
                              geo_loc_name     TEXT,
                              isolation_source TEXT,
                              lat_lon          TEXT,
                              serovar          TEXT,
                              sub_species      TEXT,
                              species          TEXT,
                              genus            TEXT,
                              host             TEXT,
                              host_disease     TEXT,
                              outbreak         TEXT,
                              srr              TEXT,
                              PDT_acc          TEXT,
                              PDS_acc          TEXT
                        );"""
    # a description table to keep a record of the dabtabase type, (species and PD version if it is a standard database)
    description_table = """CREATE TABLE IF NOT EXISTS DESC (
                              name    TEXT PRIMARY KEY, 
                              value            TEXT
                        );"""
    # create tables
    create_table(conn, metadata_table)
    create_table(conn, description_table)
    conn.commit()

class PDGParser(HTMLParser):
    list = []
    def handle_data(self, data):
        if str(data).endswith('tsv'):
            self.list.append(data)

def download_from_PD(pathogen_name):
    response = urllib.request.urlopen('https://ftp.ncbi.nlm.nih.gov/pathogen/Results/'+pathogen_name+'/latest_snps/Metadata/')
    metadata_html = response.read().decode('utf-8')
    parser = PDGParser()
    parser.feed(metadata_html)
    metadata_file = parser.list[0]
    pdg_acc = metadata_file.strip('.metadata.tsv')
    parser.close()
    # TODO: progress bar for downloading the files
    urllib.request.urlretrieve('https://ftp.ncbi.nlm.nih.gov/pathogen/Results/'+pathogen_name+'/latest_snps/Metadata/'+metadata_file,metadata_file)
    urllib.request.urlretrieve('https://ftp.ncbi.nlm.nih.gov/pathogen/Results/'+pathogen_name+'/latest_snps/Clusters/'+pdg_acc+'.reference_target.SNP_distances.tsv',pdg_acc+'.reference_target.SNP_distances.tsv')
    
    # return the metadata file name
    return metadata_file

def calculate_centroid(df_raw,pdg_acc,cwd):
    # Calculate the centroid for each cluster
    subprocess.run('cat '+pdg_acc+'.reference_target.SNP_distances.tsv | cut -f1,5,9,12 > '+pdg_acc+'_selected_distance.tsv',shell=True,check=True)
    distance_file = pdg_acc+'_selected_distance.tsv'
    df_distance = pd.read_csv(distance_file,header=0,sep='\t')
    cluster_list = df_distance.groupby('PDS_acc').size().index.tolist()
    cluster_center_dict = {'PDS_acc':[],'target_acc':[]}
    # TODO: a custom folder name for the skesa assemblies maybe
    if os.path.exists(os.path.join(cwd,'fasta')):
        pass
    else:
        os.mkdir(os.path.join(cwd,'fasta'))
    for cluster in cluster_list:
        # get all isolates in this cluster
        df_cluster =  df_distance.loc[df_distance['PDS_acc'] == cluster]
        # append the dataset with switched columns to get a "full" pairwise distance matrix
        df_append = df_cluster.append(df_cluster.rename(columns={"target_acc_1":"target_acc_2","target_acc_2":"target_acc_1"}))
        # add up all distances for one target and try downloading the skesa assembly for the one with minimum distances
        for target in df_append.groupby('target_acc_1')['delta_positions_unambiguous'].sum().sort_values().index.tolist():
            SRR =  str(df_raw.loc[df_raw['target_acc'] == target]['Run'].iloc[0])
            if SRR == 'nan':
                continue
            try:
                # try downloading the skesa assembly, if failed turn to next genome in this cluster
                subprocess.check_call(
                    "dump-ref-fasta http://sra-download.ncbi.nlm.nih.gov/srapub_files/" + SRR + "_" + SRR + ".realign > fasta/" +
                    SRR + "_skeasa.fasta", shell=True, stderr=subprocess.DEVNULL)
                cluster_center_dict['PDS_acc'].append(cluster)
                cluster_center_dict['target_acc'].append(target)
                break
            except subprocess.CalledProcessError:
                continue
    # delete all the empty fasta files
    os.system('find . -size 0 -delete')
    df_cluster_center = pd.from_dict(cluster_center_dict)
    df_cluster_center.to_csv(pdg_acc+'_PDS_center.csv')

    return df_cluster_center

# build a standard database based on a Pathogen Detection metadata file
def build_standard(args,conn,cwd):
    # download the latest PD metadata files and get the metadata file name
    pathogen_name = args.species
    metadata_file = download_from_PD(pathogen_name)
    
    df_raw = pd.read_csv(metadata_file,header=0,sep='\t')
    pdg_acc = metadata_file.strip('.metadata.tsv')
    # Calculate the centroid for each cluster
    df_cluster_center = calculate_centroid(metadata_file,cwd)

    df_cluster_center_metadata = df_cluster_center['target_acc'].to_frame().join(df_raw.set_index('target_acc'),on='target_acc')
    import_metadata(df_cluster_center_metadata,conn)
    c = conn.cursor()
    c.execute("INSERT INTO DESC VALUES ('Type','Standard');")
    c.execute("INSERT INTO DESC VALUES ('Species',"+args.species+");")
    c.execute("INSERT INTO DESC VALUES ('Version',"+pdg_acc+");")
    conn.commit()


def insert_metadata(conn, info):
    sql = '''INSERT OR IGNORE INTO METADATA(
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
             species,
             genus,
             host,
             host_disease,
             outbreak,
             srr,
             PDT_acc,
             PDS_acc) 
             VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?) '''
    c = conn.cursor()
    c.execute(sql, info)
# get metadata from NCBI according to the biosample record id and the download the skesa assembly
def metadata_by_biosample_id(id, conn):
    # get the xml formatted information for the biosample
    handle_link_sra = Entrez.elink(db='sra', dbfrom='biosample', id=id)
    record_link_sra = Entrez.read(handle_link_sra)
    try:
        handle_fetch_sra = Entrez.efetch(db='sra', id=record_link_sra[0]['LinkSetDb'][0]['Link'][0]['Id'])
    except IndexError:
        raise Exception

    # get the xml formatted sra result
    xml_result_sra = handle_fetch_sra.read()
    root_sra = ET.fromstring(xml_result_sra)
    # check the sequencing library layout and source
    lib_layout = root_sra.find('EXPERIMENT_PACKAGE').find('EXPERIMENT').find('DESIGN').find('LIBRARY_DESCRIPTOR').find(
        'LIBRARY_LAYOUT')[0].tag
    lib_source = root_sra.find('EXPERIMENT_PACKAGE').find('EXPERIMENT').find('DESIGN').find('LIBRARY_DESCRIPTOR').find(
        'LIBRARY_SOURCE').text
    if lib_source != 'GENOMIC':
        raise Exception
    elif lib_layout != 'PAIRED':
        raise Exception
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
            'species': None,
            'genus': None,
            'host': None,
            'host_disease': None,
            'outbreak': None,
            'srr':None,
            'PDT_acc':None,
            'PDS_acc':None}

    taxid = root_sra.find('EXPERIMENT_PACKAGE').find('SAMPLE').find('SAMPLE_NAME').find('TAXON_ID').text
    info['taxid'] = taxid
    srr = root_sra.find('EXPERIMENT_PACKAGE').find('RUN_SET').find('RUN')
    info['srr'] = srr.attrib['accession']

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
        elif item.find('TAG').text == 'serovar':
            info['serovar'] = item.find('VALUE').text
        elif item.find('TAG').text == 'sub_species':
            info['sub_species'] = item.find('VALUE').text
        elif item.find('TAG').text == 'species':
            info['species'] = item.find('VALUE').text
        elif item.find('TAG').text == 'genus':
            info['genus'] = item.find('VALUE').text
        elif item.find('TAG').text == 'host':
            info['host'] = item.find('VALUE').text
        elif item.find('TAG').text == 'host_disease':
            info['host_disease'] = item.find('VALUE').text
        elif item.find('TAG').text == 'outbreak':
            info['outbreak'] = item.find('VALUE').text

    handle_fetch_biosample = Entrez.efetch(db='biosample', id=id)
    xml_result_biosample = handle_fetch_biosample.read()
    root_biosample = ET.fromstring(xml_result_biosample)
    acc = root_biosample.find('BioSample')
    info['biosample_acc'] = acc.attrib['accession']
    # check whether the information is missing
    for key in info:
        if info[key] is None:
            info[key] = 'missing'

    info_list = list(info.values())
    with conn:
        insert_metadata(conn, info_list)
    
    SRR = info['srr']
    if SRR == 'missing':
        # TODO: keep a log of isolates without SRR accessions
        return
    try:
        subprocess.check_call(
                    "dump-ref-fasta http://sra-download.ncbi.nlm.nih.gov/srapub_files/" + SRR + "_" + SRR + ".realign > fasta/" +
                    SRR + "_skeasa.fasta", shell=True, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        return
# import metadata to database directly from pathogen detection metadata dataframe
def import_metadata(metadata_df,conn):
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
            'species': None,
            'genus': None,
            'host': None,
            'host_disease': None,
            'outbreak': None,
            'srr':None,
            'PDT_acc':None,
            'PDS_acc':None,}

    for index, row in metadata_df.iterrows():
        info['biosample_acc'] = row['biosample_acc']
        info['taxid'] = row['taxid']
        info['strain'] = row['strain']
        info['collected_by'] = row['collected_by']
        info['collection_date'] = row['collection_date']
        info['geo_loc_name'] = row['geo_loc_name']
        info['isolation_source'] = row['isolation_source']
        info['lat_lon'] = row['lat_lon']
        info['serovar'] = row['serovar']
        info['host'] = row['host']
        info['host_disease'] = row['host_disease']
        info['outbreak'] = row['outbreak']
        info['srr'] = row['Run']
        info['PDT_acc'] = row['target_acc']
        info['PDS_acc'] = row['PDS_acc']
        # check whether the information is missing
        for key in info:
            if info[key] is None:
                info[key] = 'missing'

        info_list = list(info.values())
 
        with conn:
            insert_metadata(conn, info_list)

    return


def build(args):
    #TODO: allow a custom directory in addition to default cwd
    cwd = os.getcwd()
    db_path = os.path.join(cwd, args.name + '.db')

    #TODO: need methods to check if database already exists
    conn = create_connection(db_path)
    create_database(conn)
    
    # two types of database: standard/custom
    if args.type == 'standard':
        build_standard(args,conn,cwd)
    # Build a custom database 
    else:
        load_dotenv('.env')
        if os.environ.get('ENTREZ_EMAIL') == None:
            print("Error! Entrez email is required. Run mashpit config first.")
            exit(0)
        Entrez.email = os.environ.get('ENTREZ_EMAIL')
        if os.environ.get('ENTREZ_KEY')!= None:
            Entrez.api_key = os.environ.get('ENTREZ_KEY')
    
        if not os.path.exists("biosample.error"):
            os.system("touch biosample.error")
        cwd = os.getcwd()
        # TODO: a custom folder name for the skesa assemblies maybe
        db_path = os.path.join(cwd, args.name + '.db') 
        if os.path.exists(os.path.join(cwd,'fasta')):
            pass
        else:
            os.mkdir(os.path.join(cwd,'fasta'))
        if args.type== "biosample_list":
            f = open(args.list, 'r')
            biosample_list = f.readlines()
            for biosample in biosample_list:
                biosample = biosample.strip('\n')
                # find the ncbi-id for the specific biosample
                handle_search = Entrez.esearch(db="biosample", term=biosample)
                record_search = Entrez.read(handle_search)
                id_list = record_search['IdList']
                try:
                    metadata_by_biosample_id(id_list[0], conn)
                except:
                    f_error_log = open("biosample.error", 'a')
                    f_error_log.write(biosample + '\n')
                    f_error_log.close()
                    continue
                conn.commit()
            f.close()
            c = conn.cursor()
            c.execute("INSERT INTO DESC VALUES ('Type','List');")
            conn.commit()
        elif args.type == "keyword":
            handle_search = Entrez.esearch(db="biosample", term=args.term, retmax=100000)
            record_search = Entrez.read(handle_search)
            id_list = record_search['IdList']
            for id in id_list:
                try:
                    metadata_by_biosample_id(id, conn)
                except:
                    handle_fetch_biosample = Entrez.efetch(db='biosample', id=id)
                    xml_result_biosample = handle_fetch_biosample.read()
                    root_biosample = ET.fromstring(xml_result_biosample)
                    acc = root_biosample.find('BioSample')
                    f_error_log = open("biosample.error", 'a')
                    f_error_log.write(acc + '\n')
                    f_error_log.close()
                    continue
            conn.commit()
            c = conn.cursor()
            c.execute("INSERT INTO DESC VALUES ('Type','Keyword');")
            conn.commit()
