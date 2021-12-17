#!/usr/bin/env python3

import os
import ntpath
import screed
import sourmash
import multiprocessing
import pandas as pd
from mashpit.build import create_connection
from multiprocessing import Process
from operator import itemgetter
from sourmash import SourmashSignature, save_signatures, load_one_signature, load_file_as_signatures

def get_target_sig(sample_name):
    genome = sample_name
    mh = sourmash.MinHash(n=1000, ksize=31)
    for record in screed.open(genome):
        mh.add_sequence(record.sequence, True)
    sig = SourmashSignature(mh, name=genome)
    with open(sample_name + '.sig', 'wt') as fp:
        save_signatures([sig], fp)

def select_by_srr(conn, srr):
    c = conn.cursor()
    cursor = c.execute("SELECT * FROM METADATA WHERE srr=?", (srr,))
    return cursor.fetchone()[0]

def calculate_similarity(i,similarity_dict,target_sig,database):
    database_sig = load_file_as_signatures(database + '_' + str(i) + '.sig')
    for sig in database_sig:
        similarity = target_sig.jaccard(sig)
        similarity_dict[str(sig)]=similarity
    return

def query(args):
    sample_name = args.sample
    sample_path = ntpath.basename(sample_name)
    cwd = os.getcwd()
    db_path = os.path.join(cwd, args.database + '.db')
    database_sig_path = os.path.join(cwd, args.database + '.sig')
    target_sig_path = os.path.join(cwd, sample_path + '.sig')

    # check if database and signature file exists
    if os.path.exists(db_path):
        pass
    else:
        print("Database not found.")
        exit(0)
    if os.path.exists(database_sig_path):
        pass
    else:
        print("Database signature file not found")
        exit(0)


    conn = create_connection(db_path)
    c = conn.cursor()
    # sketch the query sample and load the signature
    get_target_sig(sample_path)
    target_sig = load_one_signature(target_sig_path)

    # manager dict: a shared variable for multiprocessing but slow in iteration
    manager = multiprocessing.Manager()
    srr_similarity_manager_dict = manager.dict()
    # check if the signature file has been splited (need a more elegant way)
    if os.path.exists(args.database + '_1.sig'):
        proc_list = []
        for i in range(1,args.number):
            proc = Process(target=calculate_similarity, args=(i,srr_similarity_manager_dict,target_sig,args.database))
            proc.start()
            proc_list.append(proc)
        for i in proc_list:
            i.join()
    else:
        database_sig=load_file_as_signatures(database_sig_path)
        for sig in database_sig:
            similarity = target_sig.jaccard(sig)
            srr_similarity_manager_dict[str(sig)]=similarity

    srr_similarity_dict = {}
    srr_similarity_dict.update(srr_similarity_manager_dict)
    
    # get the top 50 results
    res_srr_similarity_dict = dict(sorted(srr_similarity_dict.items(), key=itemgetter(1),reverse=True)[:50])
    c.execute('SELECT * FROM METADATA')
    output_df = pd.DataFrame([])
    names = [description[0] for description in c.description]
    for i in res_srr_similarity_dict:
        sql_query = pd.read_sql_query("select * from METADATA where srr = '" +str(i)+"'", conn)
        df_query = pd.DataFrame(sql_query,columns=names)
        df_query['similarity_score'] = res_srr_similarity_dict[i]
        output_df = output_df.append(df_query,ignore_index=True)
    
    # if it is a standard database, add the link of the snp cluster to the output
    c.execute("SELECT value FROM DESC where name = 'Type';")
    db_type = c.fetchone()[0]
    if db_type == 'Standard':
        pds_list = output_df['PDS_acc'].to_list()
        cluster_link = []
        for pds in pds_list:
            cluster_link.append('https://www.ncbi.nlm.nih.gov/pathogens/isolates/#'+pds)
        output_df['link'] = cluster_link
    print(output_df)
    output_df.to_csv(sample_name+'_output.csv',index=True)

    