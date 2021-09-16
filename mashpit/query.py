#!/usr/bin/env python3

import csv
import os
import ntpath
import screed
import sourmash
import multiprocessing
import pandas as pd
from mashpit.create import create_connection
from multiprocessing import Process
from operator import itemgetter
from sourmash import SourmashSignature, save_signatures, load_one_signature, load_signatures

def get_target_sig(sample_name):
    genome = sample_name
    mh = sourmash.MinHash(n=1000, ksize=21)
    for record in screed.open(genome):
        mh.add_sequence(record.sequence, True)
    sig = SourmashSignature(mh, name=genome)
    with open(sample_name + '.sig', 'wt') as fp:
        save_signatures([sig], fp)

def select_by_srr(conn, srr):
    c = conn.cursor()
    cursor = c.execute("SELECT biosample_acc FROM SRA WHERE srr=?", (srr,))
    return cursor.fetchone()[0]

def calculate_dist(i,distance_dict,target_sig):
    database_sig = load_signatures('outbreak_' + str(i) + '.sig')
    for sig in database_sig:
        distance = target_sig.jaccard(sig)
        distance_dict[sig.name()]=distance
    return

def query(args):
    sample_name = args.sample
    sample_path = ntpath.basename(sample_name)
    cwd = os.getcwd()
    db_path = os.path.join(cwd, args.database + '.db')
    database_sig_path = os.path.join(cwd, args.database + '.sig')
    target_sig_path = os.path.join(cwd, sample_path + '.sig')

    # check if database and tables exist
    if os.path.exists(db_path):
        pass
    else:
        print("Database does not exist. Please make sure the name is correct or run mashpit create")
        exit(0)
    conn = create_connection(db_path)
    c = conn.cursor()
    c.execute(''' SELECT count(name) FROM sqlite_master WHERE type='table' AND name='BIOSAMPLE' ''')
    if c.fetchone()[0] == 0:
        print("No BIOSAMPLE table found in the database. Please make sure the name is correct or run mashpit metadata")
        exit(0)
    
    # sketch the query sample and load the signature
    get_target_sig(sample_path)
    target_sig = load_one_signature(target_sig_path)

    # manager dict: a shared variable for multiprocessing but slow in iteration
    manager = multiprocessing.Manager()
    srr_similarity_manager_dict = manager.dict()
    if os.path.exists(args.database + '_1.sig'):
        proc_list = []
        for i in range(1,args.number):
            proc = Process(target=calculate_dist, args=(i,srr_similarity_manager_dict,target_sig,))
            proc.start()
            proc_list.append(proc)
        for i in proc_list:
            i.join()
    else:
        database_sig=load_signatures(database_sig_path)
        for sig in database_sig:
            similarity = target_sig.jaccard(sig)
            srr_similarity_manager_dict[sig.name()]=similarity

    srr_similarity_dict = {}
    srr_similarity_dict.update(srr_similarity_manager_dict)

    biosample_similarity_dict = {}
    for i in srr_similarity_dict:
        biosample_acc = select_by_srr(conn,i)
        biosample_similarity_dict[biosample_acc] = srr_similarity_dict[i]
    
    # get the top 1000 results
    res_biosample_similarity_dict = dict(sorted(biosample_similarity_dict.items(), key=itemgetter(1),reverse=True)[:1000])
    c.execute('SELECT * FROM biosample')
    output_df = pd.DataFrame([])
    names = [description[0] for description in c.description]
    for i in res_biosample_similarity_dict:
        sql_query = pd.read_sql_query("select * from biosample where biosample_acc = '" +str(i)+"'", conn)
        df_query = pd.DataFrame(sql_query,columns=names)
        df_query['similarity_score'] = res_biosample_similarity_dict[i]
        output_df = output_df.append(df_query,ignore_index=True)
    print(output_df) 
    output_df.to_csv(sample_name+'_output.csv',index=True)

    