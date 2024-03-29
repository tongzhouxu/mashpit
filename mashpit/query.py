#!/usr/bin/env python3
import os
import logging
import ntpath
import screed
import sourmash
import heapq
import time
import pandas as pd
import matplotlib.pyplot as plt

from skbio import DistanceMatrix
from skbio.tree import nj
from Bio import Phylo
from mashpit.build import create_connection
from operator import itemgetter
from sourmash import SourmashSignature, save_signatures, load_one_signature, load_file_as_signatures

def get_query_sig(query_path, query_name,hash_number,kmer_size):
    genome = query_path
    mh = sourmash.MinHash(n=hash_number, ksize=kmer_size)
    for record in screed.open(genome):
        mh.add_sequence(record.sequence, True)
    sig = SourmashSignature(mh, name=genome)
    query_sig_path = os.path.join(os.path.dirname(query_path),f'{query_name}.sig')
    with open(query_sig_path, 'wt') as fp:
        save_signatures([sig], fp)

def generate_query_table(conn, sorted_asm_similarity_dict):
    c = conn.cursor()
    c.execute('SELECT * FROM METADATA')
    names = [description[0] for description in c.description]

    # List to collect DataFrames
    df_list = []

    for i in sorted_asm_similarity_dict:
        sql_query = pd.read_sql_query(f"SELECT * FROM METADATA WHERE asm_acc = '{i}'", conn)
        df_query = pd.DataFrame(sql_query, columns=names)
        df_query['similarity_score'] = sorted_asm_similarity_dict[i]
        df_list.append(df_query)

    # Concatenate all DataFrames at once
    output_df = pd.concat(df_list, ignore_index=True)

    # Additional processing (if it is a taxon database, add the link of the SNP cluster to the output)
    c.execute("SELECT value FROM DESC WHERE name = 'Type';")
    db_type = c.fetchone()[0]
    if db_type == 'Taxonomy':
        pds_list = output_df['PDS_acc'].to_list()
        cluster_link = [f'https://www.ncbi.nlm.nih.gov/pathogens/isolates/#{pds}' for pds in pds_list]
        output_df['link'] = cluster_link

    return output_df

def generate_mashtree(output_df,min_similarity,query_name,sig_path,added_annotation,database_sig):
    # check if the top query similarity is smaller than the threshold 
    if float(output_df['similarity_score'].iloc[0]) < min_similarity:
        logging.error('Top query similarity is smaller than the threshold. Mashtree can not be generated. ')
        exit(1)
    # select top results that are above the threshold
    acc_list = output_df[output_df['similarity_score'] >= min_similarity]['asm_acc'].to_list()
    # if the number of top results is smaller than 2, mashtree can not be generated
    if len(acc_list) < 2:
        logging.error('Number of top results is smaller than 2. Mashtree can not be generated. ')
        exit(1)
    leaves = [query_name] + acc_list
    # select database signatures that are present in acc_list
    sigs = []
    for sig in database_sig:
        if sig.name in acc_list:
            sigs.append(sig)
    # sort sigs based on the order of acc_list
    sigs = sorted(sigs, key=lambda x: acc_list.index(x.name))
    matrix = []
    # first row 
    first_row = [0]
    for similarity_score in output_df[output_df['similarity_score'] >= 0.85]['similarity_score'].to_list():
        first_row.append(1-similarity_score)
    matrix.append(first_row)
    for acc in acc_list:
        distance_list = []
        distance_list.append(1-(output_df[output_df['asm_acc'] == acc]['similarity_score'].iloc[0]))
        for sig in sigs:
            if sig.name == acc:
                query_sig = sig
        for sig in sigs:
            distance_list.append(1-(query_sig.jaccard(sig)))
        matrix.append(distance_list)

    dm = DistanceMatrix(matrix, leaves)
    newick_str = nj(dm, result_constructor=str)
    with open(f'{query_name}_tree.newick','w') as f: 
        f.write(newick_str)

    # add annotation
    if added_annotation is not None:
        annotated_leaves = []
        for leaf in leaves:
            if leaf in acc_list:
                added_annotation = str(output_df[output_df['asm_acc'] == leaf][added_annotation].iloc[0])
                leaf = leaf+' '+added_annotation
            annotated_leaves.append(leaf)
        dm = DistanceMatrix(matrix, annotated_leaves)
        newick_str = nj(dm, result_constructor=str)
        with open(f'{query_name}_tree.newick','w') as f: 
            f.write(newick_str)
    tree = Phylo.read(f'{query_name}_tree.newick', "newick")
    n = len(tree.get_terminals())
    fig = plt.figure(figsize=(10, n*0.35), dpi=300)
    axes = fig.add_subplot(1, 1, 1)
    # disable the axes and borders
    axes.set_frame_on(False)
    # remove ticks and labels
    axes.set_xticks([])
    axes.set_yticks([])
    axes.set_xticklabels([])
    axes.set_yticklabels([])
    # add white background
    fig.patch.set_facecolor('white')
    Phylo.draw(tree, axes=axes,do_show=False)
    plt.savefig(f'{query_name}_tree.png')

def query(args):
    t = time.localtime()
    current_time = time.strftime("%Y%m%d%H%M%S", t)
    print_handler = logging.StreamHandler()
    file_handler = logging.FileHandler(f"mashpit-{current_time}.log")
    print_handler.setLevel(logging.INFO)
    file_handler.setLevel(logging.DEBUG)
    logging.basicConfig(format='%(asctime)s [%(levelname)s] %(message)s', 
                        level=logging.DEBUG,
                        datefmt='%d-%b-%y %H:%M:%S', 
                        handlers=[print_handler,file_handler])
    logging.info('Start querying')
    time_start = time.time()
    number_results = args.number
    min_similarity = args.threshold
    added_annotation = args.annotation
    query_path = os.path.abspath(args.sample)
    query_name = ntpath.basename(query_path)
    
    # remove the extension if any
    if '.' in query_name:
        query_name = query_name.split('.')[0]
    db_folder = os.path.abspath(args.database)
    if not os.path.exists(db_folder):
        logging.error('Database path not found.')
        exit(1)
    folder_name = os.path.basename(db_folder)
    sql_path = os.path.join(db_folder,f'{folder_name}.db')
    sig_path = os.path.join(db_folder, f'{folder_name}.sig')
    if not (os.path.exists(sql_path) and os.path.exists(sig_path)):
        logging.error('Database incomplete.')
        exit(1)

    # check the hash number and kmer size in the database
    conn = create_connection(sql_path)
    c = conn.cursor()
    c.execute("SELECT value FROM DESC where name = 'Hash_number';")
    hash_number = int(c.fetchone()[0])
    c.execute("SELECT value FROM DESC where name = 'Kmer_size';")
    kmer_size = int(c.fetchone()[0])

    # sketch the query sample and load the signature
    get_query_sig(query_path,query_name, hash_number,kmer_size)
    query_sig_path = os.path.join(os.path.dirname(query_path),f'{query_name}.sig')
    query_sig = load_one_signature(query_sig_path)

    time_finish_sketch = time.time()
    logging.info(f'Query sample sketched in {time_finish_sketch-time_start:.2f} seconds')

    asm_similarity_dict = {}
    database_sig=list(load_file_as_signatures(sig_path))

    time_load_database_sig = time.time()   
    logging.info(f'Database signatures loaded in {time_load_database_sig-time_finish_sketch:.2f} seconds')

    for sig in database_sig:
        similarity = query_sig.jaccard(sig)
        asm_similarity_dict[str(sig)]=similarity

    time_calculate_similarity = time.time()
    logging.info(f'Jaccard similarity calculated in {time_calculate_similarity-time_load_database_sig:.2f} seconds')

    # get the top results
    top_items = heapq.nlargest(number_results, asm_similarity_dict.items(), key=itemgetter(1))
    sorted_asm_similarity_dict = dict(top_items)

    time_sort = time.time()
    logging.info(f'Top {number_results} results sorted in {time_sort-time_calculate_similarity:.2f} seconds')

    output_df = generate_query_table(conn,sorted_asm_similarity_dict)
    output_df.to_csv(query_name+'_output.csv',index=True)
    generate_mashtree(output_df,min_similarity,query_name,sig_path,added_annotation,database_sig)
    time_mashtree = time.time()
    logging.info(f'Mashtree generated in {time_mashtree-time_sort:.2f} seconds')