#!/usr/bin/env python3

import argparse
import os
import screed
import threading
import pickle
import glob
import subprocess
from scripts.create_db import create_connection
from sourmash import SourmashSignature, save_signatures, load_signatures, MinHash


def parse_args():
    parser = argparse.ArgumentParser(usage='sketch_db.py <database name>')
    parser.add_argument("database", help="<string>: name of the database")
    return parser.parse_args()


# method to download assemblies and get signature files
def signature(lk, database):
    global assembly_num
    cwd = os.getcwd()
    if len(to_be_sketched) == 0:
        exit(0)
    SRR = to_be_sketched[0]
    
    # all the paths that are needed
    cwd = os.getcwd()
    skesa_path = os.path.join(cwd,"tmp/" + SRR + "_skesa.fa")
    tmp_path = os.path.join(cwd, 'tmp')
    error_list_path = os.path.join(cwd,"srr_no_assembly")
    all_skesa_path = os.path.join(tmp_path,"*_skesa.fa")
    pickle_path = os.path.join(cwd,database + '.pickle')
    sig_path = os.path.join(cwd,database + '.sig')
    
    # download the skesa assembly
    try:
        subprocess.check_call(
            "dump-ref-fasta http://sra-download.ncbi.nlm.nih.gov/srapub_files/" + SRR + "_" + SRR + ".realign > " +
            skesa_path, shell=True, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        print("No SKESA assembly for " + SRR)
        f = open(error_list_path, 'a')
        f.write(SRR + '\n')
        f.close()
    
    # check if the downloading is successful (there will be an empty file if failed)
    if os.stat(skesa_path).st_size == 0:
        os.system("rm " + skesa_path)
        lk.acquire()
        del to_be_sketched[0]
        lk.release()
        return
    else:
        lk.acquire()
        assembly_num = assembly_num + 1
        print("Assembly downloaded for " + SRR)
        del to_be_sketched[0]
        lk.release()
    
    # sketching files in a batch number of 1000 to save storage space
    if len(to_be_sketched) == 0 or assembly_num >= 1000:
        genomes = glob.glob(all_skesa_path)
        minhashes = []
        for g in genomes:
            mh = MinHash(n=1000, ksize=21)
            for record in screed.open(g):
                mh.add_sequence(record.sequence, True)
            minhashes.append(mh)
        lk.acquire()
        if os.path.exists(pickle_path):
            with open(pickle_path, 'rb') as sig_pickle:
                siglist = pickle.load(sig_pickle)
        else:
            siglist = []
        for i in range(len(minhashes)):
            siglist.append(SourmashSignature(minhashes[i], name=genomes[i].strip(tmp_path).strip('_skesa.fa')))
        with open(sig_path, 'wt') as sig_new:
            save_signatures(siglist, sig_new)
        with open(pickle_path, 'wb') as sig_pickle:
            pickle.dump(siglist, sig_pickle)
        os.system("rm "+all_skesa_path)
        assembly_num = 0
        lk.release()


def main():
    args = parse_args()
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
        print("No BIOSAMPLE table found in the database. Please make sure the name is correct or run create_db.py and "
              "metadata_sra_db.py first")
        exit(0)

    # create a temp folder for skesa assemblies
    tmp_path = os.path.join(cwd, 'tmp')
    if os.path.exists(tmp_path):
        pass
    else:
        os.mkdir(tmp_path)
    
    # get a list of all the sra accessions
    c = conn.cursor()
    cursor = c.execute("SELECT srr from SRA")
    for row in cursor:
        total_srr.append(row[0])
    
    # check if there is already a signature file
    sig_path = os.path.join(cwd, args.database + ".sig")
    if os.path.exists(sig_path):
        database_sig = load_signatures(sig_path)
        for sig in database_sig:
            sketched_list.append(sig.name())
    
    # check if there is a list of srr accesions that have no assembly available
    error_list_path = os.path.join(cwd,"srr_no_assembly")
    if os.path.exists(error_list_path):
        f = open(error_list_path, 'r')
        srr_list = f.readlines()
        for srr in srr_list:
            error_sra_list.append(srr.strip('\n'))
        f.close()
    else:
        os.system("touch srr_no_assembly")
    
    # create a list of srr accesions that needed to be downloaded and sketched
    for i in total_srr:
        if (i not in sketched_list) and (i not in error_sra_list):
            to_be_sketched.append(i)

    while len(to_be_sketched) >= 1:
        for i in range(6):
            t = threading.Thread(target=signature, args=(lock, args.database,))
            t.start()
            t.join()

    conn.commit()


if __name__ == '__main__':
    global sketched_list
    global to_be_sketched
    global assembly_num
    total_srr = []
    to_be_sketched = []
    sketched_list = []
    error_sra_list = []
    assembly_num = 0
    lock = threading.Lock()
    main()