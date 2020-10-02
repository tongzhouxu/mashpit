#!/usr/bin/env python3

import argparse
import os
import screed
import sourmash
import threading
import pickle
from create_db import create_connection
from sourmash import SourmashSignature, save_signatures, load_signatures


def parse_args():
    parser = argparse.ArgumentParser(usage='sketch_db.py <database name>')
    parser.add_argument("database", help="<string>: name of the database")
    return parser.parse_args()


def signature(lk,database):
    global sketched_list
    global to_be_sketched
    if len(to_be_sketched) >= 1:
        SRR = to_be_sketched[0]
    else:
        exit(0)
    try:
        os.system(
            "dump-ref-fasta http://sra-download.ncbi.nlm.nih.gov/srapub_files/" + SRR + "_" + SRR + ".realign > " +
            "tmp/" + SRR + "_skesa.fa")
    except:
        print("Can't download SKESA assembly for " + SRR)
    if os.stat("tmp/" + SRR + "_skesa.fa").st_size == 0:
        os.system("rm tmp/" + SRR + "_skesa.fa")
        lk.acquire()
        del to_be_sketched[0]
        lk.release()
        return
    genome = "tmp/" + SRR + "_skesa.fa"
    minhashes = []
    mh = sourmash.MinHash(n=1000, ksize=21)
    for record in screed.open(genome):
        mh.add_sequence(record.sequence, True)
    minhashes.append(mh)
    lk.acquire()
    if os.path.exists(database + '.pickle'):
        with open(database + '.pickle', 'rb') as sig_pickle:
            siglist = pickle.load(sig_pickle)
    else:
        siglist = []
    siglist.append(SourmashSignature(minhashes[0], name=SRR))
    with open(database + '.sig', 'wt') as sig_new:
        save_signatures(siglist, sig_new)
    with open(database + '.pickle', 'wb') as sig_pickle:
        pickle.dump(siglist, sig_pickle)
    os.system("rm tmp/" + SRR + "_skesa.fa")
    print("Signature created!")
    del to_be_sketched[0]
    lk.release()


def main():
    args = parse_args()

    if os.path.exists(args.database + '.db'):
        pass
    else:
        print("Database does not exist. Please make sure the name is correct or run create_db.py and "
              "metadata_sra_db.py first")
        exit()

    conn = create_connection(args.database + '.db')

    if os.path.exists("tmp"):
        pass
    else:
        os.mkdir("tmp")

    c = conn.cursor()
    cursor = c.execute("SELECT srr from SRA")
    for row in cursor:
        total_srr.append(row[0])

    if os.path.exists(args.database + ".sig"):
        database_sig = load_signatures(args.database + ".sig")
        for sig in database_sig:
            sketched_list.append(sig.name())
    for i in total_srr:
        if i not in sketched_list:
            to_be_sketched.append(i)

    while len(to_be_sketched) >= 1:
        for i in range(6):
            t = threading.Thread(target=signature, args=(lock,args.database,))
            t.start()
            t.join()

    conn.commit()


if __name__ == '__main__':
    total_srr = []
    to_be_sketched = []
    sketched_list = []
    lock = threading.Lock()
    main()
