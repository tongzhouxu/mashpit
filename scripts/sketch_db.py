#!/usr/bin/env python3

import argparse
import os
import screed
import sourmash
import threading
import pickle
import glob
import subprocess
from scripts.create_db import create_connection
from sourmash import SourmashSignature, save_signatures, load_signatures


def parse_args():
    parser = argparse.ArgumentParser(usage='sketch_db.py <database name>')
    parser.add_argument("database", help="<string>: name of the database")
    return parser.parse_args()


def signature(lk, database):
    global sketched_list
    global to_be_sketched
    global assembly_num
    if len(to_be_sketched) == 0:
        exit(0)
    SRR = to_be_sketched[0]
    try:
        subprocess.check_call(
            "dump-ref-fasta http://sra-download.ncbi.nlm.nih.gov/srapub_files/" + SRR + "_" + SRR + ".realign > " +
            "tmp/" + SRR + "_skesa.fa", shell=True, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        print("No SKESA assembly for " + SRR)
        f = open("srr_no_assembly", 'a')
        f.write(SRR+'\n')
        f.close()

    if os.stat("tmp/" + SRR + "_skesa.fa").st_size == 0:
        os.system("rm tmp/" + SRR + "_skesa.fa")
        lk.acquire()
        del to_be_sketched[0]
        lk.release()
        return
    else:
        lk.acquire()
        assembly_num = assembly_num + 1
        print("Assembly downloaded for "+SRR)
        del to_be_sketched[0]
        lk.release()

    if len(to_be_sketched) == 0 or assembly_num >= 100:
        genomes = glob.glob('tmp/*_skesa.fa')
        minhashes = []
        for g in genomes:
            mh = sourmash.MinHash(n=1000, ksize=21)
            for record in screed.open(g):
                mh.add_sequence(record.sequence, True)
            minhashes.append(mh)
        lk.acquire()
        if os.path.exists(database + '.pickle'):
            with open(database + '.pickle', 'rb') as sig_pickle:
                siglist = pickle.load(sig_pickle)
        else:
            siglist = []
        for i in range(len(minhashes)):
            siglist.append(SourmashSignature(minhashes[i], name=genomes[i].strip('tmp/').strip('_skesa.fa')))
        with open(database + '.sig', 'wt') as sig_new:
            save_signatures(siglist, sig_new)
        with open(database + '.pickle', 'wb') as sig_pickle:
            pickle.dump(siglist, sig_pickle)
        os.system("rm tmp/*")
        assembly_num = 0
        lk.release()


def main():
    args = parse_args()

    if os.path.exists(args.database + '.db'):
        pass
    else:
        print("Database does not exist. Please make sure the name is correct or run create_db.py and "
              "metadata_sra_db.py first")
        exit(0)

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

    if os.path.exists("srr_no_assembly"):
        f = open("srr_no_assembly", 'r')
        srr_list = f.readlines()
        for srr in srr_list:
            error_sra_list.append(srr.strip('\n'))
        f.close()
    else:
        os.system("touch srr_no_assembly")

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
    total_srr = []
    to_be_sketched = []
    sketched_list = []
    error_sra_list = []
    assembly_num = 0
    lock = threading.Lock()
    main()
