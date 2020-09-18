#!/usr/bin/env python3

import argparse
import os
import screed
import sourmash
import threading
from create_db import create_connection
from sourmash import SourmashSignature, save_signatures


def parse_args():
    parser = argparse.ArgumentParser(usage='sketch_db.py <database name>')
    parser.add_argument("database", help="<string>: name of the database")
    return parser.parse_args()


def signature(lk):
    global sketched_list
    global to_be_sketched
    SRR = to_be_sketched[0]
    try:
        os.system(
            "dump-ref-fasta http://sra-download.ncbi.nlm.nih.gov/srapub_files/" + SRR + "_" + SRR + ".realign > " +
            "database/" + SRR + "_skesa.fa")
    except :
        print("Can't download SKESA assembly for " + SRR)
    if os.stat("database/" + SRR + "_skesa.fa").st_size == 0:
        os.system("rm database/" + SRR + "_skesa.fa")
        return
    genome = "database/" + SRR + "_skesa.fa"
    minhashes = []
    mh = sourmash.MinHash(n=1000, ksize=21)
    for record in screed.open(genome):
        mh.add_sequence(record.sequence, True)
    minhashes.append(mh)
    siglist = [SourmashSignature(minhashes[0], name=SRR)]
    with open('database.sig', 'a') as fp:
        save_signatures(siglist, fp)
    os.system("rm database/" + SRR + "_skesa.fa")
    print("Signature created!")
    lk.acquire()
    sketched_list.append(SRR)
    del to_be_sketched[0]
    lk.release()
    print(threading.current_thread().getName(), sketched_list)




def main():
    args = parse_args()

    if os.path.exists(args.database + '.db'):
        pass
    else:
        print("Database does not exist. Please make sure the name is correct or run create_db.py and "
              "metadata_sra_db.py first")
        exit()

    conn = create_connection(args.database + '.db')

    if os.path.exists("database"):
        pass
    else:
        os.mkdir("database")

    c = conn.cursor()
    cursor = c.execute("SELECT srr from sra")
    for row in cursor:
        total_srr.append(row[0])

    for i in total_srr:
        c = conn.cursor()
        c.execute("SELECT SRR,Software FROM SIG WHERE SRR=?", (i,))
        if c.fetchone():
            print("Sketched. Pass.")
        else:
            to_be_sketched.append(i)

    while len(to_be_sketched) >= 1:
        for i in range(6):
            t = threading.Thread(target=signature,args=(lock,))
            t.start()
            t.join()

    for i in sketched_list:
        c = conn.cursor()
        c.execute("INSERT INTO SIG (SRR,Software) VALUES ('" + i + "','sourmash')")
    conn.commit()


if __name__ == '__main__':
    total_srr = []
    to_be_sketched = []
    sketched_list = []
    lock = threading.Lock()
    main()
