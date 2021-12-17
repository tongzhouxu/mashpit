#!/usr/bin/env python3

import os
import screed
import glob
from sourmash import SourmashSignature, save_signatures, MinHash


def sketch(args):
    cwd = os.getcwd()
    db_path = os.path.join(cwd, args.name + '.db')
    # check for the existence of the database and tables
    if os.path.exists(db_path):
        pass
    else:
        print("Database file not found. Please make sure the name is correct or run mashpit build.")
        exit(0)

    fasta_folder = os.path.join(cwd,'fasta')
    if os.path.exists(fasta_folder):
        pass
    else:
        print("Fasta folder not found.")
        exit(0)

    sig_file_name = args.name + 'sig'
    
    all_fasta_path = os.path.join(fasta_folder,"*_skeasa.fasta")
    genomes_list = glob.glob(all_fasta_path)
    minhashes = []
    for genome in genomes_list:
        mh = MinHash(n=1000,ksize=31)
        for record in screed.open(genome):
            mh.add_sequence(record.sequence, True)
        minhashes.append(mh)
    siglist = []

    for i in range(len(minhashes)):
        signame = genomes_list[i].strip(fasta_folder).strip('_skesa.fasta')
        siglist.append(SourmashSignature(minhashes[i], name=signame))
    with open(sig_file_name,'w') as f:
        save_signatures(siglist,fp=f)
