#!/usr/bin/env python3

import os
import screed
import glob
import subprocess
from mashpit.create import create_connection
from sourmash import SourmashSignature, save_signatures, load_signatures, MinHash


# method to download skesa assembly from NCBI
def download(SRR):
    cwd = os.getcwd()
    skesa_path = os.path.join(cwd,"tmp/" + SRR + "_skesa.fa")
    error_list_path = os.path.join(cwd,"srr_no_assembly")
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
        return
    else:
        print("Assembly downloaded for " + SRR)
    return


# method to sketch all the assemblies in a batch
def assemblyToSketch(database):
    cwd = os.getcwd()
    tmp_path = os.path.join(cwd, 'tmp')
    all_skesa_path = os.path.join(tmp_path,"*_skesa.fa")
    sig_path = os.path.join(cwd,database + '.sig')
    genomes = glob.glob(all_skesa_path)
    minhashes = []
    for g in genomes:
        mh = MinHash(n=1000, ksize=21)
        for record in screed.open(g):
            mh.add_sequence(record.sequence, True)
        minhashes.append(mh)
    siglist = []
    if os.path.exists(sig_path):
        sigs = load_signatures(sig_path)
        for i in sigs:
            siglist.append(i)
    for i in range(len(minhashes)):
        # the method to get the signame may need to be modified
        signame = genomes[i].strip(tmp_path).strip('_skesa.fa')
        if signame.startswith('SRR') or signame.startswith('ERR'):
            siglist.append(SourmashSignature(minhashes[i], name=signame))
    with open(sig_path, 'w') as f:
        save_signatures(siglist, fp=f)
    os.system("rm "+all_skesa_path)
    return


def sketch(args):
    cwd = os.getcwd()
    db_path = os.path.join(cwd, args.database + '.db')
    # check for the existence of the database and tables
    if os.path.exists(db_path):
        pass
    else:
        print("Database does not exist. Please make sure the name is correct or run mashpit create and mashpit metadata first.")
        exit(0)
    conn = create_connection(db_path)
    c = conn.cursor()
    c.execute(''' SELECT count(name) FROM sqlite_master WHERE type='table' AND name='BIOSAMPLE' ''')
    if c.fetchone()[0] == 0:
        print("No BIOSAMPLE table found in the database. Please make sure the name is correct or run mashpit create and mashpit metadata first.")
        exit(0)

    # create a temp folder for skesa assemblies
    tmp_path = os.path.join(cwd, 'tmp')
    if os.path.exists(tmp_path):
        pass
    else:
        os.mkdir(tmp_path)
    
    # get a list of all the sra accessions in the database
    total_srr = []
    c = conn.cursor()
    cursor = c.execute("SELECT srr from SRA")
    for row in cursor:
        total_srr.append(row[0])
    
    # check if there is already a signature file
    sketched_list = []
    sig_path = os.path.join(cwd, args.database + ".sig")
    if os.path.exists(sig_path):
        database_sig = load_signatures(sig_path)
        for sig in database_sig:
            sketched_list.append(sig.name())
    
    # check if there is a list of srr accesions that have no assembly available
    error_sra_list = []
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
    to_be_sketched = []
    for i in total_srr:
        if (i not in sketched_list) and (i not in error_sra_list):
            to_be_sketched.append(i)
    
    total = len(to_be_sketched)      
    batch = 0
    while total > 0:
        if total < args.number:
            for SRR in to_be_sketched:
                download(SRR)
            assemblyToSketch(args.database)
            exit(0)
        elif batch < args.number:
            download(to_be_sketched[0])
            del to_be_sketched[0]
            total = total - 1 
            batch = batch + 1
        elif batch >= args.number:
            batch = 0
            assemblyToSketch(args.database)

    

    conn.commit()
