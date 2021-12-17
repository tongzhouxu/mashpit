#!/usr/bin/env python3
import os
from sourmash import load_file_as_signatures, save_signatures

def split(args):
    cwd = os.getcwd()
    database_sig_path = os.path.join(cwd, args.database + '.sig')

    # check the existence of the signature file
    if os.path.exists(database_sig_path):
        pass
    else:
        print("Database does not exist. Please make sure the name is correct or run mashpit build")
        exit(0)
    sig_list=[database_sig_path]
    database_sig = load_file_as_signatures(database_sig_path)
    for i in database_sig:
        sig_list.append(i)

    # the number of sigs in a single separated file
    number_sig = len(sig_list)/args.number
    for file in range(1,args.number+1):
        if file == args.number:
            with open(args.database+'_'+str(file)+'.sig','w') as f:
                save_signatures(sig_list[int((file-1)*number_sig+1):],fp=f)
        with open(args.database+'_'+str(file)+'.sig','w') as f:
            save_signatures(sig_list[int((file-1)*number_sig+1):int(file*number_sig)],fp=f)



    return