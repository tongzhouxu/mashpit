#!/usr/bin/env python3

def config(args):
    f = open('.env','w')
    f.write('ENTREZ_EMAIL='+args.email+'\n')
    if args.key is not None:
        f.write('ENTREZ_KEY='+args.key)
 
    return
