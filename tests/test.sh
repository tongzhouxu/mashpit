#!/usr/bin/env bash

../create_db.py
../metadata_sra_db.py -source 1 -list biosample_list.txt -email tongzhouxu97@gmail.com -key 44c9216e8c0c24c97fa7871093da74808908
../sketch_db.py
../mash_distance.py