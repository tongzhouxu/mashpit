#!/usr/bin/env python3
import sqlite3
import subprocess
from sqlite3 import Error
from create_db import create_connection
from create_db import create_table

## extract all the available sketch file name and conduct a pairwise distance calculation
def main():
    conn = create_connection('mashpit.db')
    c = conn.cursor()
    cursor = c.execute("SELECT path from sketch")
    msh=[]
    for row in cursor:
        msh.append(row[0])
    for i in msh:
        for j in msh:
            res = subprocess.check_call("mash dist ./skesa_assem"+i+" ./skesa_assem"+j,shell=True)
            print(res)



if __name__ == '__main__':
  main()