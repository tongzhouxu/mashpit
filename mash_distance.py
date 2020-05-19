#!/usr/bin/env python3
import csv
import subprocess

from create_db import create_connection


# extract all the available sketch file name and conduct a pairwise distance calculation
def main():
    conn = create_connection('mashpit.db')
    c = conn.cursor()
    cursor = c.execute("SELECT path from sketch")
    msh = []
    for row in cursor:
        msh.append(row[0])
    with open('distance.csv', mode='w', encoding='utf8') as distance:
        distance_writer = csv.writer(distance)
        for i in msh:
            for j in msh:
                res = subprocess.check_output(['mash', 'dist', i, j])
                res = res.decode('utf-8')
                line = res.split()
                distance_writer.writerow(line)

    distance.close()


if __name__ == '__main__':
    main()
