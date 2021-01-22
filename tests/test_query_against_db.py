#!/usr/bin/env python3

import unittest
import subprocess
from scripts.create_db import create_connection
from subprocess import PIPE
import sys


class MyTests(unittest.TestCase):
    def test_script(self):
        subprocess.run('query_against_db.py tests/test_sample_fasta test', shell=True)
        conn = create_connection('test.db')
        c = conn.cursor()
        c.execute(''' SELECT COUNT(*) from test_sample_fasta_output ''')
        self.assertEqual(c.fetchone()[0], 3)

    def test_script_failure(self):
        if sys.version_info[0] <= 3.7:
            result_no_args = subprocess.run(['query_against_db.py'], stdout=PIPE, stderr=PIPE)
        else:
            result_no_args = subprocess.run(['query_against_db.py'], capture_output=True)
        self.assertEqual(result_no_args.returncode, 2)


if __name__ == '__main__':
    unittest.main()