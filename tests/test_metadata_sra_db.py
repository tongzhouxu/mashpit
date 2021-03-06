#!/usr/bin/env python3

import unittest
import subprocess
from scripts.create_db import create_connection
from subprocess import PIPE
import sys


class MyTests(unittest.TestCase):
    def test_script(self):
        subprocess.run('metadata_sra_db.py test biosample_list tongzhouxu97@gmail.com --list tests/test_biosample_list.txt', shell=True)
        conn = create_connection('test.db')
        c = conn.cursor()
        c.execute(''' SELECT COUNT(*) from SRA ''')
        self.assertEqual(c.fetchone()[0], 3)
        c.execute(''' SELECT COUNT(*) from BIOSAMPLE ''')
        self.assertEqual(c.fetchone()[0], 3)

    def test_script_failure(self):
        if sys.version_info[0] <= 3.7:
            result = subprocess.run(['metadata_sra_db.py', 'test', 'list', 'tongzhouxu97@gmail.com', '--list', 'tests/test_biosample_list'],  stdout=PIPE, stderr=PIPE)
        else:
            result = subprocess.run(['metadata_sra_db.py', 'test', 'list', 'tongzhouxu97@gmail.com', '--list', 'tests/test_biosample_list'], capture_output=True)
        
        if sys.version_info[0] <= 3.7:
            result_no_args = subprocess.run(['metadata_sra_db.py'], stdout=PIPE, stderr=PIPE)
        else:
            result_no_args = subprocess.run(['metadata_sra_db.py'], capture_output=True)
        self.assertEqual(result.returncode, 2)
        self.assertEqual(result_no_args.returncode, 2)


if __name__ == '__main__':
    unittest.main()