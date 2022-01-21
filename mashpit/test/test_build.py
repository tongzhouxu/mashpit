#!/usr/bin/env python3

import unittest
import subprocess
from mashpit.build import create_connection


class MyTests(unittest.TestCase):
    
    def test_biosample_list(self):
        subprocess.run('mashpit config tongzhouxu97@gmail.com', shell=True)
        subprocess.run('mashpit build biosample_list test -l test_biosample_list.txt', shell=True)
        conn = create_connection('test.db')
        c = conn.cursor()
        c.execute(''' SELECT COUNT(*) from METADATA ''')
        self.assertEqual(c.fetchone()[0], 3)

    def test_biosample_failure(self):
        result_wrong_args = subprocess.run(['mashpit', 'build', 'biosample','test', '-l','test_biosample_list.txt'], capture_output=True)
        result_no_args = subprocess.run(['mashpit','build','biosample_list'], capture_output=True)
        self.assertEqual(result_wrong_args.returncode, 2)
        self.assertEqual(result_no_args.returncode, 2)


if __name__ == '__main__':
    unittest.main()