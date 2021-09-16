#!/usr/bin/env python3

import unittest
import subprocess
import pandas as pd
from subprocess import PIPE


class MyTests(unittest.TestCase):
    def test_script(self):
        subprocess.run('mashpit query test_sample_fasta test', shell=True)
        df = pd.read_csv('test_sample_fasta_output.csv')      
        self.assertEqual(len(df.index),3)

    def test_script_failure(self):
        result_no_args = subprocess.run(['mashpit','query'], capture_output=True)
        self.assertEqual(result_no_args.returncode, 2)


if __name__ == '__main__':
    unittest.main()