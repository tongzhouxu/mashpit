#!/usr/bin/env python3
from pandas.testing import assert_frame_equal
import pandas as pd
import unittest
import subprocess



class MyTests(unittest.TestCase):
    def assertDataframeEqual(self, a, b, msg):
        try:
            assert_frame_equal(a, b)
        except AssertionError as e:
            raise self.failureException(msg) from e

    def setUp(self):
        self.addTypeEqualityFunc(pd.DataFrame, self.assertDataframeEqual)

    def test_script(self):
        subprocess.run('mashpit query test_sample_fasta test', shell=True)
        df_generated = pd.read_csv('test_sample_fasta_output.csv')
        df_expected  = pd.read_csv('expected_query_output.csv')
        self.assertEqual(df_generated, df_expected)

    def test_script_failure(self):
        result_no_args = subprocess.run(['mashpit','query'], capture_output=True)
        self.assertEqual(result_no_args.returncode, 2)


if __name__ == '__main__':
    unittest.main()