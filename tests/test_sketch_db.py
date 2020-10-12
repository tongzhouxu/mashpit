#!/usr/bin/env python3

import unittest
import subprocess


class MyTests(unittest.TestCase):
    def test_script(self):
        subprocess.run('sketch_db.py test', shell=True)
        f_expected = open('tests/test.sig', 'r')
        sig_expected = f_expected.read()
        f_generated = open('test.sig', 'r')
        sig_generated = f_generated.read()
        self.assertEqual(sig_expected, sig_generated)
        f_generated.close()
        f_expected.close()

    def test_script_failure(self):
        result_no_args = subprocess.run(['sketch_db.py'], capture_output=True)
        self.assertEqual(result_no_args.returncode, 2)


if __name__ == '__main__':
    unittest.main()
