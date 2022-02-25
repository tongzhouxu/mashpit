#!/usr/bin/env python3

import unittest
import subprocess
from sourmash import load_file_as_signatures

class MyTests(unittest.TestCase):

    def test_script(self):
        subprocess.run('mashpit sketch test', shell=True)
        sig_dict_expected = {}
        sig_expected = load_file_as_signatures('expected_test.sig')
        for sig in sig_expected:
            sig_dict_expected[str(sig)]=str(sig.md5sum())

        sig_dict_generated = {}    
        sig_generated = load_file_as_signatures('test.sig')
        for sig in sig_generated:
            sig_dict_generated[str(sig)]=str(sig.md5sum())

        self.assertDictEqual(dict(sorted(sig_dict_expected.items())), dict(sorted(sig_dict_generated.items())))


    def test_script_failure(self):
        result_no_args = subprocess.run(['mashpit','sketch'], capture_output=True)
        self.assertEqual(result_no_args.returncode, 2)


if __name__ == '__main__':
    unittest.main()