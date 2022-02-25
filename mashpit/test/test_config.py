#!/usr/bin/env python3

import unittest
import subprocess
from dotenv import load_dotenv
import os


class MyTests(unittest.TestCase):
    
    def test_config(self):
        subprocess.run('mashpit config tongzhouxu97@gmail.com', shell=True)
        load_dotenv('.env')
        self.assertEqual(os.environ.get('ENTREZ_EMAIL'), "tongzhouxu97@gmail.com")

if __name__ == '__main__':
    unittest.main()