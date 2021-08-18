#!/usr/bin/env python3

from unittest.mock import MagicMock
from mashpit.create import create_connection
from subprocess import PIPE
import sys
import subprocess
import unittest
import sqlite3


class MyTests(unittest.TestCase):

    def test_sqlite3_connect_success(self):
        sqlite3.connect = MagicMock(return_value='connection succeeded')
        cc = create_connection('test.db')
        sqlite3.connect.assert_called_with('test.db', check_same_thread=False)
        self.assertEqual(cc, 'connection succeeded')

    def test_sqlite3_connect_fail(self):
        sqlite3.connect = MagicMock(return_value='connection failed')
        cc = create_connection('test.db')
        sqlite3.connect.assert_called_with('test.db', check_same_thread=False)
        self.assertEqual(cc, 'connection failed')

    def test_script(self):
        subprocess.run('mashpit create test', shell=True)
        conn = create_connection('test.db')
        c = conn.cursor()
        c.execute(''' SELECT count(name) FROM sqlite_master WHERE type='table' AND name='BIOSAMPLE' ''')
        self.assertIsNot(c.fetchone()[0], 0)
        c.execute(''' SELECT count(name) FROM sqlite_master WHERE type='table' AND name='SRA' ''')
        self.assertIsNot(c.fetchone()[0], 0)

    def test_script_failure(self):
        result = subprocess.run(['mashpit', 'create'], capture_output=True)
        self.assertEqual(result.returncode, 2)


if __name__ == '__main__':
    unittest.main()
