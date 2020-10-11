#!/usr/bin/env python3

from unittest.mock import MagicMock
from scripts.create_db import create_connection
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
        result_help = subprocess.run(['create_db.py', '-h'], capture_output=True)
        subprocess.run("create_db.py test", shell=True)
        conn = create_connection('test.db')
        c = conn.cursor()
        c.execute(''' SELECT count(name) FROM sqlite_master WHERE type='table' AND name='BIOSAMPLE' ''')
        self.assertIsNot(c.fetchone()[0], 0)
        c.execute(''' SELECT count(name) FROM sqlite_master WHERE type='table' AND name='SRA' ''')
        self.assertIsNot(c.fetchone()[0], 0)
        self.assertIn(b'usage: create_db.py <database name>\n\npositional arguments:\n  database    <string>: name of '
                      b'the database\n\noptional arguments:\n  -h, --help  show this help message and exit\n',
                      result_help.stdout)

    def test_script_failure(self):
        result = subprocess.run(['create_db.py', '--database', 'test'], capture_output=True)
        result_no_args = subprocess.run(['create_db.py'], capture_output=True)
        self.assertEqual(result.returncode, 2)
        self.assertEqual(result_no_args.returncode, 2)


if __name__ == '__main__':
    unittest.main()
