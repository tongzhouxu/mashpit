#!/usr/bin/env python3

from unittest.mock import MagicMock, Mock
from scripts.create_db import create_connection
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


if __name__ == '__main__':
    unittest.main()
