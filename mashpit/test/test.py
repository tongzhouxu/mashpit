import os
import unittest
import pandas as pd
import sqlite3
import shutil

from unittest.mock import MagicMock, patch
from mashpit.build import create_connection
from mashpit.build import create_database
from mashpit.build import download_metadata
from mashpit.build import calculate_centroid
from mashpit.build import download_and_sketch_assembly
from mashpit.update import compare_tables
from mashpit.build import build_accession

class TestCreateConnection(unittest.TestCase):
    def test_connection_is_not_none(self):
        conn = create_connection(":memory:")
        self.assertIsNotNone(conn)

    def test_connection_is_sqlite_connection(self):
        conn = create_connection(":memory:")
        self.assertIsInstance(conn, sqlite3.Connection)

class TestCreateDatabase(unittest.TestCase):
    def setUp(self):
        # Create an in-memory database and connect to it before each test
        self.conn = sqlite3.connect(":memory:")

    def tearDown(self):
        # Close the connection and destroy the database after each test
        self.conn.close()

    def test_metadata_table_exists(self):
        create_database(self.conn)
        c = self.conn.cursor()
        c.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='METADATA';")
        result = c.fetchone()
        self.assertIsNotNone(result)

    def test_description_table_exists(self):
        create_database(self.conn)
        c = self.conn.cursor()
        c.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='DESC';")
        result = c.fetchone()
        self.assertIsNotNone(result)

    def test_metadata_table_structure(self):
        create_database(self.conn)
        c = self.conn.cursor()
        c.execute("PRAGMA table_info(METADATA);")
        columns = [row[1] for row in c.fetchall()]
        expected_columns = [
            'biosample_acc', 'taxid', 'strain', 'collected_by', 'collection_date', 'geo_loc_name',
            'isolation_source', 'lat_lon', 'serovar', 'sub_species', 'species', 'genus', 'host',
            'host_disease', 'outbreak', 'srr', 'PDT_acc', 'PDS_acc', 'asm_acc'
        ]
        self.assertListEqual(columns, expected_columns)

    def test_description_table_structure(self):
        create_database(self.conn)
        c = self.conn.cursor()
        c.execute("PRAGMA table_info(DESC);")
        columns = [row[1] for row in c.fetchall()]
        expected_columns = ['name', 'value']
        self.assertListEqual(columns, expected_columns)

class TestDownloadMetadata(unittest.TestCase):
    def setUp(self):
        self.pathogen_name = 'Kluyvera_intermedia'
        self.pd_version = None
        os.mkdir('tmp')
        self.tmpfolder = 'tmp'

    def tearDown(self):
        shutil.rmtree(self.tmpfolder)
    
    def test_download_metadata(self):
        try:
            download_metadata(self.pathogen_name, self.pd_version, self.tmpfolder)
        except Exception as e:
            self.fail(f"download_metadata raised an exception: {e}")

class TestCalculateCentroid(unittest.TestCase):
    def setUp(self):
        # Define test inputs and expected output
        self.metadata_file = 'PDG000000040.52.metadata.tsv'
        self.distance_file = 'PDG000000040.52.reference_target.SNP_distances.tsv'
        self.expected_output = pd.DataFrame({'PDS_acc': ['PDS000010723.7', 'PDS000053359.1','PDS000049646.3','PDS000098715.1'], 
                                              'target_acc': ['PDT000136707.1', 'PDT000627411.1','PDT000849087.1','PDT000122974.1']})
        # sort the expected output by PDS_acc
        self.expected_output.sort_values(by=['PDS_acc'],inplace=True)
        self.tmpfolder = 'tmp'
        os.mkdir(self.tmpfolder)
        # Copy test inputs to temporary directory
        shutil.copy(self.metadata_file, self.tmpfolder)
        shutil.copy(self.distance_file, self.tmpfolder)
    

    def tearDown(self):
        # Remove temporary directory
        shutil.rmtree(self.tmpfolder)

    def test_calculate_centroid_returns_expected_output(self):
        # Call calculate_centroid with test inputs
        df_metadata = pd.read_csv(self.metadata_file, sep='\t')
        calculate_centroid(df_metadata[~df_metadata['asm_acc'].isnull()], 'PDG000000040.52', self.tmpfolder)
        # Compare actual and expected output
        actual_output = pd.read_csv(os.path.join(self.tmpfolder,'PDG000000040.52_cluster_center.tsv'),sep='\t')
        # sort the actual output by PDS_acc
        actual_output.sort_values(by=['PDS_acc'],inplace=True)
        pd.testing.assert_frame_equal(actual_output, self.expected_output)
    

class TestDownloadAndSketchAssembly(unittest.TestCase):
    def setUp(self):
        # Create a temporary directory to store downloaded and generated files
        self.hash_number = 1000
        self.kmer_size = 31
        self.gca_acc_list = ['GCA_015929625.1','GCA_009649915.1','GCA_001598315.1']
        self.tmpfolder = 'tmp'
        os.mkdir(self.tmpfolder)

    def tearDown(self):
        shutil.rmtree(self.tmpfolder)

    def test_download_and_sketch_assembly(self):
        download_and_sketch_assembly(self.gca_acc_list, self.hash_number, self.kmer_size,self.tmpfolder)
        self.assertTrue(os.path.isfile(os.path.join(self.tmpfolder,'signature', 'GCA_015929625.1.sig')))

class TestCompareTables(unittest.TestCase):

    def test_compare_tables(self):
        # Define input data
        database_cluster_center = pd.DataFrame({
            'PDS_acc_no_suffix': ['PDS_001', 'PDS_002', 'PDS_003'],
            'PDS_acc': ['acc1', 'acc2.1', 'acc3'],
            'asm_acc': ['asm1', 'asm2', 'asm3']
        })
        latest_cluster_center_metadata = pd.DataFrame({
            'PDS_acc_no_suffix': ['PDS_001', 'PDS_002', 'PDS_004'],
            'PDS_acc': ['acc1', 'acc2.2', 'acc4'],
            'asm_acc': ['asm1', 'asm2', 'asm4']
        })
        latest_cluster_center = pd.DataFrame({
            'PDS_acc_no_suffix': ['PDS_001', 'PDS_002', 'PDS_004'],
            'PDS_acc': ['acc1', 'acc2.2', 'acc4'],
            'asm_acc': ['asm1', 'asm2', 'asm4']
        })
        
        # Call the function to compare the tables
        pds_to_remove, asm_acc_remove, pds_to_add, asm_acc_add = compare_tables(database_cluster_center, latest_cluster_center_metadata)
        
        # Assert that the expected changes were identified
        self.assertEqual(pds_to_remove, ['PDS_003','PDS_002'])
        self.assertEqual(asm_acc_remove, ['asm3', 'asm2'])
        self.assertEqual(pds_to_add, ['PDS_002', 'PDS_004'])
        self.assertEqual(asm_acc_add, ['asm2', 'asm4'])

if __name__ == '__main__':
    unittest.main()