import hashlib
import io
import os
import shutil
import signal
import sqlite3
import subprocess
import tarfile
import tempfile
import time
import types
import unittest
import zipfile
from pathlib import Path
from unittest.mock import Mock, patch

import pandas as pd
from sourmash import MinHash, SourmashSignature, load_file_as_signatures, save_signatures

from mashpit import build as build_module
from mashpit import gui as gui_module
from mashpit import query as query_module
from mashpit.build import (
    create_connection,
    create_database,
    download_release_files,
    download_representatives,
    fetch_accession_metadata,
    insert_accession_metadata,
    insert_metadata,
    load_metadata,
    load_tree_graph,
    resolve_release,
    safe_extract_tar,
    select_all_representatives,
    select_tree_representatives,
    sketch_assemblies,
    target_key,
    validate_pathogen_name,
)
from mashpit.mashpit import commandToArgs
from mashpit.query import generate_cluster_table, generate_mashtree, generate_query_table


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
        c.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name='METADATA';"
        )
        result = c.fetchone()
        self.assertIsNotNone(result)

    def test_description_table_exists(self):
        create_database(self.conn)
        c = self.conn.cursor()
        c.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='DESC';")
        result = c.fetchone()
        self.assertIsNotNone(result)

    def test_representative_table_exists(self):
        create_database(self.conn)
        c = self.conn.cursor()
        c.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name='REPRESENTATIVE';"
        )
        result = c.fetchone()
        self.assertIsNotNone(result)

    def test_metadata_table_structure(self):
        create_database(self.conn)
        c = self.conn.cursor()
        c.execute("PRAGMA table_info(METADATA);")
        columns = [row[1] for row in c.fetchall()]
        expected_columns = [
            "biosample_acc",
            "taxid",
            "strain",
            "collected_by",
            "collection_date",
            "geo_loc_name",
            "isolation_source",
            "lat_lon",
            "serovar",
            "sub_species",
            "species",
            "genus",
            "host",
            "host_disease",
            "outbreak",
            "srr",
            "PDT_acc",
            "PDS_acc",
            "asm_acc",
        ]
        self.assertListEqual(columns, expected_columns)

    def test_description_table_structure(self):
        create_database(self.conn)
        c = self.conn.cursor()
        c.execute("PRAGMA table_info(DESC);")
        columns = [row[1] for row in c.fetchall()]
        expected_columns = ["name", "value"]
        self.assertListEqual(columns, expected_columns)

    def test_representative_table_structure(self):
        create_database(self.conn)
        c = self.conn.cursor()
        c.execute("PRAGMA table_info(REPRESENTATIVE);")
        columns = [row[1] for row in c.fetchall()]
        expected_columns = [
            "asm_acc",
            "PDT_acc",
            "PDS_acc",
            "tree_radius",
            "selection_round",
            "fasta_path",
            "download_attempts",
        ]
        self.assertListEqual(columns, expected_columns)


class TestCommandToArgs(unittest.TestCase):
    # No test previously covered mashpit.py's own argument parsing at all -
    # notable given both real bugs found in this codebase were wiring bugs,
    # not algorithmic ones. This also regression-guards the webserver
    # subcommand removal: mashpit.py used to import a module (webserver.py)
    # that a later commit deleted, so `mashpit` crashed with
    # ModuleNotFoundError before parsing any arguments at all.
    def test_build_requires_type_and_name(self):
        with self.assertRaises(SystemExit):
            commandToArgs(["build"])

    def test_build_rejects_invalid_type(self):
        with self.assertRaises(SystemExit):
            commandToArgs(["build", "bogus", "some_db"])

    def test_build_defaults(self):
        args = commandToArgs(["build", "taxon", "some_db"])
        self.assertEqual(args.type, "taxon")
        self.assertEqual(args.name, "some_db")
        self.assertEqual(args.number, 1000)
        self.assertEqual(args.ksize, 31)
        self.assertEqual(args.radius, 20.0)
        self.assertEqual(args.download_attempts, 3)
        self.assertEqual(args.download_batch_size, 500)
        self.assertEqual(args.retry_delay, 5.0)
        self.assertEqual(args.max_reselection_rounds, 5)
        self.assertFalse(args.quiet)
        self.assertIs(args.func, build_module.build)

    def test_build_rejects_non_positive_number(self):
        with self.assertRaises(SystemExit):
            commandToArgs(["build", "taxon", "some_db", "--number", "0"])

    def test_build_rejects_negative_radius(self):
        with self.assertRaises(SystemExit):
            commandToArgs(["build", "taxon", "some_db", "--radius", "-1"])

    def test_query_defaults(self):
        args = commandToArgs(["query", "sample.fna", "some_db"])
        self.assertEqual(args.sample, "sample.fna")
        self.assertEqual(args.database, "some_db")
        self.assertEqual(args.number, 200)
        self.assertEqual(args.threshold, 0.85)
        self.assertIsNone(args.annotation)
        self.assertEqual(args.tie_tolerance_hashes, 2)
        self.assertIs(args.func, query_module.query)

    def test_webserver_subcommand_is_rejected(self):
        with self.assertRaises(SystemExit):
            commandToArgs(["webserver"])

    def test_gui_defaults(self):
        args = commandToArgs(["gui"])
        self.assertIsNone(args.port)
        self.assertIs(args.func, gui_module.gui)

    def test_gui_accepts_port(self):
        args = commandToArgs(["gui", "--port", "9999"])
        self.assertEqual(args.port, 9999)

    def test_gui_rejects_non_positive_port(self):
        with self.assertRaises(SystemExit):
            commandToArgs(["gui", "--port", "0"])


class TestSafeExtractTar(unittest.TestCase):
    # Regression test: every real NCBI SNP-tree archive (e.g.
    # SNP_trees/PDS000110997.1.tar.gz) contains a leading "." entry for the
    # archive's own top-level directory. safe_extract_tar's path-traversal
    # guard used to reject that entry, because `destination / "."` resolves
    # to `destination` itself, which does not start with
    # `str(destination) + os.sep` (no trailing path component) even though
    # it is exactly the safe extraction root. That made every real archive
    # fail with "Unsafe path in archive: .", so no SNP tree ever downloaded
    # successfully and build_taxon always raised "No representatives could
    # be selected" - only caught by actually running a real taxon build.
    def setUp(self):
        self.tmpfolder = Path("tmp")
        self.tmpfolder.mkdir()

    def tearDown(self):
        shutil.rmtree(self.tmpfolder)

    def _make_tar(self, member_names):
        buffer = io.BytesIO()
        with tarfile.open(fileobj=buffer, mode="w:gz") as archive:
            root = tarfile.TarInfo(name=".")
            root.type = tarfile.DIRTYPE
            root.mode = 0o755
            archive.addfile(root)
            for name in member_names:
                data = b"content"
                info = tarfile.TarInfo(name=name)
                info.size = len(data)
                archive.addfile(info, io.BytesIO(data))
        buffer.seek(0)
        return buffer

    def test_extracts_archive_with_root_dot_entry(self):
        buffer = self._make_tar(["./cluster.newick"])
        destination = self.tmpfolder / "extracted"
        with tarfile.open(fileobj=buffer, mode="r:gz") as archive:
            safe_extract_tar(archive, destination)
        self.assertTrue((destination / "cluster.newick").is_file())

    def test_rejects_path_traversal_member(self):
        buffer = self._make_tar(["../escaped.txt"])
        destination = self.tmpfolder / "extracted"
        with tarfile.open(fileobj=buffer, mode="r:gz") as archive:
            with self.assertRaises(RuntimeError):
                safe_extract_tar(archive, destination)


class TestDownloadReleaseFiles(unittest.TestCase):
    # Replaces the old TestDownloadMetadata: `download_metadata` no longer
    # exists. build_taxon now resolves the release first, then downloads
    # metadata, cluster, and SNP-tree files together.
    def setUp(self):
        self.pathogen_name = "Kluyvera_intermedia"
        self.pd_version = None
        self.tmpfolder = Path("tmp")
        self.tmpfolder.mkdir()

    def tearDown(self):
        shutil.rmtree(self.tmpfolder)

    def test_download_release_files(self):
        try:
            pathogen_name = validate_pathogen_name(self.pathogen_name)
            release_url, pdg_acc = resolve_release(pathogen_name, self.pd_version)
            download_release_files(release_url, pdg_acc, self.tmpfolder)
        except Exception as e:
            self.fail(f"download_release_files raised an exception: {e}")


class TestSelectTreeRepresentatives(unittest.TestCase):
    # Replaces TestCalculateCentroid: `calculate_centroid` (one representative
    # per cluster, picked from a precomputed SNP-distance matrix) was removed
    # in favor of the tree-radius method, which reads the cluster's Newick
    # SNP tree directly and adaptively selects one-or-more representatives.
    #
    # Fixture is a synthetic 5-tip caterpillar tree with known pairwise path
    # distances, so expected outputs are exact rather than NCBI ground truth:
    #
    #   (PDT0000001:10,(PDT0000002:10,(PDT0000003:10,(PDT0000004:10,PDT0000005:10):10):10):10);
    #
    # Pairwise distances: tip1-tip2=30, tip1-tip3=40, tip1-{tip4,tip5}=50,
    # tip2-tip3=30, tip2-{tip4,tip5}=40, tip3-{tip4,tip5}=30, tip4-tip5=20.
    # Diameter = 50.
    def setUp(self):
        self.tmpfolder = Path("tmp")
        self.tmpfolder.mkdir()
        self.tree_path = self.tmpfolder / "cluster.nwk"
        self.tree_path.write_text(
            "(PDT0000001.1:10,(PDT0000002.1:10,(PDT0000003.1:10,"
            "(PDT0000004.1:10,PDT0000005.1:10):10):10):10);"
        )
        self.rows = pd.DataFrame(
            {
                "asm_acc": ["GCA_1", "GCA_2", "GCA_3", "GCA_4", "GCA_5"],
                "target_key": [
                    "PDT0000001",
                    "PDT0000002",
                    "PDT0000003",
                    "PDT0000004",
                    "PDT0000005",
                ],
                "PDS_acc": ["PDSTEST"] * 5,
            }
        )

    def tearDown(self):
        shutil.rmtree(self.tmpfolder)

    def test_radius_zero_selects_every_tip(self):
        selected = select_tree_representatives(self.tree_path, self.rows, 0, set())
        self.assertEqual(
            sorted(r["asm_acc"] for r in selected),
            ["GCA_1", "GCA_2", "GCA_3", "GCA_4", "GCA_5"],
        )

    def test_radius_above_diameter_selects_one_center(self):
        selected = select_tree_representatives(self.tree_path, self.rows, 50, set())
        self.assertEqual([r["asm_acc"] for r in selected], ["GCA_2"])

    def test_intermediate_radius_is_adaptive(self):
        selected = select_tree_representatives(self.tree_path, self.rows, 25, set())
        self.assertEqual(
            sorted(r["asm_acc"] for r in selected),
            ["GCA_1", "GCA_2", "GCA_3", "GCA_5"],
        )

    def test_duplicate_target_key_is_dropped_deterministically(self):
        rows = pd.DataFrame(
            {
                "asm_acc": ["GCA_2b", "GCA_1", "GCA_2a", "GCA_3"],
                "target_key": [
                    "PDT0000002",
                    "PDT0000001",
                    "PDT0000002",
                    "PDT0000003",
                ],
                "PDS_acc": ["PDSTEST"] * 4,
            }
        )
        # GCA_2a and GCA_2b map to the same tree tip (PDT0000002); rows are
        # sorted by asm_acc before selection so the lower accession always
        # wins regardless of input row order. A single center suffices at
        # this radius, and it lands on that tip, which is why GCA_2a (never
        # GCA_2b) is the sole representative returned.
        selected = select_tree_representatives(self.tree_path, rows, 1000, set())
        selected_accs = [r["asm_acc"] for r in selected]
        self.assertNotIn("GCA_2b", selected_accs)
        self.assertEqual(selected_accs, ["GCA_2a"])

    def test_select_all_representatives_groups_by_cluster(self):
        reps, summary = select_all_representatives(
            self.rows, {"PDSTEST": self.tree_path}, 25, set(), round_number=1
        )
        self.assertEqual(
            sorted(reps["asm_acc"]), ["GCA_1", "GCA_2", "GCA_3", "GCA_5"]
        )
        self.assertEqual(summary.loc[0, "representatives"], 4)
        self.assertEqual(summary.loc[0, "status"], "complete")


class TestRealClusterRepresentativeSelection(unittest.TestCase):
    # Same selection logic as TestSelectTreeRepresentatives, but against a
    # real, tiny NCBI Pathogen Detection SNP cluster instead of a synthetic
    # tree: Listeria_innocua PDG000000091.9, cluster PDS000110997.1 (3
    # genomes). Fixtures:
    #   test_trees/PDS000110997.1.newick - the real SNP tree, as downloaded
    #     from SNP_trees/PDS000110997.1.tar.gz
    #   test_real_cluster.tsv - the matching rows from
    #     Clusters/PDG000000091.9.reference_target.cluster_list.tsv (also
    #     contains PDS000111058.1, used by
    #     TestRealClusterMultiRepresentativeSelection below)
    #
    # Tree: ('PDT000641938.1':4,'PDT000378859.2':2,'PDT000641923.1':8)'':0;
    # Pairwise distances: 378859-641938=6, 378859-641923=10, 641938-641923=12
    # (diameter=12), all confirmed by actually running select_tree_
    # representatives against this fixture rather than hand-derived.
    @classmethod
    def setUpClass(cls):
        test_dir = Path(__file__).resolve().parent
        cls.tree_path = test_dir / "test_trees" / "PDS000110997.1.newick"
        cluster = pd.read_csv(
            test_dir / "test_real_cluster.tsv", sep="\t", dtype=str
        )
        cluster = cluster[cluster["PDS_acc"] == "PDS000110997.1"]
        cls.rows = pd.DataFrame(
            {
                "asm_acc": cluster["gencoll_acc"],
                "target_key": cluster["target_acc"].map(target_key),
                "PDS_acc": cluster["PDS_acc"],
            }
        )

    def test_radius_zero_selects_every_genome(self):
        selected = select_tree_representatives(self.tree_path, self.rows, 0, set())
        self.assertEqual(
            sorted(r["asm_acc"] for r in selected),
            ["GCA_004769845.1", "GCA_010090495.1", "GCA_010091595.1"],
        )

    def test_radius_above_diameter_selects_one_genome(self):
        selected = select_tree_representatives(self.tree_path, self.rows, 15, set())
        self.assertEqual([r["asm_acc"] for r in selected], ["GCA_004769845.1"])

    def test_intermediate_radius_selects_two_genomes(self):
        selected = select_tree_representatives(self.tree_path, self.rows, 7, set())
        self.assertEqual(
            sorted(r["asm_acc"] for r in selected),
            ["GCA_004769845.1", "GCA_010091595.1"],
        )

    def test_download_and_sketch_selected_representatives(self):
        # Mimics the real build pipeline end-to-end for this one small
        # cluster - real tree, real selection, real NCBI download, real
        # sourmash sketch - without paying for a full taxon build.
        selected = select_tree_representatives(self.tree_path, self.rows, 7, set())
        accessions = sorted(r["asm_acc"] for r in selected)

        with tempfile.TemporaryDirectory(prefix="mashpit-real-cluster-") as tmp:
            tmp = Path(tmp)
            verified, failed, attempts, errors = download_representatives(
                accessions,
                tmp / "assemblies",
                3,
                500,
                5.0,
            )
            self.assertEqual(failed, set())
            representatives = pd.DataFrame({"asm_acc": accessions})
            sketch_assemblies(
                representatives,
                verified,
                tmp / "signature",
                1000,
                31,
            )
            for accession in accessions:
                sig_path = tmp / "signature" / f"{accession}.sig"
                self.assertTrue(sig_path.is_file())
                self.assertGreater(sig_path.stat().st_size, 0)


class TestRealClusterMultiRepresentativeSelection(unittest.TestCase):
    # A second real fixture, chosen to be structurally different from
    # PDS000110997.1's flat 3-tip star: PDS000111058.1 is a real,
    # genuinely nested cluster (11 genomes, one internal clade of 9 nested
    # under a 3-way root), and several sibling tips share 0-length branches
    # - both are completely ordinary in real outbreak-cluster SNP trees, and
    # neither is exercisable from a hand-built synthetic tree alone.
    #
    # Tree: (('PDT001207112.1':2,'PDT001207104.1':1,'PDT001207107.1':1,
    # 'PDT001207111.1':1,'PDT001207120.1':0,'PDT001207113.1':0,
    # 'PDT001207109.1':0,'PDT001207121.1':0,'PDT001207108.1':0)'':7,
    # 'PDT001206974.1':4,'PDT001207141.1':6)'':0;
    #
    # Selected counts below were obtained by actually running
    # select_tree_representatives against this fixture, not derived by hand.
    @classmethod
    def setUpClass(cls):
        test_dir = Path(__file__).resolve().parent
        cls.tree_path = test_dir / "test_trees" / "PDS000111058.1.newick"
        cluster = pd.read_csv(
            test_dir / "test_real_cluster.tsv", sep="\t", dtype=str
        )
        cluster = cluster[cluster["PDS_acc"] == "PDS000111058.1"]
        cls.rows = pd.DataFrame(
            {
                "asm_acc": cluster["gencoll_acc"],
                "target_key": cluster["target_acc"].map(target_key),
                "PDS_acc": cluster["PDS_acc"],
            }
        )

    def test_radius_zero_does_not_require_every_genome(self):
        # Several tips share 0-length branches to the same parent, so they
        # sit at distance 0 from each other; selecting one covers the rest
        # even at radius 0, so fewer than all 11 genomes are needed.
        selected = select_tree_representatives(self.tree_path, self.rows, 0, set())
        self.assertEqual(len(selected), 7)

    def test_intermediate_radius_selects_three_genomes(self):
        selected = select_tree_representatives(self.tree_path, self.rows, 8, set())
        self.assertEqual(
            sorted(r["asm_acc"] for r in selected),
            ["GCA_021240885.1", "GCA_021244185.1", "GCA_021251585.1"],
        )

    def test_default_radius_selects_one_genome(self):
        # 20 is build_taxon's actual default --radius.
        selected = select_tree_representatives(self.tree_path, self.rows, 20, set())
        self.assertEqual([r["asm_acc"] for r in selected], ["GCA_021251585.1"])


class TestRealMetadataPipeline(unittest.TestCase):
    # Exercises the parts of the pipeline the selection-only tests above
    # never touch: load_metadata's real-file merge of a metadata table with
    # a cluster/isolate table, and insert_metadata writing genuine NCBI
    # field values (not synthetic placeholders) into METADATA/REPRESENTATIVE.
    #
    # Fixtures: test_real_metadata.tsv (14 real rows, full metadata.tsv
    # schema, for the biosamples in both PDS000110997.1 and PDS000111058.1)
    # and test_real_cluster.tsv (the matching PDS_acc/target_acc mapping).
    @classmethod
    def setUpClass(cls):
        test_dir = Path(__file__).resolve().parent
        cls.metadata, cls.eligible = load_metadata(
            test_dir / "test_real_metadata.tsv",
            test_dir / "test_real_cluster.tsv",
        )
        cls.tree_paths = {
            "PDS000110997.1": test_dir / "test_trees" / "PDS000110997.1.newick",
            "PDS000111058.1": test_dir / "test_trees" / "PDS000111058.1.newick",
        }

    def test_load_metadata_keeps_every_fixture_row(self):
        self.assertEqual(len(self.metadata), 14)
        self.assertEqual(len(self.eligible), 14)

    def test_select_all_representatives_at_default_radius(self):
        # Both real clusters collapse to a single representative at
        # build_taxon's actual default radius of 20.
        reps, summary = select_all_representatives(
            self.eligible, self.tree_paths, 20, set(), round_number=1
        )
        self.assertEqual(
            sorted(reps["asm_acc"]), ["GCA_004769845.1", "GCA_021251585.1"]
        )
        self.assertEqual(summary["representatives"].tolist(), [1, 1])

    def test_insert_metadata_writes_real_field_values(self):
        reps, _ = select_all_representatives(
            self.eligible, self.tree_paths, 20, set(), round_number=1
        )
        conn = create_connection(":memory:")
        create_database(conn)
        verified = {acc: Path(f"/fake/{acc}.fna") for acc in reps["asm_acc"]}
        attempts = {acc: 1 for acc in reps["asm_acc"]}
        insert_metadata(conn, self.metadata, reps, 20, verified, attempts)

        cursor = conn.cursor()
        cursor.execute(
            "SELECT collection_date, geo_loc_name FROM METADATA "
            "WHERE asm_acc = 'GCA_004769845.1'"
        )
        collection_date, geo_loc_name = cursor.fetchone()
        conn.close()
        # Real values from NCBI's PDG000000091.9 metadata table, not
        # synthetic or "missing" placeholders.
        self.assertEqual(collection_date, "2018-07")
        self.assertEqual(geo_loc_name, "United Kingdom: United Kingdom")


class TestFetchAccessionMetadata(unittest.TestCase):
    # Isolated coverage of fetch_accession_metadata, decoupled from the full
    # accession-build E2E test. SAMN20822594 is the same biosample
    # TestBuildAccession builds a whole database around; its real resolved
    # assembly accession is asserted here directly instead of only being
    # exercised incidentally inside that larger test.
    def test_resolves_known_biosample(self):
        accession = fetch_accession_metadata(
            "SAMN20822594", "test@example.com", None
        )
        self.assertEqual(accession, "GCA_019647415.1")

    def test_returns_none_for_unknown_biosample(self):
        accession = fetch_accession_metadata(
            "SAMN00000000000nonexistent", "test@example.com", None
        )
        self.assertIsNone(accession)


class TestInsertAccessionMetadata(unittest.TestCase):
    # Isolated coverage of insert_accession_metadata (the fix for the bug
    # where build_accession never populated METADATA, so querying an
    # accession-built database always crashed with IndexError). No network
    # needed: verified/accession_to_biosample are just plain dicts.
    def test_writes_biosample_and_asm_acc_with_missing_placeholders(self):
        conn = create_connection(":memory:")
        create_database(conn)
        accession_to_biosample = {"GCA_000000001.1": "SAMN00000001"}
        verified = {"GCA_000000001.1": Path("/fake/GCA_000000001.1.fna")}

        insert_accession_metadata(conn, accession_to_biosample, verified)

        cursor = conn.cursor()
        cursor.execute("SELECT * FROM METADATA")
        columns = [description[0] for description in cursor.description]
        row = dict(zip(columns, cursor.fetchone()))
        conn.close()

        self.assertEqual(row["biosample_acc"], "SAMN00000001")
        self.assertEqual(row["asm_acc"], "GCA_000000001.1")
        other_columns = set(columns) - {"biosample_acc", "asm_acc"}
        for column in other_columns:
            self.assertEqual(row[column], "missing")

    def test_unmapped_accession_gets_missing_biosample(self):
        # Defensive path: an accession present in `verified` but absent from
        # accession_to_biosample (shouldn't happen in practice, since
        # build_accession populates both from the same loop, but the
        # function shouldn't KeyError if it ever does).
        conn = create_connection(":memory:")
        create_database(conn)
        verified = {"GCA_000000002.1": Path("/fake/GCA_000000002.1.fna")}

        insert_accession_metadata(conn, {}, verified)

        cursor = conn.cursor()
        cursor.execute(
            "SELECT biosample_acc FROM METADATA WHERE asm_acc = 'GCA_000000002.1'"
        )
        biosample_acc = cursor.fetchone()[0]
        conn.close()
        self.assertEqual(biosample_acc, "missing")


def make_test_signature(name, seed_offset):
    mh = MinHash(n=50, ksize=21)
    bases = "ACGT"
    sequence = "".join(bases[(i * 7 + seed_offset) % 4] for i in range(500))
    mh.add_sequence(sequence, force=True)
    return SourmashSignature(mh, name=name)


class TestGenerateQueryTable(unittest.TestCase):
    # Query's own module (query.py) previously had zero direct tests - every
    # bit of coverage was incidental to the two full E2E build+query tests
    # below, which only ever run at the default --threshold 0.85. No test
    # needs the network for this: METADATA/DESC can be populated directly.
    def _make_conn(self, db_type):
        conn = create_connection(":memory:")
        create_database(conn)
        conn.execute(
            "INSERT INTO METADATA (biosample_acc, PDS_acc, asm_acc) VALUES (?, ?, ?)",
            ("SAMN1", "PDS000000001.1", "acc1"),
        )
        conn.execute(
            "INSERT INTO METADATA (biosample_acc, PDS_acc, asm_acc) VALUES (?, ?, ?)",
            ("SAMN2", "PDS000000002.1", "acc2"),
        )
        conn.execute("INSERT INTO DESC (name, value) VALUES ('Type', ?)", (db_type,))
        conn.commit()
        return conn

    def test_taxonomy_database_adds_snp_tree_link(self):
        conn = self._make_conn("Taxonomy")
        result = generate_query_table(conn, {"acc1": 0.9, "acc2": 0.7})
        conn.close()

        self.assertEqual(result["asm_acc"].tolist(), ["acc1", "acc2"])
        self.assertEqual(result["similarity_score"].tolist(), [0.9, 0.7])
        self.assertIn("SNP_tree_link", result.columns)
        self.assertIn("PDS000000001.1", result["SNP_tree_link"].iloc[0])

    def test_accession_database_has_no_snp_tree_link(self):
        conn = self._make_conn("Accession")
        result = generate_query_table(conn, {"acc1": 0.9, "acc2": 0.7})
        conn.close()

        self.assertNotIn("SNP_tree_link", result.columns)


class TestGenerateClusterTable(unittest.TestCase):
    # Coverage for the cluster-candidate grouping: multiple representatives
    # supporting the same SNP cluster (PDS_acc) collapse into one row with a
    # hit count, instead of appearing as scattered duplicate rows in the flat
    # per-representative table. Ranking must stay on best_similarity_score -
    # representative count is context for interpreting confidence, never a
    # substitute ranking key (a cluster with more hits should not outrank a
    # cluster with a single much closer hit).
    def setUp(self):
        self.conn = create_connection(":memory:")
        create_database(self.conn)
        # PDS_A has 3 total representatives in the database; PDS_B has 1.
        for asm_acc, pds_acc in [
            ("GCA_1", "PDS_A"),
            ("GCA_2", "PDS_A"),
            ("GCA_3", "PDS_A"),
            ("GCA_4", "PDS_B"),
        ]:
            self.conn.execute(
                "INSERT INTO REPRESENTATIVE (asm_acc, PDS_acc) VALUES (?, ?)",
                (asm_acc, pds_acc),
            )
        self.conn.commit()

    def tearDown(self):
        self.conn.close()

    def test_groups_by_cluster_and_distinguishes_hits_from_total(self):
        # Only 2 of PDS_A's 3 total representatives appear in these results.
        representative_df = pd.DataFrame(
            {
                "asm_acc": ["GCA_1", "GCA_2", "GCA_4"],
                "PDS_acc": ["PDS_A", "PDS_A", "PDS_B"],
                "similarity_score": [0.99, 0.5, 0.2],
            }
        )
        cluster_df = generate_cluster_table(self.conn, representative_df, hash_number=1000)

        pds_a = cluster_df[cluster_df["PDS_acc"] == "PDS_A"].iloc[0]
        self.assertEqual(pds_a["hits_in_results"], 2)
        self.assertEqual(pds_a["total_representatives"], 3)
        self.assertAlmostEqual(pds_a["best_similarity_score"], 0.99)
        self.assertAlmostEqual(pds_a["mean_similarity_score"], (0.99 + 0.5) / 2)
        self.assertAlmostEqual(pds_a["min_similarity_score"], 0.5)

    def test_ranks_by_best_similarity_not_hit_count(self):
        # PDS_B has only 1 hit but it is the closest match overall; PDS_A has
        # 2 hits but neither is as close. best_similarity_score must win.
        representative_df = pd.DataFrame(
            {
                "asm_acc": ["GCA_1", "GCA_2", "GCA_4"],
                "PDS_acc": ["PDS_A", "PDS_A", "PDS_B"],
                "similarity_score": [0.5, 0.4, 0.99],
            }
        )
        cluster_df = generate_cluster_table(self.conn, representative_df, hash_number=1000)
        self.assertEqual(cluster_df["PDS_acc"].tolist(), ["PDS_B", "PDS_A"])

    def test_near_top_uses_sketch_resolution_not_a_fixed_cutoff(self):
        # 0.999 vs 0.998 at hash_number=1000 is a single-hash gap, within
        # tie_tolerance_hashes=2 (the default) - both should read as
        # statistically tied for best, unlike a fixed threshold such as 0.85
        # which says nothing about whether two high scores are distinguishable.
        representative_df = pd.DataFrame(
            {
                "asm_acc": ["GCA_1", "GCA_2", "GCA_4"],
                "PDS_acc": ["PDS_A", "PDS_A", "PDS_B"],
                "similarity_score": [0.999, 0.999, 0.998],
            }
        )
        cluster_df = generate_cluster_table(
            self.conn, representative_df, hash_number=1000, tie_tolerance_hashes=2
        )
        self.assertTrue(cluster_df["near_top"].all())

    def test_near_top_excludes_genuinely_distant_clusters(self):
        representative_df = pd.DataFrame(
            {
                "asm_acc": ["GCA_1", "GCA_4"],
                "PDS_acc": ["PDS_A", "PDS_B"],
                "similarity_score": [0.99, 0.5],
            }
        )
        cluster_df = generate_cluster_table(self.conn, representative_df, hash_number=1000)
        near_top_by_pds = dict(zip(cluster_df["PDS_acc"], cluster_df["near_top"]))
        self.assertTrue(near_top_by_pds["PDS_A"])
        self.assertFalse(near_top_by_pds["PDS_B"])


class TestGenerateMashtree(unittest.TestCase):
    # Regression test: generate_mashtree used a hardcoded 0.85 (instead of
    # the min_similarity/--threshold parameter) when building the query's
    # row of the distance matrix, while acc_list/leaves were sized using the
    # real threshold. Any --threshold != 0.85 that changed which hits
    # qualified made the matrix row length disagree with the leaf count,
    # crashing DistanceMatrix construction with a ValueError - dormant at
    # the default threshold, so no existing test (which only ever runs at
    # 0.85) had ever exercised it. Fixed in query.py's generate_mashtree.
    def setUp(self):
        self.query_name = "mashpit_test_mashtree_" + self.id().rsplit(".", 1)[-1]
        self.sigs = [make_test_signature(f"cand{i}", i) for i in range(4)]

    def tearDown(self):
        for path in Path.cwd().glob(f"{self.query_name}*"):
            path.unlink()

    def test_non_default_threshold_builds_tree(self):
        output_df = pd.DataFrame(
            {
                "asm_acc": ["cand0", "cand1", "cand2", "cand3"],
                "similarity_score": [0.9, 0.7, 0.6, 0.2],
            }
        )
        generate_mashtree(output_df, 0.5, self.query_name, "unused", None, self.sigs)
        self.assertTrue(Path(f"{self.query_name}_tree.newick").is_file())
        self.assertTrue(Path(f"{self.query_name}_tree.png").is_file())

    def test_below_threshold_top_hit_exits(self):
        output_df = pd.DataFrame(
            {
                "asm_acc": ["cand0", "cand1", "cand2", "cand3"],
                "similarity_score": [0.9, 0.7, 0.6, 0.2],
            }
        )
        with self.assertRaises(SystemExit):
            generate_mashtree(output_df, 0.95, self.query_name, "unused", None, self.sigs)

    def test_fewer_than_two_qualifying_hits_exits(self):
        output_df = pd.DataFrame({"asm_acc": ["cand0"], "similarity_score": [0.9]})
        with self.assertRaises(SystemExit):
            generate_mashtree(output_df, 0.5, self.query_name, "unused", None, self.sigs)


class TestNegativeBranchLength(unittest.TestCase):
    def setUp(self):
        self.tmpfolder = Path("tmp")
        self.tmpfolder.mkdir()

    def tearDown(self):
        shutil.rmtree(self.tmpfolder)

    def test_negative_branch_length_is_clamped(self):
        tree_path = self.tmpfolder / "negative.nwk"
        tree_path.write_text("(PDT0000001.1:-5,PDT0000002.1:5);")
        _, adjacency, _ = load_tree_graph(tree_path)
        lengths = [length for edges in adjacency.values() for _, length in edges]
        self.assertTrue(all(length >= 0 for length in lengths))


def hash_metadata_table(db_path):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute("PRAGMA table_info(METADATA);")
    columns = [row[1] for row in cursor.fetchall()]
    columns.sort()  # Ensure column order is deterministic
    sorted_columns = ", ".join(columns)
    cursor.execute(f"SELECT {sorted_columns} FROM METADATA ORDER BY {sorted_columns}")
    rows = cursor.fetchall()
    conn.close()

    hasher = hashlib.sha256()
    for row in rows:
        hasher.update(str(row).encode("utf-8"))
    return hasher.hexdigest()


def fake_download_factory(failing_accessions):
    # Returns a download_representatives-shaped stand-in that fails exactly
    # the given accessions on their first appearance and succeeds on any
    # later retry - used to simulate a transient NCBI download failure
    # without touching the network.
    seen = set()

    def fake_download(accessions, assembly_root, attempts, batch_size, retry_delay, api_key=None):
        verified, failed, attempt_counts, errors = {}, set(), {}, {}
        for accession in accessions:
            attempt_counts[accession] = 1
            if accession in failing_accessions and accession not in seen:
                seen.add(accession)
                failed.add(accession)
                errors[accession] = "simulated failure"
            else:
                verified[accession] = Path(f"/fake/{accession}.fna")
        return verified, failed, attempt_counts, errors

    return fake_download


def always_fail_download(accessions, assembly_root, attempts, batch_size, retry_delay, api_key=None):
    return (
        {},
        set(accessions),
        {accession: 1 for accession in accessions},
        {accession: "simulated failure" for accession in accessions},
    )


class TestBuildTaxonReselection(unittest.TestCase):
    # Regression test for the reselection-round off-by-one bug in
    # build_taxon: the retry loop used to allow only
    # (max_reselection_rounds - 1) reselections despite the flag's name and
    # help text promising max_reselection_rounds. Real NCBI download
    # failures can't be reproduced on demand, so this mocks every
    # network-touching helper build_taxon calls and lets the real
    # selection/reselection loop run against a real, tiny synthetic tree.
    #
    # Tree: (PDT0000001.1:5,PDT0000002.1:5); - two tips, distance 10 apart.
    # At radius=20 (>= diameter), select_tree_representatives deterministically
    # picks GCA_A first (confirmed by actually running it, same as
    # TestSelectTreeRepresentatives); GCA_B is the only alternate once GCA_A
    # is excluded.
    def setUp(self):
        self.name = "mashpit_test_reselection_" + self.id().rsplit(".", 1)[-1]
        self.db_folder = Path.cwd() / self.name
        self.tmp_folder = Path.cwd() / ("tmp_" + self.name)
        shutil.rmtree(self.db_folder, ignore_errors=True)
        shutil.rmtree(self.tmp_folder, ignore_errors=True)

        self.tree_path = Path.cwd() / (self.name + ".nwk")
        self.tree_path.write_text("(PDT0000001.1:5,PDT0000002.1:5);")

        self.eligible = pd.DataFrame(
            {
                "asm_acc": ["GCA_A", "GCA_B"],
                "target_key": ["PDT0000001", "PDT0000002"],
                "PDS_acc": ["PDSTEST", "PDSTEST"],
            }
        )
        self.args = types.SimpleNamespace(
            name=self.name,
            quiet=True,
            species="Test_species",
            pd_version=None,
            radius=20.0,
            download_attempts=3,
            download_batch_size=500,
            retry_delay=0.0,
            max_reselection_rounds=3,
            number=1000,
            ksize=31,
            key=None,
        )
        self.patches = [
            patch("mashpit.build.validate_pathogen_name", return_value="Test_species"),
            patch(
                "mashpit.build.resolve_release",
                return_value=("http://fake/", "PDGFAKE"),
            ),
            patch(
                "mashpit.build.download_release_files",
                return_value=(
                    Path("meta.tsv"),
                    Path("iso.tsv"),
                    {"PDSTEST": self.tree_path},
                ),
            ),
            patch(
                "mashpit.build.load_metadata",
                return_value=(self.eligible, self.eligible),
            ),
            patch(
                "mashpit.build.sketch_assemblies",
                side_effect=lambda reps, verified, sigdir, h, k: {
                    acc: Path(f"/fake/{acc}.sig") for acc in verified
                },
            ),
            patch("mashpit.build.merge_signatures", return_value=None),
        ]
        for p in self.patches:
            p.start()

    def tearDown(self):
        for p in self.patches:
            p.stop()
        shutil.rmtree(self.db_folder, ignore_errors=True)
        shutil.rmtree(self.tmp_folder, ignore_errors=True)
        self.tree_path.unlink(missing_ok=True)
        for log_file in Path.cwd().glob("*.log"):
            log_file.unlink()

    def test_reselection_picks_alternate_after_one_failure(self):
        self.args.max_reselection_rounds = 3
        fake_download = fake_download_factory({"GCA_A"})
        with patch(
            "mashpit.build.download_representatives", side_effect=fake_download
        ) as mock_download:
            build_module.build_taxon(self.args)

        # Exactly one reselection happened: initial attempt (GCA_A, fails)
        # + one retry (GCA_B, succeeds). Zero reselections (the pre-fix
        # behavior) would leave call_count at 1 and raise instead.
        self.assertEqual(mock_download.call_count, 2)

        conn = sqlite3.connect(str(self.db_folder / f"{self.name}.db"))
        cursor = conn.cursor()
        cursor.execute("SELECT asm_acc FROM REPRESENTATIVE")
        representatives = [row[0] for row in cursor.fetchall()]
        conn.close()
        self.assertEqual(representatives, ["GCA_B"])

    def test_reselection_exhausts_budget_and_raises(self):
        self.args.max_reselection_rounds = 1
        with patch(
            "mashpit.build.download_representatives", side_effect=always_fail_download
        ) as mock_download:
            with self.assertRaises(RuntimeError) as context:
                build_module.build_taxon(self.args)

        self.assertIn("1 reselection rounds", str(context.exception))
        # max_reselection_rounds=1 must permit exactly one reselection
        # attempt (2 download call-sets total) before giving up - the old
        # off-by-one bug raised after the very first failure, using zero
        # reselections despite the flag promising one.
        self.assertEqual(mock_download.call_count, 2)


class TestBuildTaxonAndQuery(unittest.TestCase):
    def setUp(self):
        self.pathogen_name = "Listeria_innocua"
        self.pd_version = "PDG000000091.9"

    def tearDown(self):
        shutil.rmtree("test_listeria_innocua")
        shutil.rmtree("ncbi_dataset")
        # remove all files starting with GCA_022617975
        for file in os.listdir():
            if file.startswith("GCA_022617975"):
                os.remove(file)
        os.remove("ncbi_dataset.zip")
        for file in os.listdir():
            if file.endswith(".log"):
                os.remove(file)

    def test_build_taxonomy(self):
        subprocess.run(
            [
                "mashpit",
                "build",
                "taxon",
                "test_listeria_innocua",
                "--species",
                self.pathogen_name,
                "--pd_version",
                self.pd_version,
            ]
        )

        # Regenerated by actually running this build end-to-end after the
        # tree-radius rewrite (98/98 clusters covered, 135 representatives
        # selected, 0 unavailable) - the old centroid-era hashes no longer
        # applied since selection changed which assemblies get included.
        expected_sqlite_sha = (
            "4cb2a693eb89d1f338b23732d1dc8dff76f2cd59835179d2154f50996adce435"
        )
        expected_signature_sha = (
            "d35c4cee32be2a0dd45346a104175f27b92b6e142a25fa7231d3d5229d81d208"
        )

        actual_sqlite_sha = hash_metadata_table(
            "test_listeria_innocua/test_listeria_innocua.db"
        )
        hasher = hashlib.sha256()
        database_sig = load_file_as_signatures(
            "test_listeria_innocua/test_listeria_innocua.sig"
        )
        database_sig_sorted = sorted(database_sig, key=lambda x: x.name)
        with open("test_listeria_innocua/test_listeria_innocua.sig.sorted", "wb") as f:
            save_signatures(database_sig_sorted, fp=f)
        with open("test_listeria_innocua/test_listeria_innocua.sig.sorted", "rb") as f:
            buf = f.read()
            hasher.update(buf)
        actual_signature_sha = hasher.hexdigest()
        self.assertEqual(actual_sqlite_sha, expected_sqlite_sha)

        self.assertEqual(actual_signature_sha, expected_signature_sha)
        subprocess.run(
            ["datasets", "download", "genome", "accession", "GCA_022617975.1"]
        )
        # unzip the downloaded file
        with zipfile.ZipFile("ncbi_dataset.zip", "r") as zip_ref:
            zip_ref.extractall("ncbi_dataset")
        subprocess.run(
            [
                "mashpit",
                "query",
                "ncbi_dataset/ncbi_dataset/data/GCA_022617975.1/GCA_022617975.1_PDT001269761.1_genomic.fna",
                "test_listeria_innocua",
            ]
        )

        # Regenerated the same way, against the representative-level query
        # result produced by this build's actual database. Renamed from
        # _output.csv when the cluster-candidate table was added below.
        expected_sha = (
            "aa6b95fa5e7df951dafe99f0c2629cba4a5c65d2fdde76ffdf811695305409a7"
        )
        hasher = hashlib.sha256()
        output_file = pd.read_csv("GCA_022617975_representative_matches.csv")
        # sort the output file by PDS_acc
        output_file.sort_values(by=["PDS_acc"], inplace=True)
        # drop the first column
        output_file.drop(output_file.columns[0], axis=1, inplace=True)
        output_file.to_csv(
            "GCA_022617975_representative_matches.csv.sorted", index=False
        )
        with open("GCA_022617975_representative_matches.csv.sorted", "rb") as f:
            buf = f.read()
            hasher.update(buf)
        actual_sha = hasher.hexdigest()
        self.assertEqual(actual_sha, expected_sha)

        # Cluster-candidate table: groups representative hits by PDS_acc so
        # multiple representatives supporting the same cluster show up as
        # one row with a hit count, ranked by best_similarity_score rather
        # than by how many representatives happened to be returned.
        expected_cluster_sha = (
            "bc1f618220b7d409ffbc7434869a45664838eda6a288a626284b5314c82b6078"
        )
        hasher = hashlib.sha256()
        cluster_file = pd.read_csv("GCA_022617975_cluster_candidates.csv")
        cluster_file.sort_values(by=["PDS_acc"], inplace=True)
        cluster_file.drop(cluster_file.columns[0], axis=1, inplace=True)
        cluster_file.to_csv(
            "GCA_022617975_cluster_candidates.csv.sorted", index=False
        )
        with open("GCA_022617975_cluster_candidates.csv.sorted", "rb") as f:
            buf = f.read()
            hasher.update(buf)
        actual_cluster_sha = hasher.hexdigest()
        self.assertEqual(actual_cluster_sha, expected_cluster_sha)


class TestBuildAccession(unittest.TestCase):
    def setUp(self):
        return

    def tearDown(self):
        shutil.rmtree("test_accession")
        os.remove("test_accession_list")
        for file in os.listdir():
            if file.endswith(".log"):
                os.remove(file)

    def test_build_accession(self):
        # generate the test accession file
        accession_list = [
            "SAMN20822594",
        ]
        with open("test_accession_list", "w") as f:
            for accession in accession_list:
                f.write(accession + "\n")
        subprocess.run(
            [
                "mashpit",
                "build",
                "accession",
                "test_accession",
                "--list",
                "test_accession_list",
                "--email",
                "test@example.com",
            ]
        )

        # build_accession never inserted anything into METADATA (only
        # build_taxon did, via insert_metadata), so METADATA was always
        # empty and querying an accession-built database crashed with
        # IndexError in generate_mashtree. Fixed by insert_accession_metadata,
        # which writes one row per verified assembly (biosample_acc + asm_acc
        # real, remaining fields "missing" since this build mode has no
        # richer metadata source). Hashes below are regenerated post-fix by
        # actually running this build end-to-end.
        expected_sqlite_sha = (
            "585a9b2cd509037829d52cc6a2ac9c04fa527856847f0270acd7138d025c283a"
        )
        expected_signature_sha = (
            "2f15d2e3e662205c4dfb1777a33ebd08cc1efadea6f798a80143d45c8aac43d8"
        )
        actual_sqlite_sha = hash_metadata_table("test_accession/test_accession.db")
        hasher = hashlib.sha256()
        with open("test_accession/test_accession.sig", "rb") as f:
            buf = f.read()
            hasher.update(buf)
        actual_signature_sha = hasher.hexdigest()
        self.assertEqual(actual_sqlite_sha, expected_sqlite_sha)
        self.assertEqual(actual_signature_sha, expected_signature_sha)


class TestMashpitGui(unittest.TestCase):
    def test_gui_starts_without_crashing(self):
        # A machine with no ~/.streamlit/credentials.toml yet blocks forever
        # on Streamlit's interactive "enter your email" onboarding prompt
        # before the server even starts - gui.suppress_first_run_prompt
        # exists specifically to avoid that hang. This only checks the
        # process survives past startup; page contents are covered manually
        # via a real browser session, similar to the old webserver smoke test
        # this replaces.
        #
        # Port is picked dynamically (bind to 0, read back the assigned
        # port) rather than hardcoded, so back-to-back runs of this test -
        # or a leftover TIME_WAIT socket from manual testing - never collide.
        import socket

        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as probe:
            probe.bind(("localhost", 0))
            free_port = probe.getsockname()[1]

        process = subprocess.Popen(
            ["mashpit", "gui", "--port", str(free_port)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        try:
            time.sleep(3)
            return_code = process.poll()
            if return_code is not None:
                stdout, stderr = process.communicate()
                self.fail(
                    f"`mashpit gui` exited early with code {return_code}:\n"
                    f"STDOUT: {stdout.decode()}\n"
                    f"STDERR: {stderr.decode()}"
                )
        finally:
            os.kill(process.pid, signal.SIGTERM)
            process.wait()


if __name__ == "__main__":
    unittest.main()
