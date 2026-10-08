"""Core offline regressions. Run with: python -m pytest test/test.py"""

import hashlib
import io
import shutil
import sqlite3
import subprocess
import tarfile
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import patch
import pandas as pd
from sourmash import MinHash, SourmashSignature
from mashpit import build as build_module
from mashpit.build import (
    create_connection,
    create_database,
    find_local_fasta_files,
    insert_accession_metadata,
    insert_metadata,
    load_all_signatures,
    load_metadata,
    load_tree_graph,
    read_annotation_log,
    reshard_database,
    safe_extract_tar,
    select_all_representatives,
    select_tree_representatives,
    sketch_custom_sequences,
    validate_column_name,
    write_signature_shards,
    write_signatures,
)
from mashpit.query import (
    MashtreeSkipped,
    compare_sharded_database,
    generate_cluster_table,
    generate_mashtree,
    generate_query_table,
)
import sys
import pytest


@pytest.fixture(autouse=True)
def isolated_workspace(tmp_path, monkeypatch):
    """Keep generated databases, trees, and logs out of the checkout."""
    monkeypatch.chdir(tmp_path)


class TestSafeExtractTar(unittest.TestCase):
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


class TestSelectTreeRepresentatives(unittest.TestCase):
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

    def test_selection_across_radii(self):
        for radius, expected in [
            (0, ["GCA_1", "GCA_2", "GCA_3", "GCA_4", "GCA_5"]),
            (25, ["GCA_1", "GCA_2", "GCA_3", "GCA_5"]),
            (50, ["GCA_2"]),
        ]:
            with self.subTest(radius=radius):
                selected = select_tree_representatives(
                    self.tree_path, self.rows, radius, set()
                )
                self.assertEqual(sorted(r["asm_acc"] for r in selected), expected)

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
        selected = select_tree_representatives(self.tree_path, rows, 1000, set())
        selected_accs = [r["asm_acc"] for r in selected]
        self.assertNotIn("GCA_2b", selected_accs)
        self.assertEqual(selected_accs, ["GCA_2a"])


class TestRealMetadataPipeline(unittest.TestCase):
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
        self.assertEqual(collection_date, "2018-07")
        self.assertEqual(geo_loc_name, "United Kingdom: United Kingdom")


class TestInsertAccessionMetadata(unittest.TestCase):
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


def make_test_signature(name, seed_offset):
    mh = MinHash(n=50, ksize=21)
    bases = "ACGT"
    sequence = "".join(bases[(i * 7 + seed_offset) % 4] for i in range(500))
    mh.add_sequence(sequence, force=True)
    return SourmashSignature(mh, name=name)


class TestGenerateQueryTable(unittest.TestCase):
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


class TestGenerateMashtree(unittest.TestCase):
    def setUp(self):
        self.query_name = "mashpit_test_mashtree_" + self.id().rsplit(".", 1)[-1]
        self.sigs = [make_test_signature(f"cand{i}", i) for i in range(4)]

    def tearDown(self):
        for path in Path.cwd().glob(f"{self.query_name}*"):
            path.unlink()

    def test_tree_skipped_for_insufficient_hits(self):
        output_df = pd.DataFrame(
            {"asm_acc": ["cand0", "cand1"], "similarity_score": [0.9, 0.7]}
        )
        for threshold in (0.95, 0.8):
            with self.subTest(threshold=threshold), self.assertRaises(MashtreeSkipped):
                generate_mashtree(
                    output_df, threshold, self.query_name, "unused", None, self.sigs
                )

    def test_annotation_with_multiple_hits(self):
        output_df = pd.DataFrame(
            {
                "asm_acc": ["cand0", "cand1", "cand2", "cand3"],
                "similarity_score": [0.9, 0.7, 0.6, 0.2],
                "isolation_source": ["water", "soil", "clinical", "food"],
            }
        )
        generate_mashtree(
            output_df, 0.5, self.query_name, "unused", "isolation_source", self.sigs
        )
        newick_text = Path(f"{self.query_name}_tree.newick").read_text()
        self.assertIn("cand0_water", newick_text)
        self.assertIn("cand1_soil", newick_text)
        self.assertIn("cand2_clinical", newick_text)
        self.assertNotIn("cand3", newick_text)
        for suffix in ("png", "svg"):
            self.assertGreater(
                Path(f"{self.query_name}_tree.{suffix}").stat().st_size, 0
            )


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


def fake_download_factory(failing_accessions):
    seen = set()

    def fake_download(
        accessions, assembly_root, attempts, batch_size, retry_delay, api_key=None
    ):
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


def always_fail_download(
    accessions, assembly_root, attempts, batch_size, retry_delay, api_key=None
):
    return (
        {},
        set(accessions),
        {accession: 1 for accession in accessions},
        {accession: "simulated failure" for accession in accessions},
    )


class TestBuildTaxonReselection(unittest.TestCase):
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
        prepare = build_module.prepare
        self.patches = [
            patch(
                "mashpit.build.prepare",
                side_effect=lambda args: prepare(args, require_datasets=False),
            ),
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
        self.assertEqual(mock_download.call_count, 2)


class TestFindLocalFastaFiles(unittest.TestCase):
    def test_rejects_duplicate_sample_id_across_extensions(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "sample1.fasta").write_text(">c1\nACGT\n")
            (tmp / "sample1.fa").write_text(">c1\nACGT\n")
            with self.assertRaises(ValueError):
                find_local_fasta_files(tmp)


class TestSketchCustomSequences(unittest.TestCase):
    def test_one_malformed_file_does_not_prevent_others_from_sketching(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            good = tmp / "good_sample.fasta"
            good.write_text(">c1\n" + "ACGTACGTAC" * 10 + "\n")
            empty = tmp / "empty_sample.fasta"
            empty.write_text("")
            fasta_paths = {"good_sample": good, "empty_sample": empty}
            signature_paths, errors = sketch_custom_sequences(
                fasta_paths, tmp / "sigs", 100, 21
            )
            self.assertEqual(set(signature_paths), {"good_sample"})
            self.assertIn("empty_sample", errors)


class TestReshardDatabase(unittest.TestCase):
    def setUp(self):
        self.db_folder = Path(tempfile.mkdtemp(prefix="mashpit_reshard_test_"))
        self.db_name = "reshardtest"
        conn = create_connection(str(self.db_folder / (self.db_name + ".db")))
        create_database(conn)
        conn.executemany(
            "INSERT OR REPLACE INTO DESC(name, value) VALUES (?, ?)",
            [("Type", "Custom"), ("Hash_number", "50"), ("Kmer_size", "21")],
        )
        conn.commit()
        conn.close()

        self.signatures = [make_test_signature(f"sig{i}", i) for i in range(6)]
        self.sig_path = self.db_folder / (self.db_name + ".sig")
        write_signatures(self.signatures, self.sig_path, nshards=1)

    def tearDown(self):
        shutil.rmtree(self.db_folder, ignore_errors=True)

    def test_original_survives_an_interruption_installing_the_replacement(self):
        original_rename = Path.rename

        def failing_rename(self_path, target):
            if self_path.name.endswith(".reshard-tmp"):
                raise OSError("simulated failure installing replacement")
            return original_rename(self_path, target)

        with patch.object(Path, "rename", failing_rename):
            with self.assertRaises(OSError):
                reshard_database(self.db_folder, nshards=3)

        self.assertTrue(self.sig_path.is_file())
        loaded_names = {str(sig) for sig in load_all_signatures(self.sig_path)}
        self.assertEqual(loaded_names, {str(sig) for sig in self.signatures})


class TestCompareShardedDatabase(unittest.TestCase):
    def setUp(self):
        self.sig_dir = Path(tempfile.mkdtemp(prefix="mashpit_compare_shard_test_"))
        self.signatures = [make_test_signature(f"sig{i:02d}", i) for i in range(12)]
        write_signature_shards(self.signatures, self.sig_dir, nshards=4)

        self.query_dir = Path(tempfile.mkdtemp(prefix="mashpit_compare_shard_query_"))
        self.query_sig = make_test_signature("query", seed_offset=1)
        self.query_sig_path = self.query_dir / "query.sig"
        with self.query_sig_path.open("wt") as handle:
            from sourmash import save_signatures

            save_signatures([self.query_sig], fp=handle)

    def tearDown(self):
        shutil.rmtree(self.sig_dir, ignore_errors=True)
        shutil.rmtree(self.query_dir, ignore_errors=True)

    def test_matches_single_threaded_jaccard_over_the_same_signatures(self):
        expected = {str(sig): self.query_sig.jaccard(sig) for sig in self.signatures}

        result, nworkers = compare_sharded_database(
            self.sig_dir, str(self.query_sig_path), top_n=12, max_workers=2
        )

        self.assertEqual(nworkers, 2)
        self.assertEqual(set(result), set(expected))
        for name, similarity in expected.items():
            self.assertAlmostEqual(result[name], similarity)


class TestClusterReport(unittest.TestCase):
    def test_dates_and_conflicting_serotypes(self):
        from mashpit.report import collection_year, computed_serotype, source_group

        self.assertEqual(collection_year("2020-02-29"), 2020)
        for value in ("2020-02-31", "2018/2020", "missing", None):
            self.assertIsNone(collection_year(value))
        result = computed_serotype("serotype=Enteritidis; serovar=Typhimurium; MLST=11")
        self.assertIn("Enteritidis", result)
        self.assertIn("Typhimurium", result)
        self.assertNotIn("MLST", result)
        self.assertEqual(
            computed_serotype("serotype=4,[5],12:i:-; MLST=19"),
            "serotype: 4,[5],12:i:-",
        )
        self.assertEqual(computed_serotype("serotype=missing"), "Unknown")
        self.assertEqual(source_group("missing"), "Unknown / other")
        self.assertEqual(source_group("nonclinical"), "Unknown / other")

    def test_tree_height_and_exports(self):
        from mashpit.tree_plot import render_tree
        from io import BytesIO
        from PIL import Image

        png, svg, tips = render_tree("(query:0.1,(a:0.1,b:0.2):0.1);", "query")
        self.assertEqual(tips, 3)
        self.assertIn(b"<svg", svg)
        width, height = Image.open(BytesIO(png)).size
        self.assertGreater(width, height)
        many = "(" + ",".join(f"tip_{i}:.1" for i in range(100)) + ");"
        large_png, _, count = render_tree(many)
        self.assertEqual(count, 100)
        self.assertGreater(Image.open(BytesIO(large_png)).height, height * 5)

    def test_streamlit_report_and_controls(self):
        from streamlit.testing.v1 import AppTest

        script = """
import pandas as pd
import streamlit as st
from mashpit.report import summarize_clusters
from mashpit.streamlit_app import display_results
st.set_page_config(layout="wide")
members = pd.DataFrame({"target_acc": ["PDT1", "PDT2"], "PDS_acc": ["PDS1", "PDS1"],
    "computed_types": ["serotype=Enteritidis", None], "collection_date": ["2020", "2021-02"],
    "epi_type": ["clinical", "environmental/other"]})
clusters = summarize_clusters(members)
clusters["best_similarity_score"] = .99
clusters["near_top"] = True
representatives = pd.DataFrame({"asm_acc": ["a", "b"], "similarity_score": [.99, .98]})
results = {"query_name": "query", "representative_df": representatives,
    "cluster_df": clusters, "members": members, "cluster_csv": b"csv", "representative_csv": b"csv",
    "tree_newick": b"(query:0.1,(a:0.1,b:0.2):0.1);", "log": ""}
display_results(results, {"Type": "Taxonomy"})
"""
        app = AppTest.from_string(script).run(timeout=30)
        self.assertEqual(len(app.exception), 0, str(app.exception))
        self.assertEqual(app.metric[3].value, "2")
        app.slider[0].set_value(14).run(timeout=30)
        self.assertEqual(len(app.exception), 0, str(app.exception))
        legacy = script.replace(
            '"members": members', '"members": pd.DataFrame()'
        ).replace(
            "clusters = summarize_clusters(members)",
            'clusters = pd.DataFrame({"PDS_acc": ["PDS1"]})',
        )
        app = AppTest.from_string(legacy).run(timeout=30)
        self.assertEqual(len(app.exception), 0, str(app.exception))
        custom = script.replace('"cluster_df": clusters', '"cluster_df": None')
        app = AppTest.from_string(custom).run(timeout=30)
        self.assertEqual(len(app.exception), 0, str(app.exception))


def test_cluster_membership_summary_and_ranking(tmp_path):
    """Join full membership, persist it, then rank by similarity rather than size."""
    from mashpit.report import read_cluster_members

    metadata = tmp_path / "metadata.tsv"
    membership = tmp_path / "membership.tsv"
    metadata.write_text(
        "target_acc\tasm_acc\tbiosample_acc\tepi_type\tcollection_date\tcomputed_types\n"
        "PDT000000001.1\tGCA_000000001.1\tSAMN1\tclinical\t2019\tserotype=4,[5],12:i:-\n"
        "PDT000000002.1\tGCA_000000002.1\tSAMN2\tenvironmental/other\t2021-02\tserotype=Enteritidis\n"
        "PDT000000003.1\tGCA_000000003.1\tSAMN3\tclinical\t2020\tserotype=Newport\n"
    )
    membership.write_text(
        "target_acc\tPDS_acc\nPDT000000001.1\tA\nPDT000000002.1\tA\n"
        "PDT000000004.1\tA\nPDT000000004.1\tA\nPDT000000003.1\tB\n"
    )
    full, eligible = load_metadata(metadata, membership)
    assert len(full) == 4 and len(eligible) == 3
    conn = create_connection(":memory:")
    try:
        create_database(conn)
        hits = pd.DataFrame(
            {"PDS_acc": ["A", "A", "B"], "similarity_score": [0.990, 0.980, 0.999]}
        )
        # An old database must not substitute representatives for cluster size.
        assert read_cluster_members(conn, ["A"]).empty
        assert generate_cluster_table(conn, hits, 1000)["cluster_size"].isna().all()
        insert_metadata(
            conn,
            full,
            eligible,
            20,
            {acc: Path("fake") for acc in eligible.asm_acc},
            {},
        )
        assert (
            conn.execute(
                "SELECT computed_types FROM METADATA WHERE biosample_acc='SAMN1'"
            ).fetchone()[0]
            == "serotype=4,[5],12:i:-"
        )
        summary = generate_cluster_table(conn, hits, 1000)
        assert summary.PDS_acc.tolist() == ["B", "A"]
        assert summary.near_top.tolist() == [True, False]
        a = summary.set_index("PDS_acc").loc["A"]
        assert (a.cluster_size, a.hits_in_results, a.total_representatives) == (3, 2, 2)
        assert (a.env_cli_ratio, a.unknown_source_count) == ("1:1", 1)
        assert (a.sampling_year_start, a.sampling_year_end, a.dated_isolates) == (
            2019,
            2021,
            2,
        )
        assert a.computed_serotype_known == 2
        hits.loc[0, "similarity_score"] = 0.998
        assert generate_cluster_table(conn, hits, 1000).near_top.all()
    finally:
        conn.close()


def test_annotation_rejects_unsafe_and_reserved_columns():
    for name in ("bad;DROP TABLE", "bad name", "1bad", "", "ASM_ACC", "PDS_acc"):
        with pytest.raises(ValueError):
            validate_column_name(name)


def run_cli(*args):
    result = subprocess.run(
        [sys.executable, "-m", "mashpit.mashpit", *map(str, args)],
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode == 0, result.stderr
    return result


def test_local_build_query_annotate_reshard_workflow(tmp_path):
    """One real offline CLI workflow replaces repeated builds in separate classes."""
    from random import Random

    inputs = tmp_path / "inputs"
    inputs.mkdir()
    for i in range(3):
        rng = Random(i)
        (inputs / f"sample{i}.fasta").write_text(
            ">contig\n" + "".join(rng.choice("ACGT") for _ in range(1000)) + "\n"
        )
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        "sample_id\tasm_acc\tstrain\tcustom_batch\nsample0\tGCA_override\tStrainA\tBatch1\n"
    )
    run_cli("build", "custom", "db", "--input-dir", inputs, "--metadata", metadata)
    db = tmp_path / "db" / "db.db"
    sig = tmp_path / "db" / "db.sig"
    before_hashes = {s.name: dict(s.minhash.hashes) for s in load_all_signatures(sig)}
    assert set(before_hashes) == {"GCA_override", "sample1", "sample2"}
    with sqlite3.connect(db) as conn:
        assert (
            conn.execute("SELECT value FROM DESC WHERE name='Type'").fetchone()[0]
            == "Custom"
        )
        assert conn.execute(
            "SELECT strain, custom_batch FROM METADATA WHERE asm_acc='sample1'"
        ).fetchone() == ("missing", "missing")
    values = tmp_path / "values.tsv"
    values.write_text("asm_acc\tproject_batch\nGCA_override\tUpdated\n")
    run_cli("annotate", "db", "--values", values)
    with sqlite3.connect(db) as conn:
        assert conn.execute(
            "SELECT strain, custom_batch, project_batch FROM METADATA WHERE asm_acc='GCA_override'"
        ).fetchone() == ("StrainA", "Batch1", "Updated")
        assert (
            conn.execute(
                "SELECT project_batch FROM METADATA WHERE asm_acc='sample1'"
            ).fetchone()[0]
            == "missing"
        )
    log = read_annotation_log(db)
    assert (
        len(log) == 1 and log[0][5] == hashlib.sha256(values.read_bytes()).hexdigest()
    )

    def query():
        run_cli("query", inputs / "sample0.fasta", "db", "--threshold", "0")
        frame = pd.read_csv("sample0_representative_matches.csv", index_col=0)
        assert frame.iloc[0]["asm_acc"] == "GCA_override"
        assert frame.iloc[0]["similarity_score"] == 1
        assert frame.iloc[0]["project_batch"] == "Updated"
        return frame.sort_values("asm_acc").reset_index(drop=True)

    expected = query()
    for shards in (2, 1):
        run_cli("reshard", "db", "--shards", shards)
        pd.testing.assert_frame_equal(query(), expected)
        assert {
            s.name: dict(s.minhash.hashes) for s in load_all_signatures(sig)
        } == before_hashes


def test_custom_build_rejects_duplicate_identity_and_cleans_up(tmp_path):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    for sample in ("a", "b"):
        (inputs / f"{sample}.fa").write_text(">c\n" + "ACGT" * 100 + "\n")
    metadata = tmp_path / "metadata.tsv"
    for key in ("asm_acc", "biosample_acc"):
        metadata.write_text(f"sample_id\t{key}\na\tshared\nb\tshared\n")
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "mashpit.mashpit",
                "build",
                "custom",
                "db",
                "--input-dir",
                str(inputs),
                "--metadata",
                str(metadata),
            ],
            capture_output=True,
            text=True,
            timeout=60,
        )
        assert result.returncode != 0 and "Duplicate identity" in result.stderr
        assert not (tmp_path / "db").exists()
