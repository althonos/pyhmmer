import io
import itertools
import os
import pickle
import shutil
import unittest
import tempfile
import pkg_resources

import pyhmmer
from pyhmmer.errors import EaselError
from pyhmmer.easel import SequenceFile, TextMSA, DigitalMSA
from pyhmmer.plan7 import HMMFile, Pipeline, TopHits


class TestTopHits(unittest.TestCase):

    def assertAlignmentEqual(self, a1, a2):
        for attr in (
            "hmm_from",
            "hmm_to",
            "hmm_name",
            "hmm_accession",
            "hmm_sequence",
            "target_from",
            "target_to",
            "target_name",
            "target_sequence",
            "identity_sequence",
        ):
            self.assertEqual(
                getattr(a1, attr),
                getattr(a2, attr),
                "attribute {!r} differs".format(attr)
            )

    def assertDomainEqual(self, d1, d2):
        for attr in (
            "env_from",
            "env_to",
            "score",
            "bias",
            "correction",
            "envelope_score",
            "c_evalue",
            "i_evalue",
            "pvalue"
        ):
            self.assertEqual(
                getattr(d1, attr),
                getattr(d2, attr),
                "attribute {!r} differs".format(attr)
            )
        self.assertAlignmentEqual(d1.alignment, d2.alignment)

    def assertHitEqual(self, h1, h2):
        for attr in (
            "name",
            "accession",
            "description",
            "score",
            "pre_score",
            "sum_score",
            "bias",
            "evalue",
            "pvalue",
        ):
            self.assertEqual(
                getattr(h1, attr),
                getattr(h2, attr),
                "attribute {!r} differs".format(attr)
            )
        self.assertEqual(len(h1.domains), len(h2.domains))
        for d1, d2 in zip(h1.domains, h2.domains):
            self.assertDomainEqual(d1, d2)

    def assertHitsEqual(self, hits1, hits2):
        self.assertEqual(hits1.query_name, hits2.query_name)
        self.assertEqual(hits1.query_accession, hits2.query_accession)
        self.assertEqual(len(hits1), len(hits2))
        for h1, h2 in zip(hits1, hits2):
            self.assertHitEqual(h1, h2)

    def setUp(self):
        hmm_path = pkg_resources.resource_filename("pyhmmer.tests", "data/hmms/txt/PF02826.hmm")
        seqs_path = pkg_resources.resource_filename("pyhmmer.tests", "data/seqs/938293.PRJEB85.HG003687.faa")

        with HMMFile(hmm_path) as hmm_file:
            self.hmm = hmm_file.read()
        with SequenceFile(seqs_path, digital=True, alphabet=self.hmm.alphabet) as seqs_file:
            self.seqs = seqs_file.read_block()

        self.pipeline = Pipeline(alphabet=self.hmm.alphabet)
        self.hits = self.pipeline.search_hmm(self.hmm, self.seqs)
        self.pipeline.clear()

    def test_bool(self):
        self.assertFalse(pyhmmer.plan7.TopHits())
        self.assertTrue(self.hits)

    def test_index_error(self):
        # check out of bounds indexing is caught
        with self.assertRaises(IndexError):
            self.hits[len(self.hits) + 1]
        # check out of bounds negative indexing works
        with self.assertRaises(IndexError):
            self.hits[-len(self.hits) - 1]

    def test_index_negative(self):
        dom = self.hits[len(self.hits) - 1]
        dom_last = self.hits[-1]
        self.assertEqual(dom.name, dom_last.name)

    def test_Z(self):
        empty = TopHits()
        self.assertEqual(empty.Z, 0)
        self.assertEqual(self.hits.Z, len(self.seqs))

    def test_strand(self):
        empty = TopHits()
        self.assertIs(empty.strand, None)

    def test_searched_sequences(self):
        empty = TopHits()
        self.assertEqual(empty.searched_sequences, 0)
        self.assertEqual(self.hits.searched_sequences, len(self.seqs))

    def test_searched_nodes(self):
        empty = TopHits()
        self.assertEqual(empty.searched_nodes, 0)
        self.assertEqual(self.hits.searched_nodes, self.hmm.M)

    def test_searched_residues(self):
        empty = TopHits()
        self.assertEqual(empty.searched_residues, 0)
        self.assertEqual(self.hits.searched_residues, sum(map(len, self.seqs)))

    def test_searched_models(self):
        empty = TopHits()
        self.assertEqual(empty.searched_sequences, 0)
        self.assertEqual(self.hits.searched_models, 1)

    def test_merge_empty(self):
        empty = TopHits()
        self.assertFalse(empty.long_targets)
        self.assertEqual(empty.Z, 0.0)
        self.assertEqual(empty.domZ, 0.0)

        empty2 = empty.merge()
        self.assertFalse(empty2.long_targets)
        self.assertEqual(empty2.Z, 0.0)
        self.assertEqual(empty2.domZ, 0.0)

        merged_empty = empty.merge(TopHits())
        self.assertHitsEqual(merged_empty, empty)
        self.assertEqual(merged_empty.searched_residues, 0)
        self.assertEqual(merged_empty.searched_sequences, 0)
        self.assertFalse(merged_empty.long_targets)
        self.assertFalse(self.hits.long_targets)
        self.assertEqual(merged_empty.Z, 0.0)
        self.assertEqual(merged_empty.domZ, 0.0)

        merged = empty.merge(self.hits)
        self.assertHitsEqual(merged, self.hits)

    def test_merge_pipeline(self):
        hits1 = self.pipeline.search_hmm(self.hmm, self.seqs[:1000])
        hits2 = self.pipeline.search_hmm(self.hmm, self.seqs[1000:2000])
        hits3 = self.pipeline.search_hmm(self.hmm, self.seqs[2000:])

        self.assertEqual(hits1.Z, 1000)
        self.assertEqual(hits2.Z, 1000)
        self.assertEqual(hits3.Z, len(self.seqs) - 2000)
        self.assertEqual(len(hits1) + len(hits2) + len(hits3), len(self.hits))

        merged = hits1.merge(hits2, hits3)
        self.assertEqual(merged.searched_sequences, self.hits.searched_sequences)
        self.assertEqual(merged.searched_models, self.hits.searched_models)
        self.assertEqual(merged.Z, self.hits.Z)
        self.assertEqual(merged.domZ, self.hits.domZ)
        self.assertHitsEqual(merged, self.hits)

    def test_merged_pipeline_fixed_Z(self):
        pipeline = Pipeline(alphabet=self.hmm.alphabet, Z=200.0)
        hits1 = pipeline.search_hmm(self.hmm, self.seqs[:1000])
        hits2 = pipeline.search_hmm(self.hmm, self.seqs[1000:2000])
        hits3 = pipeline.search_hmm(self.hmm, self.seqs[2000:])

        self.assertEqual(hits1.Z, 200.0)
        self.assertEqual(hits2.Z, 200.0)
        self.assertEqual(hits3.Z, 200.0)

        merged = hits1.merge(hits2, hits3)
        self.assertEqual(merged.Z, 200.0)

        hits = pipeline.search_hmm(self.hmm, self.seqs)
        self.assertHitsEqual(merged, hits)

    def test_copy(self):
        copy = self.hits.copy()
        self.assertHitsEqual(copy, self.hits)

    def test_sort(self):
        # check the hits are sorted by default
        self.assertTrue(self.hits.is_sorted())

        # check sorting
        self.hits.sort()
        self.assertTrue(self.hits.is_sorted())
        self.assertFalse(self.hits.is_sorted(by="seqidx"))

        # check sorting using a different key
        self.hits.sort(by="seqidx")
        self.assertTrue(self.hits.is_sorted(by="seqidx"))
        self.assertFalse(self.hits.is_sorted(by="key"))

        # check error otherwise
        self.assertRaises(ValueError, self.hits.sort, by="nonsense")
        self.assertRaises(ValueError, self.hits.is_sorted, by="nonsense")

    def test_sort_while_hit(self):
        # check sorting while a Hit instance exists does not crash, and that
        # the Hit still points to the right data
        hit = self.hits[-1]
        self.assertEqual(hit.name, b"938293.PRJEB85.HG003687_187")
        self.hits.sort(by="seqidx")
        if self.hits[-1].name != b"938293.PRJEB85.HG003687_187":
            self.assertEqual(hit.name, b"938293.PRJEB85.HG003687_187")

    def test_threshold(self):
        # after thresholding, one of the hits will not be included.
        self.assertEqual(self.hits.hits_included, 15)
        self.assertEqual(self.hits.hits_reported, 22)

    def test_to_msa(self):
        msa = self.hits.to_msa(self.hmm.alphabet)
        self.assertIsInstance(msa, TextMSA)
        unique_names = { seq.name.split(b"/")[0] for seq in msa.sequences }
        self.assertEqual(len(unique_names), self.hits.hits_included)

        msa_d = self.hits.to_msa(
            self.hmm.alphabet,
            trim=True,
            digitize=True,
            all_consensus_cols=True
        )
        self.assertIsInstance(msa_d, DigitalMSA)

    def test_pickle(self):
        pickled = pickle.loads(pickle.dumps(self.hits))
        self.assertHitsEqual(pickled, self.hits)

    def test_query_name(self):
        self.assertEqual(self.hits.query_name, self.hmm.name)

    def test_query_accession(self):
        self.assertEqual(self.hits.query_accession, self.hmm.accession)

    def test_write_target(self):
        buffer = io.BytesIO()
        self.hits.write(buffer, format="targets")
        lines = buffer.getvalue().decode().splitlines()

        expected = pkg_resources.resource_string("pyhmmer.tests", "data/tables/PF02826.tbl").decode().splitlines()
        while expected[-1].startswith("#"):
            expected.pop()

        self.assertMultiLineEqual("\n".join(lines), "\n".join(expected))

    def test_write_domains(self):
        buffer = io.BytesIO()
        self.hits.write(buffer, format="domains")
        lines = buffer.getvalue().decode().splitlines()

        expected = pkg_resources.resource_string("pyhmmer.tests", "data/tables/PF02826.domtbl").decode().splitlines()
        while expected[-1].startswith("#"):
            expected.pop()

        self.assertMultiLineEqual("\n".join(lines), "\n".join(expected))
