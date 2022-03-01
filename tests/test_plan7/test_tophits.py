import io
import itertools
import os
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
        self.assertEqual(len(hits1), len(hits2))
        for h1, h2 in zip(hits1, hits2):
            self.assertHitEqual(h1, h2)

    def setUp(self):
        hmm_path = pkg_resources.resource_filename("tests", "data/hmms/txt/PF02826.hmm")
        seqs_path = pkg_resources.resource_filename("tests", "data/seqs/938293.PRJEB85.HG003687.faa")

        with HMMFile(hmm_path) as hmm_file:
            self.hmm = next(hmm_file)
        with SequenceFile(seqs_path) as seqs_file:
            seqs_file.set_digital(self.hmm.alphabet)
            self.seqs = list(seqs_file)

        self.pipeline = Pipeline(alphabet=self.hmm.alphabet)
        self.hits = self.pipeline.search_hmm(self.hmm, self.seqs)

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

    def test_merge_empty(self):
        hits = TopHits()
        merged = hits.merge(self.hits)
        self.assertHitsEqual(merged, self.hits)

    def test_merge_pipeline_results(self):
        hits1 = self.pipeline.search_hmm(self.hmm, self.seqs[:1000])
        hits2 = self.pipeline.search_hmm(self.hmm, self.seqs[1000:])

        self.assertEqual(hits1.Z, 1000.0)
        self.assertEqual(hits2.Z, len(self.seqs) - 1000)
        self.assertEqual(len(hits1) + len(hits2), len(self.hits))

        merged = hits1.merge(hits2)
        self.assertHitsEqual(merged, self.hits)

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
        self.assertEqual(self.hits.included, 15)
        self.assertEqual(self.hits.reported, 22)

    def test_to_msa(self):
        msa = self.hits.to_msa(self.hmm.alphabet)
        self.assertIsInstance(msa, TextMSA)
        unique_names = { seq.name.split(b"/")[0] for seq in msa.sequences }
        self.assertEqual(len(unique_names), self.hits.included)

        msa_d = self.hits.to_msa(self.hmm.alphabet, True, True, True)
        self.assertIsInstance(msa_d, DigitalMSA)
