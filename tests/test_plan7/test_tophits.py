import io
import itertools
import os
import shutil
import unittest
import tempfile
import pkg_resources

import pyhmmer
from pyhmmer.errors import EaselError
from pyhmmer.easel import SequenceFile
from pyhmmer.plan7 import HMMFile, Pipeline


class TestTopHits(unittest.TestCase):

    def setUp(self):
        hmm_path = pkg_resources.resource_filename("tests", "data/hmms/txt/PF02826.hmm")
        seqs_path = pkg_resources.resource_filename("tests", "data/seqs/938293.PRJEB85.HG003687.faa")

        with HMMFile(hmm_path) as hmm_file:
            self.hmm = next(hmm_file)
        with SequenceFile(seqs_path) as seqs_file:
            self.seqs = [seq.digitize(self.hmm.alphabet) for seq in seqs_file]

        self.pipeline = Pipeline(alphabet=self.hmm.alphabet)
        self.hits = self.pipeline.search(self.hmm, self.seqs)

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

    def test_sort_while_hit(self):
        # check sorting while a Hit instance exists does not crash, and that
        # the Hit still points to the right data
        hit = self.hits[-1]
        self.assertEqual(hit.name, b"938293.PRJEB85.HG003687_187")
        self.hits.sort(by="seqidx")
        if self.hits[-1].name != b"938293.PRJEB85.HG003687_187":
            self.assertEqual(hit.name, b"938293.PRJEB85.HG003687_187")

    def test_threshold(self):
        # after thresholding, one of the htis will not be included.
        self.assertEqual(self.hits.included, 15)
        self.assertEqual(self.hits.reported, 22)
