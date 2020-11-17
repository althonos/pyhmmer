import os
import unittest
import tempfile
import pkg_resources

import pyhmmer
from pyhmmer.easel import SequenceFile
from pyhmmer.plan7 import HMMFile, Pipeline


class TestHMMFile(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.hmms_folder = os.path.join(os.path.dirname(__file__), "data", "hmm")

    def test_init_error_empty(self):
        with tempfile.NamedTemporaryFile() as empty:
            # the file format can't be determined on an empty file, so this
            # should raise a value error
            self.assertRaises(ValueError, HMMFile, empty.name)

    def test_init_error_filenotfound(self):
        self.assertRaises(FileNotFoundError, HMMFile, "path/to/missing/file.hmm")

    def test_iter(self):
        hmm = os.path.join(self.hmms_folder, "Thioesterase.hmm")
        with HMMFile(hmm) as f:
            thioesterase = next(f)
            self.assertRaises(StopIteration, next, f)
        self.assertEqual(thioesterase.name, b"Thioesterase")


class TestTopHits(unittest.TestCase):

    def setUp(self):
        hmm_path = pkg_resources.resource_filename(__name__, "data/hmm/PF02826.hmm")
        seqs_path = pkg_resources.resource_filename(__name__, "data/seqs/938293.PRJEB85.HG003687.faa")

        with HMMFile(hmm_path) as hmm_file:
            self.hmm = next(hmm_file)
        with SequenceFile(seqs_path) as seqs_file:
            self.seqs = list(seqs_file)

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
        # check the hits are not sorted by default
        self.assertFalse(self.hits.is_sorted())

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
        hit = self.hits[0]
        self.assertEqual(hit.name, b"938293.PRJEB85.HG003688_14")
        self.hits.sort()
        self.assertNotEqual(self.hits[0].name, b"938293.PRJEB85.HG003688_14")
        self.assertEqual(hit.name, b"938293.PRJEB85.HG003688_14")

    def test_threshold(self):
        # before thresholding, all hits are reported / included
        self.assertEqual(self.hits.included, len(self.hits))
        self.assertEqual(self.hits.reported, len(self.hits))

        # thresholding should not segfault, even without sorting
        self.hits.threshold()

        # after thresholding, one of the htis will not be included.
        self.assertEqual(self.hits.included, 15)
        self.assertEqual(self.hits.reported, 16)
