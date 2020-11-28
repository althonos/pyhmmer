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


class TestHMM(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.hmms_folder = os.path.join(os.path.dirname(__file__), "data", "hmm")
        cls.hmm_path = os.path.join(cls.hmms_folder, "Thioesterase.hmm")
        with HMMFile(cls.hmm_path) as hmm_file:
            cls.hmm = next(hmm_file)

    def test_write(self):
        buffer = io.BytesIO()
        self.hmm.write(buffer)

        self.assertNotEqual(buffer.tell(), 0)
        buffer.seek(0)

        with open(self.hmm_path, "rb") as f:
            # 1st line is skipped, cause it contains the format date, which will
            # obviously not be the same as the reference
            for line_written, line_ref in itertools.islice(zip(buffer, f), 1, None):
                self.assertEqual(line_written, line_ref)

    def test_write_error(self):
        buffer = io.StringIO()

        with self.assertRaises(EaselError) as ctx:
            self.hmm.write(buffer)
        self.assertEqual(ctx.exception.code, 27)
        self.assertIsInstance(ctx.exception.__cause__, TypeError)


class _TestHMMFile:

    @classmethod
    def setUpClass(cls):
        cls.hmms_folder = os.path.join(os.path.dirname(__file__), "data", "hmm")

    def open_hmm(self, path):
        return NotImplemented

    def test_init_empty(self):
        with tempfile.NamedTemporaryFile() as empty:
            self.assertRaises(EOFError, self.open_hmm, empty.name)

    def test_read_h3m(self):
        path = os.path.join(self.hmms_folder, "Thioesterase.hmm.h3m")
        with self.open_hmm(path) as f:
            self.assertEqual(next(f).name, b"Thioesterase")

    def test_read_hmm3(self):
        path = os.path.join(self.hmms_folder, "Thioesterase.hmm")
        with self.open_hmm(path) as f:
            self.assertEqual(next(f).name, b"Thioesterase")

    def test_read_hmm2(self):
        path = os.path.join(self.hmms_folder, "PKSI-AT.hmm2")
        with self.open_hmm(path) as f:
            self.assertEqual(next(f).name, b"PKS-AT.tcoffee")

    def test_iter(self):
        path = os.path.join(self.hmms_folder, "Thioesterase.hmm")
        with self.open_hmm(path) as f:
            thioesterase = next(f)
            self.assertRaises(StopIteration, next, f)
        self.assertEqual(thioesterase.name, b"Thioesterase")


class TestHMMFilePath(_TestHMMFile, unittest.TestCase):

    def open_hmm(self, path):
        return HMMFile(path)

    def test_init_filenotfound(self):
        self.assertRaises(FileNotFoundError, HMMFile, "path/to/missing/file.hmm")


class TestHMMFileFileobj(_TestHMMFile, unittest.TestCase):

    def setUp(self):
        self.files = []

    def tearDown(self):
        for f in self.files:
            f.close()

    def open_hmm(self, path):

        with open(path, "rb") as f:
            buffer = io.BytesIO(f.read())

        # self.files.append(open(path, "rb"))
        return HMMFile(buffer)


class TestTopHits(unittest.TestCase):

    def setUp(self):
        hmm_path = pkg_resources.resource_filename(__name__, "data/hmm/PF02826.hmm")
        seqs_path = pkg_resources.resource_filename(__name__, "data/seqs/938293.PRJEB85.HG003687.faa")

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
        self.assertEqual(hit.name, b"938293.PRJEB85.HG003687_154")
        self.hits.sort(by="seqidx")
        if self.hits[-1].name != b"938293.PRJEB85.HG003687_154":
            self.assertEqual(hit.name, b"938293.PRJEB85.HG003687_154")

    def test_threshold(self):
        # after thresholding, one of the htis will not be included.
        self.assertEqual(self.hits.included, 15)
        self.assertEqual(self.hits.reported, 16)
