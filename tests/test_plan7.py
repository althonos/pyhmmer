import os
import unittest
import tempfile

from pyhmmer import plan7


class TestHMMFile(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.hmms_folder = os.path.join(os.path.dirname(__file__), "data", "hmm")

    def test_init_error_empty(self):
        with tempfile.NamedTemporaryFile() as empty:
            # the file format can't be determined on an empty file, so this
            # should raise a value error
            self.assertRaises(ValueError, plan7.HMMFile, empty.name)

    def test_init_error_filenotfound(self):
        self.assertRaises(FileNotFoundError, plan7.HMMFile, "path/to/missing/file.hmm")

    def test_iter(self):
        hmm = os.path.join(self.hmms_folder, "Thioesterase.hmm")
        with plan7.HMMFile(hmm) as f:
            thioesterase = next(f)
            self.assertRaises(StopIteration, next, f)
        self.assertEqual(thioesterase.name, b"Thioesterase")
