import io
import itertools
import os
import shutil
import unittest
import tempfile

import pyhmmer
from pyhmmer.errors import EaselError
from pyhmmer.easel import SequenceFile
from pyhmmer.plan7 import HMMFile, HMMPressedFile, Pipeline

from .. import __name__ as __package__
from .utils import HMMER_FOLDER, resource_files

# --- Mixins -------------------------------------------------------------------

class _TestHMMFile:

    ID = NotImplemented
    NAMES = NotImplemented

    @classmethod
    def setUpClass(cls):
        cls.hmms_folder = resource_files(__package__).joinpath("data", "hmms")

    def open_hmm(self, path):
        raise NotImplementedError()

    def check_hmmfile(self, hmmfile):
        for hmm, expected in itertools.zip_longest(hmmfile, self.NAMES):
            self.assertIsNot(hmm, None, "missing HMM: {}".format(expected))
            self.assertIsNot(expected, None, "unexpected extra HMM: {}".format(hmm))
            self.assertEqual(hmm.name, expected)
            self.assertIsNot(hmm.cutoffs, None)
            self.assertIsNot(hmm.evalue_parameters, None)

    def test_empty(self):
        with tempfile.NamedTemporaryFile() as empty:
            self.assertRaises(EOFError, self.open_hmm, empty.name)

    def test_read_hmmpressed(self):
        path = os.path.join(self.hmms_folder, "db", "{}.hmm".format(self.ID))
        with self.open_hmm(path) as f:
            self.check_hmmfile(f)

    def test_read_h3m(self):
        path = os.path.join(self.hmms_folder, "bin", "{}.h3m".format(self.ID))
        with self.open_hmm(path) as f:
            self.check_hmmfile(f)

    def test_read_hmm3(self):
        path = os.path.join(self.hmms_folder, "txt", "{}.hmm".format(self.ID))
        with self.open_hmm(path) as f:
            self.check_hmmfile(f)

    def test_read_hmm2(self):
        path = os.path.join(self.hmms_folder, "txt2", "{}.hmm2".format(self.ID))
        with self.open_hmm(path) as f:
            self.check_hmmfile(f)


class _TestHMMFileobj:

    def open_hmm(self, path):
        with open(path, "rb") as f:
            buffer = io.BytesIO(f.read())
        return HMMFile(buffer)

    def test_name(self):
        path = os.path.join(self.hmms_folder, "db", "{}.hmm".format(self.ID))
        with self.open_hmm(path) as f:
            self.assertIs(f.name, None)


class _TestHMMPath:

    def open_hmm(self, path):
        return HMMFile(path)

    def open_pressed(self, path):
        return HMMPressedFile(path)

    def test_init_error_filenotfound(self):
        self.assertRaises(FileNotFoundError, HMMFile, "path/to/missing/file")

    def test_init_error_folder(self):
        folder = tempfile.gettempdir()
        self.assertRaises(IsADirectoryError, HMMFile, folder)

    def test_read_optimized_profiles(self):
        path = os.path.join(self.hmms_folder, "db", "{}.hmm".format(self.ID))
        with self.open_hmm(path) as f:
            self.check_hmmfile(f.optimized_profiles())

    def test_rewind(self):
        path = os.path.join(self.hmms_folder, "txt", "{}.hmm".format(self.ID))
        with self.open_hmm(path) as f:
            hmm1 = f.read()
            f.rewind()
            hmm2 = f.read()
            self.assertIsNot(hmm1, hmm2)
            self.assertEqual(hmm1.name, hmm2.name)

    def test_rewind_optimized_profiles(self):
        path = os.path.join(self.hmms_folder, "db", "{}.hmm".format(self.ID))
        with self.open_pressed(path) as f:
            om1 = f.read()
            f.rewind()
            om2 = f.read()
            self.assertIsNot(om1, om2)
            self.assertEqual(om1.name, om2.name)

    def test_name(self):
        path = os.path.join(self.hmms_folder, "db", "{}.hmm".format(self.ID))
        with self.open_hmm(path) as f:
            self.assertEqual(f.name, path)
            self.assertEqual(f.optimized_profiles().name, path)


class _TestThioesterase(_TestHMMFile):
    ID = "Thioesterase"
    NAMES = [b"Thioesterase"]


class _TestT2PKS(_TestHMMFile):
    ID = "t2pks"
    NAMES = [
        b"CLF", b"CLF_7", b"CLF_8|9", b"CLF_11|12", b"AT", b"CYC", b"CYC_C7-C12",
        b"CYC_C5-C14", b"CYC_C5-C14/C3-C16", b"CYC_C1-C18|C2-C19", b"CYC_C2-C19",
        b"CYC_C5-C18", b"CYC_C4-C17/C2-C19", b"CYC_C4-C21/C2-C23|C2-C19",
        b"CYC_C9-C14", b"KSIII", b"ACP", b"KR", b"KR_C9", b"KR_C11", b"KR_C15",
        b"KR_C17", b"KR_C19", b"KS", b"OXY", b"GT", b"MET", b"MET_carboxy_O",
        b"MET_C2O|C2N", b"MET_C6|C8", b"MET_C9O", b"MET_C10", b"MET_C11O",
        b"MET_C13O|C17O", b"MET_C18O", b"DIMER", b"LIG", b"HAL", b"AMIN",
        b"AMID"
    ]


# --- Test cases ---------------------------------------------------------------

@unittest.skipUnless(resource_files, "importlib.resources.files not available")
class TestFileobjSingle(_TestHMMFileobj, _TestThioesterase, unittest.TestCase):
    pass

@unittest.skipUnless(resource_files, "importlib.resources.files not available")
class TestFileObjMultiple(_TestHMMFileobj, _TestT2PKS, unittest.TestCase):
    pass

@unittest.skipUnless(resource_files, "importlib.resources.files not available")
class TestPathSingle(_TestHMMPath, _TestThioesterase, unittest.TestCase):
    pass

@unittest.skipUnless(resource_files, "importlib.resources.files not available")
class TestPathMultiple(_TestHMMPath, _TestT2PKS, unittest.TestCase):
    pass
