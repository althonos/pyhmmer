import io
import itertools
import os
import platform
import unittest
import tempfile

import pyhmmer
from pyhmmer.errors import EaselError, AlphabetMismatch
from pyhmmer.easel import Alphabet, SequenceFile
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
        try:
            fd, filename = tempfile.mkstemp(suffix=".msa")
            self.assertTrue(os.path.exists(filename))
            self.assertRaises(EOFError, self.open_hmm, filename)
        finally:
            os.close(fd)
            if os.path.exists(filename):
                os.remove(filename)

    @unittest.skipIf(platform.system() == "Windows", "deadlocks on Windows")
    def test_read_hmmpressed(self):
        path = self.hmms_folder.joinpath("db", "{}.hmm".format(self.ID))
        if not path.exists():
            self.skipTest("data files not available")
        with self.open_hmm(path) as f:
            self.check_hmmfile(f)

    def test_read_h3m(self):
        path = self.hmms_folder.joinpath("bin", "{}.h3m".format(self.ID))
        if not path.exists():
            self.skipTest("data files not available")
        with self.open_hmm(path) as f:
            self.check_hmmfile(f)

    def test_read_hmm3(self):
        path = self.hmms_folder.joinpath("txt", "{}.hmm".format(self.ID))
        if not path.exists():
            self.skipTest("data files not available")
        with self.open_hmm(path) as f:
            self.check_hmmfile(f)

    def test_read_hmm2(self):
        path = self.hmms_folder.joinpath("txt2", "{}.hmm2".format(self.ID))
        if not path.exists():
            self.skipTest("data files not available")
        with self.open_hmm(path) as f:
            self.check_hmmfile(f)


class _TestHMMFileobj:

    def open_hmm(self, path, alphabet=None):
        with open(path, "rb") as f:
            buffer = io.BytesIO(f.read())
        return HMMFile(buffer, alphabet=alphabet)

    def test_name(self):
        path = self.hmms_folder.joinpath("bin", "{}.h3m".format(self.ID))
        if not path.exists():
            self.skipTest("data files not available")
        with self.open_hmm(path) as f:
            self.assertIs(f.name, None)

    def test_implicit_alphabet(self):
        path = self.hmms_folder.joinpath("bin", "{}.h3m".format(self.ID))
        if not path.exists():
            self.skipTest("data files not available")
        with self.open_hmm(path) as f:
            hmm = next(f)
            self.assertEqual(hmm.alphabet, self.ALPHABET)

    def test_explicit_alphabet(self):
        path = self.hmms_folder.joinpath("bin", "{}.h3m".format(self.ID))
        if not path.exists():
            self.skipTest("data files not available")
        with self.open_hmm(path, alphabet=self.ALPHABET) as f:
            hmm = next(f)
            self.assertEqual(hmm.alphabet, self.ALPHABET)

    def test_explicit_alphabet_mismatch(self):
        mm = Alphabet.amino() if self.ALPHABET.is_nucleotide() else Alphabet.dna()
        path = self.hmms_folder.joinpath("bin", "{}.h3m".format(self.ID))
        if not path.exists():
            self.skipTest("data files not available")
        with self.open_hmm(path, alphabet=mm) as f:
            self.assertRaises(AlphabetMismatch, f.read)


class _TestHMMPath:

    def open_hmm(self, path, alphabet=None):
        return HMMFile(path, alphabet=alphabet)

    def open_pressed(self, path):
        return HMMPressedFile(path)

    def test_init_error_filenotfound(self):
        self.assertRaises(FileNotFoundError, HMMFile, "path/to/missing/file")

    def test_init_error_folder(self):
        folder = tempfile.gettempdir()
        self.assertRaises(IsADirectoryError, HMMFile, folder)

    def test_read_optimized_profiles(self):
        path = self.hmms_folder.joinpath("db", "{}.hmm".format(self.ID))
        if not path.exists():
            self.skipTest("data files not available")
        with self.open_hmm(path) as f:
            self.check_hmmfile(f.optimized_profiles())

    def test_optimized_profiles_length(self):
        path = self.hmms_folder.joinpath("db", "{}.hmm".format(self.ID))
        if not path.exists():
            self.skipTest("data files not available")
        with self.open_hmm(path) as f:
            profiles = f.optimized_profiles()
            self.assertEqual(len(profiles), len(self.NAMES))
            profiles.read()
            self.assertEqual(len(profiles), len(self.NAMES) - 1)
            profiles.rewind()
            self.assertEqual(len(profiles), len(self.NAMES))

    def test_rewind_hmm(self):
        path = self.hmms_folder.joinpath("txt", "{}.hmm".format(self.ID))
        if not path.exists():
            self.skipTest("data files not available")
        with self.open_hmm(path) as f:
            hmm1 = f.read()
            f.rewind()
            hmm2 = f.read()
            self.assertIsNot(hmm1, hmm2)
            self.assertEqual(hmm1.name, hmm2.name)

    def test_rewind_h3m(self):
        path = self.hmms_folder.joinpath("bin", "{}.h3m".format(self.ID))
        if not path.exists():
            self.skipTest("data files not available")
        with self.open_hmm(path) as f:
            hmm1 = f.read()
            f.rewind()
            hmm2 = f.read()
            self.assertIsNot(hmm1, hmm2)
            self.assertEqual(hmm1.name, hmm2.name)

    def test_rewind_optimized_profiles(self):
        path = self.hmms_folder.joinpath("db", "{}.hmm".format(self.ID))
        if not path.exists():
            self.skipTest("data files not available")
        with self.open_pressed(path) as f:
            om1 = f.read()
            f.rewind()
            om2 = f.read()
            self.assertIsNot(om1, om2)
            self.assertEqual(om1.name, om2.name)

    def test_name_h3m(self):
        path = self.hmms_folder.joinpath("bin", "{}.h3m".format(self.ID))
        if not path.exists():
            self.skipTest("data files not available")
        with self.open_hmm(path) as f:
            self.assertEqual(f.name, str(path))

    def test_name_hmmpressed(self):
        path = self.hmms_folder.joinpath("db", "{}.hmm".format(self.ID))
        if not path.exists():
            self.skipTest("data files not available")
        with self.open_hmm(path) as f:
            self.assertEqual(f.name, str(path))
            self.assertEqual(f.optimized_profiles().name, str(path))


class _TestThioesterase(_TestHMMFile):
    ALPHABET = Alphabet.amino()
    ID = "Thioesterase"
    NAMES = ["Thioesterase"]


class _TestRREFam(_TestHMMFile):
    ALPHABET = Alphabet.amino()
    ID = "RREFam"
    NAMES = [
        "Stand_Alone_Lasso_RRE",
        "Thiopeptide_F_RRE",
        "PqqD_RRE",
        "Proteusin_Epimerase_RRE",
        "Thurincin_rSAM_RRE",
        "Thuricin_rSAM_RRE",
        "Other_Sactipeptide_rSAM_RRE",
        "Ranthipeptide_rSAM_RRE",
        "Trifolitoxin_RRE",
        "Thiaglutamate_B_RRE",        
    ]


# --- Test cases ---------------------------------------------------------------

@unittest.skipUnless(resource_files, "importlib.resources.files not available")
class TestFileobjSingle(_TestHMMFileobj, _TestThioesterase, unittest.TestCase):
    pass

@unittest.skipUnless(resource_files, "importlib.resources.files not available")
class TestFileObjMultiple(_TestHMMFileobj, _TestRREFam, unittest.TestCase):
    pass

@unittest.skipUnless(resource_files, "importlib.resources.files not available")
class TestPathSingle(_TestHMMPath, _TestThioesterase, unittest.TestCase):
    pass

@unittest.skipUnless(resource_files, "importlib.resources.files not available")
class TestPathMultiple(_TestHMMPath, _TestRREFam, unittest.TestCase):
    pass
