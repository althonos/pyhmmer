import copy
import gc
import io
import os
import unittest
import tempfile
from itertools import zip_longest

from pyhmmer import easel

from ..utils import EASEL_FOLDER


class TestMSAFile(unittest.TestCase):

    def test_guess_alphabet_empty_sequence(self):
        buffer = io.BytesIO(b">seq1\n\n")
        self.assertRaises(ValueError, easel.MSAFile, buffer, format="afa", digital=True)

    @unittest.skipUnless(os.path.exists(EASEL_FOLDER), "test data not available")
    def test_init_error_unknownformat(self):
        stockholm = os.path.join(EASEL_FOLDER, "formats", "stockholm.1")
        self.assertRaises(ValueError, easel.MSAFile, stockholm, format="nonsense")

    def test_init_error_empty(self):
        with tempfile.NamedTemporaryFile() as f:
            self.assertRaises(ValueError, easel.MSAFile, f.name)  # cannot guess format

    def test_init_file_not_found(self):
        self.assertRaises(FileNotFoundError, easel.MSAFile, "path/to/missing/file.sto")

    @unittest.skipUnless(os.path.exists(EASEL_FOLDER), "test data not available")
    def test_author(self):
        trna_5 = os.path.join(EASEL_FOLDER, "testsuite", "trna-5.stk")
        with easel.MSAFile(trna_5) as f:
            msa = f.read()
        self.assertEqual(msa.author, b"Infernal 0.1")

    @unittest.skipUnless(os.path.exists(EASEL_FOLDER), "test data not available")
    def test_read_wrong_format(self):
        stockholm = os.path.join(EASEL_FOLDER, "formats", "stockholm.1")
        # check reading while giving another file format fails
        with easel.MSAFile(stockholm, "clustal") as f:
            self.assertRaises(ValueError, f.read)

    @unittest.skipUnless(os.path.exists(EASEL_FOLDER), "test data not available")
    def test_iter(self):
        trna_5 = os.path.join(EASEL_FOLDER, "testsuite", "trna-5.stk")
        with easel.MSAFile(trna_5) as f:
            msa = next(f)
            self.assertRaises(StopIteration, next, f)

    @unittest.skipUnless(os.path.exists(EASEL_FOLDER), "test data not available")
    def test_closed_file(self):
        trna_5 = os.path.join(EASEL_FOLDER, "testsuite", "trna-5.stk")
        with easel.MSAFile(trna_5) as f:
            pass
        self.assertRaises(ValueError, f.read)


class _TestReadFilename(object):

    def test_read_filename_guess_format(self):
        for filename, start, count, length in zip_longest(self.filenames, self.starts, self.counts, self.lengths):
            path = os.path.join(self.folder, filename)
            with easel.MSAFile(path) as f:
                msa = f.read()
                self.assertEqual(msa.sequences[0].sequence[:10], start)
                self.assertEqual(len(msa), length)
                self.assertEqual(len(msa.sequences), count)

    def test_read_filename_given_format(self):
        for filename, start, count, length in zip_longest(self.filenames, self.starts, self.counts, self.lengths):
            path = os.path.join(self.folder, filename)
            with easel.MSAFile(path, self.format) as f:
                msa = f.read()
                self.assertEqual(msa.sequences[0].sequence[:10], start)
                self.assertEqual(len(msa), length)
                self.assertEqual(len(msa.sequences), count)

    def test_guess_alphabet(self):
        for filename, alphabet in zip_longest(self.filenames, self.alphabet):
            path = os.path.join(self.folder, filename)
            with easel.MSAFile(path, self.format, digital=True) as f:
                self.assertEqual(f.alphabet, alphabet)


class _TestReadFileObject(object):

    def test_read_fileobject_guess_format(self):
        for filename, start, count, length in zip_longest(self.filenames, self.starts, self.counts, self.lengths):
            path = os.path.join(self.folder, filename)
            with open(path, "rb") as f:
                buffer = io.BytesIO(f.read())
            with easel.MSAFile(buffer) as f:
                msa = f.read()
                self.assertEqual(msa.sequences[0].sequence[:10], start)
                self.assertEqual(len(msa), length)
                self.assertEqual(len(msa.sequences), count)

    def test_read_filename_given_format(self):
        for filename, start, count, length in zip_longest(self.filenames, self.starts, self.counts, self.lengths):
            path = os.path.join(self.folder, filename)
            with open(path, "rb") as f:
                buffer = io.BytesIO(f.read())
            with easel.MSAFile(buffer, self.format) as f:
                msa = f.read()
                self.assertEqual(msa.sequences[0].sequence[:10], start)
                self.assertEqual(len(msa), length)
                self.assertEqual(len(msa.sequences), count)

    def test_guess_alphabet(self):
        for filename, alphabet in zip_longest(self.filenames, self.alphabet):
            path = os.path.join(self.folder, filename)
            with open(path, "rb") as f:
                buffer = io.BytesIO(f.read())
            with easel.MSAFile(buffer, self.format) as f:
                self.assertEqual(f.alphabet, alphabet)


@unittest.skipUnless(os.path.exists(EASEL_FOLDER), "test data not available")
class TestA2MFile(_TestReadFilename, _TestReadFileObject, unittest.TestCase):
    folder    = os.path.join(EASEL_FOLDER, "esl_msa_testfiles", "a2m")
    filenames = ["a2m.good.1", "a2m.good.2"]
    starts    = ["VLSEGEWQLV", "UCCGAUAUAG"]
    format    = "a2m"
    lengths   = [165, 73]
    counts    = [4, 5]
    alphabet  = [easel.Alphabet.amino(), easel.Alphabet.rna()]

    @unittest.expectedFailure
    def test_read_filename_guess_format(self):
        super().test_read_filename_guess_format()

    @unittest.expectedFailure
    def test_read_fileobject_guess_format(self):
        super().test_read_fileobject_guess_format()


@unittest.skipUnless(os.path.exists(EASEL_FOLDER), "test data not available")
class TestAlignedFastaFile(_TestReadFilename, _TestReadFileObject, unittest.TestCase):
    folder    = os.path.join(EASEL_FOLDER, "esl_msa_testfiles", "afa")
    filenames = ["afa.good.1", "afa.good.2", "afa.good.3"]
    starts    = ["VLSEGEWQLV", "UCCGAUAUAG", "mqifvktltg"]
    format    = "afa"
    lengths   = [165, 73, 128]
    counts    = [4, 5, 6]
    alphabet  = [easel.Alphabet.amino(), easel.Alphabet.rna(), easel.Alphabet.amino()]


@unittest.skipUnless(os.path.exists(EASEL_FOLDER), "test data not available")
class TestClustalFile(_TestReadFilename, _TestReadFileObject, unittest.TestCase):
    folder    = os.path.join(EASEL_FOLDER, "esl_msa_testfiles", "clustal")
    filenames = ["clustal.good.2"]
    starts    = ["UCCGAUAUAG"]
    format    = "clustallike"
    lengths   = [73]
    counts    = [5]
    alphabet  = [easel.Alphabet.rna()]


@unittest.skipUnless(os.path.exists(EASEL_FOLDER), "test data not available")
class TestClustalLikeFile(_TestReadFilename, _TestReadFileObject, unittest.TestCase):
    folder    = os.path.join(EASEL_FOLDER, "esl_msa_testfiles", "clustal")
    filenames = ["clustal.good.1"]
    starts    = ["VLSEGEWQLV"]
    format    = "clustallike"
    lengths   = [165]
    counts    = [4]
    alphabet  = [easel.Alphabet.amino()]


@unittest.skipUnless(os.path.exists(EASEL_FOLDER), "test data not available")
class TestPhylipFile(_TestReadFilename, _TestReadFileObject, unittest.TestCase):
    folder    = os.path.join(EASEL_FOLDER, "esl_msa_testfiles", "phylip")
    filenames = ["phylip.good.1", "phylip.good.2", "phylip.good.3"]
    starts    = ["AAGCTNGGGC", "ATGGCGAAGG", "MKVILLFVLA"]
    format    = "phylip"
    lengths   = [42, 50, 384]
    counts    = [5, 7, 3]
    alphabet  = [easel.Alphabet.dna(), easel.Alphabet.dna(), easel.Alphabet.amino()]


@unittest.skipUnless(os.path.exists(EASEL_FOLDER), "test data not available")
class TestPhylipsFile(_TestReadFilename, _TestReadFileObject, unittest.TestCase):
    folder    = os.path.join(EASEL_FOLDER, "esl_msa_testfiles", "phylips")
    filenames = ["phylips.good.1", "phylips.good.2"]
    starts    = ["AAGCTNGGGC", "MKVILLFVLA"]
    format    = "phylips"
    lengths   = [42, 384]
    counts    = [5, 3]
    alphabet  = [easel.Alphabet.dna(), easel.Alphabet.amino()]


@unittest.skipUnless(os.path.exists(EASEL_FOLDER), "test data not available")
class TestPsiblastFile(_TestReadFilename, _TestReadFileObject, unittest.TestCase):
    folder    = os.path.join(EASEL_FOLDER, "esl_msa_testfiles", "psiblast")
    filenames = ["psiblast.good.1", "psiblast.good.2"]
    starts    = ["VLSEGEWQLV", "UCCGAUAUAG"]
    format    = "psiblast"
    lengths   = [165, 73]
    counts    = [4, 5]
    alphabet  = [easel.Alphabet.amino(), easel.Alphabet.rna()]

    @unittest.expectedFailure
    def test_read_fileobject_guess_format(self):
        super().test_read_fileobject_guess_format()


@unittest.skipUnless(os.path.exists(EASEL_FOLDER), "test data not available")
class TestSelexFile(_TestReadFilename, _TestReadFileObject, unittest.TestCase):
    folder    = os.path.join(EASEL_FOLDER, "esl_msa_testfiles", "selex")
    filenames = ["selex.good.1", "selex.good.2", "selex.good.3", "selex.good.4"]
    starts    = ["ACDEFGHIKL", "ACDEFGHIKL", "ACDEFGHIKL", "gGAGUAAGAU"]
    format    = "selex"
    lengths   = [20, 40, 52, 31]
    counts    = [5, 5, 7, 11]
    alphabet  = [easel.Alphabet.amino(), easel.Alphabet.amino(), easel.Alphabet.amino(), easel.Alphabet.rna()]

    @unittest.expectedFailure
    def test_read_fileobject_guess_format(self):
        super().test_read_fileobject_guess_format()


@unittest.skipUnless(os.path.exists(EASEL_FOLDER), "test data not available")
class TestStockholmFile(_TestReadFilename, _TestReadFileObject, unittest.TestCase):
    folder    = os.path.join(EASEL_FOLDER, "esl_msa_testfiles", "stockholm")
    filenames = ["stockholm.good.1"]
    starts    = ["ACDEFGHKLM"]
    format    = "stockholm"
    lengths   = [38]
    counts    = [7]
    alphabet  = [easel.Alphabet.amino()]
