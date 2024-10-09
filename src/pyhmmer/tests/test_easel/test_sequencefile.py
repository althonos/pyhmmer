import copy
import gc
import io
import os
import unittest
import tempfile
from itertools import zip_longest

from pyhmmer import easel

from .. import __name__ as __package__
from .utils import EASEL_FOLDER, resource_files


class TestSequenceFile(unittest.TestCase):

    def test_guess_alphabet_empty_sequence(self):
        buffer = io.BytesIO(b">seq1\n\n")
        self.assertRaises(ValueError, easel.SequenceFile, buffer, format="fasta", digital=True)

    def test_init_error_unknownformat(self):
        with self.assertRaises(ValueError):
            _file = easel.SequenceFile("file.x", format="nonsense")

    def test_init_error_empty(self):
        with tempfile.NamedTemporaryFile() as empty:
            # without a format argument, we can't determine the file format
            self.assertRaises(ValueError, easel.SequenceFile, empty.name)
            # with a format argument, we don't expect the file to be empty
            self.assertRaises(EOFError, easel.SequenceFile, empty.name, "fasta")
        with io.BytesIO(b"") as buffer:
            # without a format argument, we can't determine the file format
            self.assertRaises(ValueError, easel.SequenceFile, buffer)
            # with a format argument, we don't expect the file to be empty
            self.assertRaises(EOFError, easel.SequenceFile, buffer, "fasta")

    def test_init_error_filenotfound(self):
        self.assertRaises(
            FileNotFoundError, easel.SequenceFile, "path/to/missing/file.fa"
        )

    def test_init_error_folder(self):
        folder = tempfile.gettempdir()
        self.assertRaises(IsADirectoryError, easel.SequenceFile, folder)

    @unittest.skipUnless(os.path.exists(EASEL_FOLDER), "test data not available")
    def test_iter(self):
        fasta = os.path.join(EASEL_FOLDER, "formats", "fasta")

        with easel.SequenceFile(fasta) as f:
            seq1, seq2 = list(f)

        self.assertEqual(seq1.name, b"SNRPA_DROME")
        self.assertEqual(seq2.name, b"SNRPA_HUMAN")

    def test_parse(self):
        seq = easel.SequenceFile.parse(b"> EcoRI\nGAATTC\n", format="fasta")
        self.assertEqual(seq.name, b"EcoRI")

    def test_parse_error_unknownformat(self):
        with self.assertRaises(ValueError):
            _file = easel.SequenceFile.parse(b"> EcoRI\nGAATTC\n", format="nonsense")

    @unittest.skipUnless(os.path.exists(EASEL_FOLDER), "test data not available")
    def test_read_skip_info(self):
        fasta = os.path.join(EASEL_FOLDER, "formats", "fasta")

        with easel.SequenceFile(fasta) as f:
            seq = f.read()
            self.assertEqual(seq.name, b"SNRPA_DROME")
            self.assertTrue(seq.sequence.startswith("MEMLPNQTIY"))

        with easel.SequenceFile(fasta) as f:
            seq = f.read(skip_info=True)
            self.assertEqual(seq.name, b"")
            self.assertTrue(seq.sequence.startswith("MEMLPNQTIY"))

    @unittest.skipUnless(os.path.exists(EASEL_FOLDER), "test data not available")
    def test_read_skip_sequence(self):
        fasta = os.path.join(EASEL_FOLDER, "formats", "fasta")

        with easel.SequenceFile(fasta) as f:
            seq = f.read()
            self.assertEqual(seq.name, b"SNRPA_DROME")
            self.assertTrue(seq.sequence.startswith("MEMLPNQTIY"))

        with easel.SequenceFile(fasta) as f:
            seq = f.read(skip_sequence=True)
            self.assertEqual(seq.name, b"SNRPA_DROME")
            self.assertEqual(seq.sequence, "")

    @unittest.skipUnless(resource_files, "importlib.resources.files not available")
    def test_ignore_gaps(self):
        luxc = resource_files(__package__).joinpath("data", "msa", "LuxC.faa")
        # fails if not ignoring gaps (since if contains gaps)
        with easel.SequenceFile(luxc, "fasta") as seq_file:
            self.assertRaises(ValueError, seq_file.read)
        # succeeds if ignoring gaps
        with easel.SequenceFile(luxc, "afa") as seq_file:
            sequences = list(seq_file)
            self.assertEqual(len(sequences), 13)

class _TestReadFilename(object):

    def test_read_filename_guess_format(self):
        # check reading a file without specifying the format works
        for filename, start in zip_longest(self.filenames, self.starts):
            path = os.path.join(self.folder, filename)
            with easel.SequenceFile(path) as f:
                seq = f.read()
                self.assertEqual(seq.sequence[:10], start)

    def test_read_filename_given_format(self):
        # check reading a file while specifying the format works
        for filename, start in zip_longest(self.filenames, self.starts):
            path = os.path.join(self.folder, filename)
            with easel.SequenceFile(path, self.format) as f:
                seq = f.read()
                self.assertEqual(seq.sequence[:10], start)

    def test_read_filename_count_sequences(self):
        # check reading a file while specifying the format works
        for filename, count in zip_longest(self.filenames, self.counts):
            path = os.path.join(self.folder, filename)
            with easel.SequenceFile(path, self.format) as f:
                sequences = list(f)
                self.assertEqual(len(sequences), count)

    def test_read_filename_guess_alphabet(self):
        for filename, alphabet in zip_longest(self.filenames, self.alphabet):
            path = os.path.join(self.folder, filename)
            if alphabet is not None:
                with easel.SequenceFile(path, self.format, digital=True) as f:
                    self.assertEqual(f.alphabet, alphabet)
                # FIXME: Until EddyRivasLab/easel#61 is merged, we cannot
                #        test the case where `guess_alphabet` would throw
                #        an error because of a bug in `sqascii_GuessAlphabet`
                #        causing `eslEOD` to be returned when `eslNOALPHABET`
                #        is expected.

    def test_read_filename_given_alphabet(self):
        for filename, alphabet in zip_longest(self.filenames, self.alphabet):
            path = os.path.join(self.folder, filename)
            if alphabet is not None:
                with easel.SequenceFile(path, self.format, digital=True, alphabet=alphabet) as f:
                    self.assertEqual(f.alphabet, alphabet)

    def test_rewind_filename(self):
        for filename, count, alphabet in zip_longest(self.filenames, self.counts, self.alphabet):
            if alphabet is None:
                continue
            path = os.path.join(self.folder, filename)
            with easel.SequenceFile(path, self.format, digital=True) as f:
                seqs1 = f.read_block()
                f.rewind()
                seqs2 = f.read_block()
                self.assertEqual(len(seqs1), count)
                self.assertEqual(len(seqs2), count)
                self.assertEqual(seqs1, seqs2)


class _TestReadFileObject(object):

    def test_read_fileobject_guess_format(self):
        # check reading a file while specifying the format works
        for filename, start in zip_longest(self.filenames, self.starts):
            path = os.path.join(self.folder, filename)
            with open(path, "rb") as f:
                buffer = io.BytesIO(f.read())
            with easel.SequenceFile(buffer, self.format) as f:
                seq = f.read()
                self.assertEqual(seq.sequence[:10], start)

    def test_read_fileobject_given_format(self):
        # check reading a file while specifying the format works
        for filename, start in zip_longest(self.filenames, self.starts):
            path = os.path.join(self.folder, filename)
            with open(path, "rb") as f:
                buffer = io.BytesIO(f.read())
            with easel.SequenceFile(buffer, self.format) as f:
                seq = f.read()
                self.assertEqual(seq.sequence[:10], start)

    def test_read_fileobject_count_sequences(self):
        # check reading a file while specifying the format works
        for filename, count in zip_longest(self.filenames, self.counts):
            path = os.path.join(self.folder, filename)
            with open(path, "rb") as f:
                buffer = io.BytesIO(f.read())
            with easel.SequenceFile(buffer, self.format) as f:
                sequences = list(f)
                self.assertEqual(len(sequences), count)

    def test_read_fileobject_guess_alphabet(self):
        for filename, alphabet in zip_longest(self.filenames, self.alphabet):
            path = os.path.join(self.folder, filename)
            with open(path, "rb") as f:
                buffer = io.BytesIO(f.read())
            if alphabet is not None:
                with easel.SequenceFile(buffer, self.format, digital=True) as f:
                    self.assertEqual(f.alphabet, alphabet)
                # FIXME: Until EddyRivasLab/easel#61 is merged, we cannot
                #        test the case where `guess_alphabet` would throw
                #        an error because of a bug in `sqascii_GuessAlphabet`
                #        causing `eslEOD` to be returned when `eslNOALPHABET`
                #        is expected.

    def test_rewind_fileobj(self):
        # FIXME: Rewinding a `SequenceFile` with an underlying file-like object
        #        is currently unsupported, because Easel will just reopen the
        #        file
        if self.format in easel.MSAFile._FORMATS:
            raise unittest.SkipTest("{!r} format doesn't support rewinding with file-like objects".format(self.format))
        # check reading a file while specifying the format works
        for filename, count, alphabet in zip_longest(self.filenames, self.counts, self.alphabet):
            if alphabet is None:
                continue
            path = os.path.join(self.folder, filename)
            with open(path, "rb") as f:
                buffer = io.BytesIO(f.read())
            with easel.SequenceFile(buffer, self.format) as f:
                seqs1 = f.read_block()
                f.rewind()
                seqs2 = f.read_block()
                self.assertEqual(len(seqs1), count)
                self.assertEqual(len(seqs2), count)
                self.assertEqual(seqs1, seqs2)


@unittest.skipUnless(os.path.exists(EASEL_FOLDER), "test data not available")
class TestEMBLFile(_TestReadFilename, _TestReadFileObject, unittest.TestCase):
    folder    = os.path.join(EASEL_FOLDER, "formats")
    filenames = ["embl"]
    starts    = ["gaattcctga"]
    format    = "embl"
    counts    = [1]
    alphabet  = [easel.Alphabet.dna()]


@unittest.skipUnless(os.path.exists(EASEL_FOLDER), "test data not available")
class TestFastaFile(_TestReadFilename, _TestReadFileObject, unittest.TestCase):
    folder    = os.path.join(EASEL_FOLDER, "formats")
    filenames = ["fasta", "fasta.2", "fasta.odd.1"]
    starts    = ["MEMLPNQTIY", "MEMLPNQTIY", ""]
    format    = "fasta"
    counts    = [2, 2, 2]
    alphabet  = [easel.Alphabet.amino(), easel.Alphabet.amino(), None]


@unittest.skipUnless(os.path.exists(EASEL_FOLDER), "test data not available")
class TestGenbankFile(_TestReadFilename, _TestReadFileObject, unittest.TestCase):
    folder    = os.path.join(EASEL_FOLDER, "formats")
    filenames = ["genbank", "genbank.2"]
    starts    = ["atccacggcc", "atccacggcc"]
    format    = "genbank"
    counts    = [2, 2]
    alphabet  = [easel.Alphabet.dna(), easel.Alphabet.dna()]


@unittest.skipUnless(os.path.exists(EASEL_FOLDER), "test data not available")
class TestUniprotFile(_TestReadFilename, _TestReadFileObject, unittest.TestCase):
    folder    = os.path.join(EASEL_FOLDER, "formats")
    filenames = ["uniprot"]
    starts    = ["MEMLPNQTIY"]
    format    = "uniprot"
    counts    = [1]
    alphabet  = [easel.Alphabet.amino()]


@unittest.skipUnless(os.path.exists(EASEL_FOLDER), "test data not available")
class TestA2MFile(_TestReadFilename, _TestReadFileObject, unittest.TestCase):
    folder    = os.path.join(EASEL_FOLDER, "esl_msa_testfiles", "a2m")
    filenames = ["a2m.good.1", "a2m.good.2"]
    starts    = ["VLSEGEWQLV", "UCCGAUAUAG"]
    format    = "a2m"
    counts    = [4, 5]
    alphabet  = [easel.Alphabet.amino(), easel.Alphabet.rna()]

    # @unittest.expectedFailure
    def test_read_fileobject_guess_format(self):
        # cannot guess format of A2M sequences, they get identified as AFA
        # or FASTA, and then reading fails because sequences do not have the
        # same length
        super().test_read_fileobject_guess_format()


@unittest.skipUnless(os.path.exists(EASEL_FOLDER), "test data not available")
class TestAlignedFastaFile(_TestReadFilename, _TestReadFileObject, unittest.TestCase):
    folder    = os.path.join(EASEL_FOLDER, "esl_msa_testfiles", "afa")
    filenames = ["afa.good.1", "afa.good.2", "afa.good.3"]
    starts    = ["VLSEGEWQLV", "UCCGAUAUAG", "mqifvktltg"]
    format    = "afa"
    counts    = [4, 5, 6]
    alphabet  = [easel.Alphabet.amino(), easel.Alphabet.rna(), easel.Alphabet.amino()]

    @unittest.expectedFailure
    def test_read_filename_guess_format(self):
        super().test_read_filename_guess_format()


@unittest.skipUnless(os.path.exists(EASEL_FOLDER), "test data not available")
class TestClustalFile(_TestReadFilename, _TestReadFileObject, unittest.TestCase):
    folder    = os.path.join(EASEL_FOLDER, "esl_msa_testfiles", "clustal")
    filenames = ["clustal.good.2"]
    starts    = ["UCCGAUAUAG"]
    format    = "clustal"
    counts    = [5]
    alphabet  = [easel.Alphabet.rna()]


@unittest.skipUnless(os.path.exists(EASEL_FOLDER), "test data not available")
class TestClustalLikeFile(_TestReadFilename, _TestReadFileObject, unittest.TestCase):
    folder    = os.path.join(EASEL_FOLDER, "esl_msa_testfiles", "clustal")
    filenames = ["clustal.good.1"]
    starts    = ["VLSEGEWQLV"]
    format    = "clustallike"
    counts    = [4]
    alphabet  = [easel.Alphabet.amino()]


@unittest.skipUnless(os.path.exists(EASEL_FOLDER), "test data not available")
class TestPhylipFile(_TestReadFilename, _TestReadFileObject, unittest.TestCase):
    folder    = os.path.join(EASEL_FOLDER, "esl_msa_testfiles", "phylip")
    filenames = ["phylip.good.1", "phylip.good.2", "phylip.good.3"]
    starts    = ["AAGCTNGGGC", "ATGGCGAAGG", "MKVILLFVLA"]
    format    = "phylip"
    counts    = [5, 7, 3]
    alphabet  = [easel.Alphabet.dna(), easel.Alphabet.dna(), easel.Alphabet.amino()]


@unittest.skipUnless(os.path.exists(EASEL_FOLDER), "test data not available")
class TestPhylipsFile(_TestReadFilename, _TestReadFileObject, unittest.TestCase):
    folder    = os.path.join(EASEL_FOLDER, "esl_msa_testfiles", "phylips")
    filenames = ["phylips.good.1", "phylips.good.2"]
    starts    = ["AAGCTNGGGC", "MKVILLFVLA"]
    format    = "phylips"
    counts    = [5, 3]
    alphabet  = [easel.Alphabet.dna(), easel.Alphabet.amino()]


@unittest.skipUnless(os.path.exists(EASEL_FOLDER), "test data not available")
class TestPsiblastFile(_TestReadFilename, _TestReadFileObject, unittest.TestCase):
    folder    = os.path.join(EASEL_FOLDER, "esl_msa_testfiles", "psiblast")
    filenames = ["psiblast.good.1", "psiblast.good.2"]
    starts    = ["VLSEGEWQLV", "UCCGAUAUAG"]
    format    = "psiblast"
    counts    = [4, 5]
    alphabet  = [easel.Alphabet.amino(), easel.Alphabet.rna()]


@unittest.skipUnless(os.path.exists(EASEL_FOLDER), "test data not available")
class TestSelexFile(_TestReadFilename, _TestReadFileObject, unittest.TestCase):
    folder    = os.path.join(EASEL_FOLDER, "esl_msa_testfiles", "selex")
    filenames = ["selex.good.1", "selex.good.2", "selex.good.3", "selex.good.4"]
    starts    = ["ACDEFGHIKL", "ACDEFGHIKL", "ACDEFGHIKL", "gGAGUAAGAU"]
    format    = "selex"
    counts    = [5, 5, 7, 11]
    alphabet  = [easel.Alphabet.amino(), easel.Alphabet.amino(), easel.Alphabet.amino(), easel.Alphabet.rna()]

    @unittest.expectedFailure
    def test_read_fileobject_given_format(self):
        # unknown error: '#=RF line isn't in expected order in block'
        super().test_read_fileobject_given_format()

    @unittest.expectedFailure
    def test_read_fileobject_guess_format(self):
        # unknown error: '#=RF line isn't in expected order in block'
        super().test_read_fileobject_guess_format()

    @unittest.expectedFailure
    def test_read_fileobject_count_sequences(self):
        # unknown error: '#=RF line isn't in expected order in block'
        super().test_read_fileobject_count_sequences()


@unittest.skipUnless(os.path.exists(EASEL_FOLDER), "test data not available")
class TestStockholmFile(_TestReadFilename, _TestReadFileObject, unittest.TestCase):
    folder    = os.path.join(EASEL_FOLDER, "esl_msa_testfiles", "stockholm")
    filenames = ["stockholm.good.1"]
    starts    = ["ACDEFGHKLM"]
    format    = "stockholm"
    counts    = [7]
    alphabet  = [easel.Alphabet.amino()]
