import copy
import gc
import io
import os
import unittest
import tempfile

from pyhmmer import easel



class _TestSequenceFileBase(object):

    @classmethod
    def setUpClass(cls):
        cls.formats_folder = os.path.realpath(
            os.path.join(
                __file__, os.pardir, os.pardir, os.pardir, "vendor", "easel", "formats"
            )
        )


class TestSequenceFile(_TestSequenceFileBase, unittest.TestCase):

    def test_guess_alphabet(self):
        embl = os.path.join(self.formats_folder, "embl")
        with easel.SequenceFile(embl, "embl") as f:
            self.assertEqual(f.guess_alphabet(), easel.Alphabet.dna())

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

    def test_iter(self):
        fasta = os.path.join(self.formats_folder, "fasta")

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

    def test_read_skip_info(self):
        fasta = os.path.join(self.formats_folder, "fasta")

        with easel.SequenceFile(fasta) as f:
            seq = f.read()
            self.assertEqual(seq.name, b"SNRPA_DROME")
            self.assertNotEqual(seq.sequence, "")

        with easel.SequenceFile(fasta) as f:
            seq = f.read(skip_info=True)
            self.assertEqual(seq.name, b"")
            self.assertNotEqual(seq.sequence, "")

    def test_read_skip_sequence(self):
        fasta = os.path.join(self.formats_folder, "fasta")

        with easel.SequenceFile(fasta) as f:
            seq = f.read()
            self.assertEqual(seq.name, b"SNRPA_DROME")
            self.assertNotEqual(seq.sequence, "")

        with easel.SequenceFile(fasta) as f:
            seq = f.read(skip_sequence=True)
            self.assertEqual(seq.name, b"SNRPA_DROME")
            self.assertEqual(seq.sequence, "")

    def test_ignore_gaps(self):
        luxc = os.path.realpath(
            os.path.join(__file__, os.pardir, os.pardir, "data", "msa", "LuxC.faa")
        )
        # fails if not ignoring gaps (since if contains gaps)
        with easel.SequenceFile(luxc, "fasta") as seq_file:
            self.assertRaises(ValueError, seq_file.read)
        # succeeds if ignoring gaps
        with easel.SequenceFile(luxc, "fasta", ignore_gaps=True) as seq_file:
            sequences = list(seq_file)
            self.assertEqual(len(sequences), 13)


class _TestReadFilename(object):

    def test_read_filename_guess_format(self):
        # check reading a file without specifying the format works
        for filename, start in zip(self.filenames, self.starts):
            path = os.path.join(self.formats_folder, filename)
            with easel.SequenceFile(path) as f:
                seq = f.read()
                self.assertEqual(seq.sequence[:10], start)

    def test_read_filename_given_format(self):
        # check reading a file while specifying the format works
        for filename, start in zip(self.filenames, self.starts):
            path = os.path.join(self.formats_folder, filename)
            with easel.SequenceFile(path, self.format) as f:
                seq = f.read()
                self.assertEqual(seq.sequence[:10], start)

    def test_read_filename_count_sequences(self):
        # check reading a file while specifying the format works
        for filename, count in zip(self.filenames, self.counts):
            path = os.path.join(self.formats_folder, filename)
            with easel.SequenceFile(path, self.format) as f:
                sequences = list(f)
                self.assertEqual(len(sequences), count)

    def test_guess_alphabet(self):
        for filename in self.filenames:
            path = os.path.join(self.formats_folder, filename)
            with easel.SequenceFile(path, self.format) as f:
                alphabet = f.guess_alphabet()
                self.assertEqual(alphabet, self.alphabet)


class _TestReadFileObject(object):

    def test_raise_errors(self):

        class WeirdError(ValueError):
            pass

        class FailibleReader(object):
            def __init__(self, handle, fault=0):
                self.handle = handle
            def read(self, n=-1):
                raise WeirdError("oops")
            def readable(self):
                return True
            def seekable(self):
                return False
            def writable(self):
                return False

        path = os.path.join(self.formats_folder, self.filenames[0])
        with open(path, "rb") as f:
            buffer = io.BytesIO(f.read())

        # this should fail in the constructor if `buffer.read` immediately fails
        self.assertRaises(WeirdError, easel.SequenceFile, FailibleReader(buffer, 0), format=self.format)

    def test_read_fileobject_given_format(self):
        # check reading a file while specifying the format works
        for filename, start in zip(self.filenames, self.starts):
            path = os.path.join(self.formats_folder, filename)
            with open(path, "rb") as f:
                buffer = io.BytesIO(f.read())
            with easel.SequenceFile(buffer, self.format) as f:
                seq = f.read()
                self.assertEqual(seq.sequence[:10], start)

    def test_read_filename_given_format(self):
        # check reading a file while specifying the format works
        for filename, start in zip(self.filenames, self.starts):
            path = os.path.join(self.formats_folder, filename)
            with open(path, "rb") as f:
                buffer = io.BytesIO(f.read())
            with easel.SequenceFile(buffer, self.format) as f:
                seq = f.read()
                self.assertEqual(seq.sequence[:10], start)

    def test_read_filename_count_sequences(self):
        # check reading a file while specifying the format works
        for filename, count in zip(self.filenames, self.counts):
            path = os.path.join(self.formats_folder, filename)
            with open(path, "rb") as f:
                buffer = io.BytesIO(f.read())
            with easel.SequenceFile(buffer, self.format) as f:
                sequences = list(f)
                self.assertEqual(len(sequences), count)


class TestFastaFile(_TestSequenceFileBase, _TestReadFilename, _TestReadFileObject, unittest.TestCase):
    filenames = ["fasta", "fasta.2"]
    starts    = ["MEMLPNQTIY", "MEMLPNQTIY"]
    format    = "fasta"
    counts    = [2]
    alphabet  = easel.Alphabet.amino()


class TestGenbankFile(_TestSequenceFileBase, _TestReadFilename, _TestReadFileObject, unittest.TestCase):
    filenames = ["genbank", "genbank.2"]
    starts    = ["atccacggcc", "atccacggcc"]
    format    = "genbank"
    counts    = [2]
    alphabet  = easel.Alphabet.dna()


class TestUniprotFile(_TestSequenceFileBase, _TestReadFilename, _TestReadFileObject, unittest.TestCase):
    filenames = ["uniprot"]
    starts    = ["MEMLPNQTIY"]
    format    = "uniprot"
    counts    = [1]
    alphabet  = easel.Alphabet.amino()
