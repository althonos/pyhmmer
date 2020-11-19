import copy
import gc
import os
import unittest
import tempfile

from pyhmmer import easel


class TestSequenceFile(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.formats_folder = os.path.realpath(os.path.join(
            __file__, os.pardir, os.pardir, "vendor", "easel", "formats"
        ))

    def test_init_error_unknownformat(self):
        with self.assertRaises(ValueError):
            _file = easel.SequenceFile("file.x", format="nonsense")

    def test_init_error_empty(self):
        with tempfile.NamedTemporaryFile() as empty:
            # without a format argument, we can't determine the file format
            self.assertRaises(ValueError, easel.SequenceFile, empty.name)
            # with a format argument, we don't expect the file to be empty
            self.assertRaises(EOFError, easel.SequenceFile, empty.name, "fasta")

    def test_init_error_filenotfound(self):
        self.assertRaises(FileNotFoundError, easel.SequenceFile, "path/to/missing/file.fa")

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

        with easel.SequenceFile(fasta) as f:
            seq = f.read(skip_info=True)
            self.assertEqual(seq.name, b"")

    def test_readformat_fasta(self):
        for filename in ["fasta", "fasta.2", "fasta.odd.1"]:
            fasta = os.path.join(self.formats_folder, filename)

            # check reading a FASTA file works
            with easel.SequenceFile(fasta, "fasta") as f:
                seq = f.read()
                self.assertTrue(seq.name)

            # check reading a FASTA file without specifying the format works too
            with easel.SequenceFile(fasta) as f:
                seq2 = f.read()
                self.assertTrue(seq2.name)

            # check reading a FASTA file while giving another format fails
            with easel.SequenceFile(fasta, "embl") as f:
                self.assertRaises(ValueError, f.read)

    def test_readformat_fasta_error(self):
        fasta = os.path.join(self.formats_folder, "fasta.bad.1")
        with easel.SequenceFile(fasta, "fasta") as f:
            self.assertRaises(ValueError, f.read)

    def test_readformat_embl(self):
        embl = os.path.join(self.formats_folder, "embl")

        # check reading an EMBL file works
        with easel.SequenceFile(embl, "embl") as f:
            seq = f.read()
            self.assertTrue(seq.accession)

        # check reading an EMBL file without specifying the format works too
        with easel.SequenceFile(embl) as f:
            seq2 = f.read()
            self.assertTrue(seq2.accession)

        # check reading an EMBL file while giving another format fails
        with easel.SequenceFile(embl, "fasta") as f:
            self.assertRaises(ValueError, f.read)


class TestTextSequence(unittest.TestCase):

    def test_copy(self):
        seq = easel.TextSequence(name=b"TEST", sequence=b"ATGC")

        # copy the sequence
        cpy = seq.copy()

        # check the copy is equal and attributes are equal
        self.assertEqual(seq, cpy)
        self.assertEqual(seq.name, cpy.name)
        self.assertEqual(seq.checksum(), cpy.checksum())

        # check attributes can be read even after the original object
        # is (hopefully) deallocated, to make sure internal strings were copied
        del seq
        gc.collect()
        self.assertEqual(cpy.name, b"TEST")

        # check `copy.copy` works too
        cpy2 = copy.copy(cpy)
        self.assertEqual(cpy, cpy2)
        self.assertEqual(cpy.name, cpy2.name)
