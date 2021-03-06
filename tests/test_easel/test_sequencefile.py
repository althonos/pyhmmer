import copy
import gc
import io
import os
import unittest
import tempfile

from pyhmmer import easel


class TestSequenceFile(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.formats_folder = os.path.realpath(
            os.path.join(
                __file__, os.pardir, os.pardir, os.pardir, "vendor", "easel", "formats"
            )
        )

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

    def test_readformat_fasta(self):
        for filename in ["fasta", "fasta.2", "fasta.odd.1"]:
            fasta = os.path.join(self.formats_folder, filename)

            # check reading a FASTA file works
            with easel.SequenceFile(fasta, "fasta") as f:
                seq = f.read()
                self.assertNotEqual(seq.name, b"")

            # check reading a FASTA file without specifying the format works too
            with easel.SequenceFile(fasta) as f:
                seq2 = f.read()
                self.assertNotEqual(seq2.name, b"")

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
            self.assertNotEqual(seq.accession, b"")

        # check reading an EMBL file without specifying the format works too
        with easel.SequenceFile(embl) as f:
            seq2 = f.read()
            self.assertNotEqual(seq2.accession, b"")

        # check reading an EMBL file while giving another format fails
        with easel.SequenceFile(embl, "fasta") as f:
            self.assertRaises(ValueError, f.read)

    def test_readformat_genbank(self):
        gbk = os.path.join(self.formats_folder, "genbank")

        # check reading an EMBL file works
        with easel.SequenceFile(gbk, "genbank") as f:
            seq = f.read()
            self.assertNotEqual(seq.accession, b"")

        # check reading an GenBank file without specifying the format works too
        with easel.SequenceFile(gbk) as f:
            seq2 = f.read()
            self.assertNotEqual(seq2.accession, b"")

        # check reading an GenBank file while giving another format fails
        with easel.SequenceFile(gbk, "fasta") as f:
            self.assertRaises(ValueError, f.read)
