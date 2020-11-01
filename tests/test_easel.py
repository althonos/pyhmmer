import os
import unittest
import tempfile

from hmmer import easel


class TestSequenceFile(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.formats_folder = os.path.realpath(os.path.join(
            __file__, os.pardir, os.pardir, "vendor", "easel", "formats"
        ))

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

    def test_read_fasta(self):
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

    def test_read_embl(self):
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

    def test_read_fasta_error(self):
        fasta = os.path.join(self.formats_folder, "fasta.bad.1")
        with easel.SequenceFile(fasta, "fasta") as f:
            self.assertRaises(ValueError, f.read)
