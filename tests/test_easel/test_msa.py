import copy
import functools
import gc
import io
import os
import unittest
import tempfile

from pyhmmer import easel


class TestMSA(object):
    @classmethod
    def setUpClass(cls):
        cls.formats_folder = os.path.realpath(
            os.path.join(
                __file__, os.pardir, os.pardir, os.pardir, "vendor", "easel", "formats"
            )
        )

    def test_write_roundtrip_stockholm(self):
        sto = os.path.join(self.formats_folder, "stockholm.1")
        msa = self.read_msa(sto)
        with io.BytesIO() as buffer:
            msa.write(buffer, "stockholm")
            actual = buffer.getvalue().decode()
        with open(sto) as f:
            expected = f.read()
        self.assertMultiLineEqual(actual, expected)

    def test_write_invalid_format(self):
        msa = easel.TextMSA()
        self.assertRaises(ValueError, msa.write, io.BytesIO(), "invalidformat")

    def test_sequences(self):
        sto = os.path.join(self.formats_folder, "stockholm.1")
        msa = self.read_msa(sto)
        self.assertEqual(len(msa.sequences), 2)
        self.assertEqual(msa.sequences[0].name, b"seq1")
        self.assertEqual(msa.sequences[1].name, b"seq2")

    def test_init_empty(self):
        msa = self.MSA()
        self.assertEqual(len(msa), 0)
        self.assertEqual(len(msa.sequences), 0)
        self.assertFalse(msa)

    def test_init_name(self):
        msa = self.MSA(name=b"ali")
        self.assertEqual(msa.name, b"ali")

    def test_init_description(self):
        d = b"an alignment made from Python"
        msa = self.MSA(description=d)
        self.assertEqual(msa.description, d)

    def test_init_author(self):
        author = b"Martin Larralde"
        msa = self.MSA(author=author)
        self.assertEqual(msa.author, author)

    def test_init_accession(self):
        acc = b"TST001"
        msa = self.MSA(accession=acc)
        self.assertEqual(msa.accession, acc)


class TestTextMSA(TestMSA, unittest.TestCase):

    MSA = easel.TextMSA

    @staticmethod
    def read_msa(sto):
        with easel.MSAFile(sto, "stockholm") as msa_file:
            return msa_file.read()

    def test_eq(self):
        s1 = easel.TextSequence(name=b"seq1", sequence="ATGC")
        s2 = easel.TextSequence(name=b"seq2", sequence="ATGG")
        msa = easel.TextMSA(sequences=[s1, s2])
        msa2 = easel.TextMSA(sequences=[s1, s2])
        self.assertEqual(msa, msa2)

        msa2.name = b"other"
        self.assertNotEqual(msa, msa2)

        self.assertNotEqual(msa, 1)
        self.assertNotEqual(msa, [s1, s2])

    def test_eq_copy(self):
        s1 = easel.TextSequence(name=b"seq1", sequence="ATGC")
        s2 = easel.TextSequence(name=b"seq2", sequence="ATGG")
        msa = easel.TextMSA(sequences=[s1, s2])
        self.assertEqual(msa, msa.copy())
        self.assertEqual(msa, copy.copy(msa))

    def test_init_sequences(self):
        s1 = easel.TextSequence(name=b"seq1", sequence="ATGC")
        msa = easel.TextMSA(sequences=[s1])
        self.assertEqual(len(msa), 4)
        self.assertEqual(len(msa.sequences), 1)
        self.assertTrue(msa)

    def test_init_length_mismatch(self):
        s1 = easel.TextSequence(name=b"seq1", sequence="ATGC")
        s2 = easel.TextSequence(name=b"seq2", sequence="AT")
        self.assertRaises(ValueError, easel.TextMSA, sequences=[s1, s2])

    def test_init_duplicate_names(self):
        s1 = easel.TextSequence(name=b"seq1", sequence="ATGC")
        s2 = easel.TextSequence(name=b"seq1", sequence="ATTC")
        self.assertRaises(ValueError, easel.TextMSA, sequences=[s1, s2])


class TestDigitalMSA(TestMSA, unittest.TestCase):

    alphabet = easel.Alphabet.dna()
    MSA = functools.partial(easel.DigitalMSA, alphabet)

    @staticmethod
    def read_msa(sto):
        with easel.MSAFile(sto, "stockholm") as msa_file:
            msa_file.set_digital(msa_file.guess_alphabet())
            return msa_file.read()

    def test_eq(self):
        s1 = easel.DigitalSequence(self.alphabet, name=b"seq1", sequence=bytearray([1, 2, 3, 4]))
        s2 = easel.DigitalSequence(self.alphabet, name=b"seq2", sequence=bytearray([1, 2, 3, 3]))
        msa = easel.DigitalMSA(self.alphabet, sequences=[s1, s2])
        msa2 = easel.DigitalMSA(self.alphabet, sequences=[s1, s2])
        self.assertEqual(msa, msa2)

        msa2.name = b"other"
        self.assertNotEqual(msa, msa2)

        self.assertNotEqual(msa, 1)
        self.assertNotEqual(msa, [s1, s2])

    def test_eq_copy(self):
        s1 = easel.DigitalSequence(self.alphabet, name=b"seq1", sequence=bytearray([1, 2, 3, 4]))
        s2 = easel.DigitalSequence(self.alphabet, name=b"seq2", sequence=bytearray([1, 2, 3, 3]))
        msa = easel.DigitalMSA(self.alphabet, sequences=[s1, s2])
        self.assertEqual(msa, msa.copy())
        self.assertEqual(msa, copy.copy(msa))

    def test_init_sequences(self):
        s1 = easel.DigitalSequence(self.alphabet, name=b"seq1", sequence=bytearray([1, 2, 3, 4]))
        msa = easel.DigitalMSA(self.alphabet, sequences=[s1])
        self.assertEqual(len(msa), 4)
        self.assertEqual(len(msa.sequences), 1)
        self.assertTrue(msa)

    def test_init_length_mismatch(self):
        s1 = easel.DigitalSequence(self.alphabet, name=b"seq1", sequence=bytearray([1, 2, 3, 4]))
        s2 = easel.DigitalSequence(self.alphabet, name=b"seq2", sequence=bytearray([1, 2]))
        self.assertRaises(ValueError, easel.DigitalMSA, self.alphabet, sequences=[s1, s2])

    def test_init_duplicate_names(self):
        s1 = easel.DigitalSequence(self.alphabet, name=b"seq1", sequence=bytearray([1, 2, 3, 4]))
        s2 = easel.DigitalSequence(self.alphabet, name=b"seq1", sequence=bytearray([1, 2, 3, 3]))
        self.assertRaises(ValueError, easel.DigitalMSA, self.alphabet, sequences=[s1, s2])
