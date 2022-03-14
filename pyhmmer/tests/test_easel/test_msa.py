import copy
import functools
import gc
import io
import os
import unittest
import tempfile

from pyhmmer import easel

from ..utils import EASEL_FOLDER


class TestMSA(object):
    @classmethod
    def setUpClass(cls):
        cls.formats_folder = os.path.join(EASEL_FOLDER, "formats")

    @unittest.skipUnless(os.path.exists(EASEL_FOLDER), "test data not available")
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

    @unittest.skipUnless(os.path.exists(EASEL_FOLDER), "test data not available")
    def test_sequences(self):
        sto = os.path.join(self.formats_folder, "stockholm.1")
        msa = self.read_msa(sto)
        self.assertEqual(len(msa.sequences), 2)
        self.assertEqual(msa.sequences[0].name, b"seq1")
        self.assertEqual(msa.sequences[1].name, b"seq2")
        self.assertEqual(msa.sequences[-2].name, b"seq1")
        self.assertEqual(msa.sequences[-1].name, b"seq2")

        with self.assertRaises(IndexError):
            msa.sequences[-3]

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

    def test_init_error(self):
        self.assertRaises(TypeError, easel.TextMSA, sequences=[1, 2, 3])

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

    def test_insert_sequences(self):
        s1 = easel.TextSequence(name=b"seq1", sequence="ATGC", description=b"first", accession=b"SQ.1")
        s2 = easel.TextSequence(name=b"seq2", sequence="ATTC", description=b"second", accession=b"SQ.2")
        msa = easel.TextMSA(sequences=[s1, s2])
        self.assertEqual(len(msa), 4)
        self.assertEqual(len(msa.sequences), 2)
        self.assertEqual(msa.sequences[0].name, b"seq1")
        self.assertEqual(msa.sequences[1].name, b"seq2")

        s3 = easel.TextSequence(name=b"seq3", sequence="AGGC")
        msa.sequences[0] = s3
        self.assertEqual(len(msa), 4)
        self.assertEqual(len(msa.sequences), 2)
        self.assertEqual(msa.sequences[0].name, b"seq3")
        self.assertEqual(msa.sequences[1].name, b"seq2")

        s4 = easel.TextSequence(name=b"seq4", sequence="TTTT")
        msa.sequences[-1] = s4
        self.assertEqual(len(msa), 4)
        self.assertEqual(len(msa.sequences), 2)
        self.assertEqual(msa.sequences[0].name, b"seq3")
        self.assertEqual(msa.sequences[1].name, b"seq4")

        s5 = easel.TextSequence(name=b"seq5", sequence="AAAA")
        with self.assertRaises(IndexError):
            msa.sequences[2] = s5
        with self.assertRaises(IndexError):
            msa.sequences[-10] = s5

        s4bis = easel.TextSequence(name=b"seq4", sequence="TTTT")
        with self.assertRaises(ValueError):
            msa.sequences[0] = s4bis

        s6 = easel.TextSequence(name=b"seq4", sequence="TT")
        with self.assertRaises(ValueError):
            msa.sequences[0] = s6

        dna = easel.Alphabet.dna()
        s7 = easel.DigitalSequence(dna, name=b"seq4", sequence=bytearray([2, 2]))
        with self.assertRaises(TypeError):
            msa.sequences[0] = s7

    def test_digitize(self):
        dna = easel.Alphabet.dna()
        s1 = easel.TextSequence(name=b"seq1", sequence="ACGT")
        msa_t = easel.TextMSA(sequences=[s1])
        msa_d = msa_t.digitize(dna)

        self.assertIsInstance(msa_d, easel.DigitalMSA)
        self.assertEqual(len(msa_d), 4)
        self.assertEqual(len(msa_d.sequences), 1)
        self.assertTrue(msa_d)
        self.assertEqual(msa_d.sequences[0].name, b"seq1")
        self.assertEqual(msa_d.sequences[0].sequence, easel.VectorU8([0, 1, 2, 3]))

    def test_digitize_invalid(self):
        dna = easel.Alphabet.dna()
        s1 = easel.TextSequence(name=b"seq1", sequence=">v3ry 1nv4l1d")
        msa_t = easel.TextMSA(sequences=[s1])
        self.assertRaises(ValueError, msa_t.digitize, dna)


class TestDigitalMSA(TestMSA, unittest.TestCase):

    alphabet = easel.Alphabet.dna()
    MSA = functools.partial(easel.DigitalMSA, alphabet)

    @staticmethod
    def read_msa(sto):
        with easel.MSAFile(sto, "stockholm", digital=True) as msa_file:
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

    def test_init_error(self):
        self.assertRaises(TypeError, easel.TextMSA, sequences=[1, 2, 3])

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

    def test_insert_sequences(self):
        s1 = easel.DigitalSequence(self.alphabet, name=b"seq1", sequence=bytearray([1, 2, 3, 4]), description=b"first", accession=b"SQ.1")
        s2 = easel.DigitalSequence(self.alphabet, name=b"seq2", sequence=bytearray([1, 3, 3, 4]), description=b"second", accession=b"SQ.2")
        msa = easel.DigitalMSA(self.alphabet, sequences=[s1, s2])
        self.assertEqual(len(msa), 4)
        self.assertEqual(len(msa.sequences), 2)
        self.assertEqual(msa.sequences[0].name, b"seq1")
        self.assertEqual(msa.sequences[1].name, b"seq2")

        s3 = easel.DigitalSequence(self.alphabet, name=b"seq3", sequence=bytearray([4, 2, 3, 4]))
        msa.sequences[0] = s3
        self.assertEqual(len(msa), 4)
        self.assertEqual(len(msa.sequences), 2)
        self.assertEqual(msa.sequences[0].name, b"seq3")
        self.assertEqual(msa.sequences[1].name, b"seq2")

        s4 = easel.DigitalSequence(self.alphabet, name=b"seq4", sequence=bytearray([2, 2, 2, 2]))
        msa.sequences[-1] = s4
        self.assertEqual(len(msa), 4)
        self.assertEqual(len(msa.sequences), 2)
        self.assertEqual(msa.sequences[0].name, b"seq3")
        self.assertEqual(msa.sequences[1].name, b"seq4")

        s5 = easel.DigitalSequence(self.alphabet, name=b"seq5", sequence=bytearray([1, 1, 1, 1]))
        with self.assertRaises(IndexError):
            msa.sequences[2] = s5
        with self.assertRaises(IndexError):
            msa.sequences[-10] = s5

        s4bis = easel.DigitalSequence(self.alphabet, name=b"seq4", sequence=bytearray([2, 2, 2, 2]))
        with self.assertRaises(ValueError):
            msa.sequences[0] = s4bis

        s6 = easel.DigitalSequence(self.alphabet, name=b"seq4", sequence=bytearray([2, 2]))
        with self.assertRaises(ValueError):
            msa.sequences[0] = s6

        s7 = easel.TextSequence(name=b"seq4", sequence="ATGC")
        with self.assertRaises(TypeError):
            msa.sequences[0] = s7

    def test_textize(self):
        dna = easel.Alphabet.dna()
        s1 = easel.DigitalSequence(self.alphabet, name=b"seq1", sequence=bytearray([0, 1, 2, 3]))
        msa_d = easel.DigitalMSA(self.alphabet, sequences=[s1])
        msa_t = msa_d.textize()

        self.assertIsInstance(msa_t, easel.TextMSA)
        self.assertEqual(len(msa_t), 4)
        self.assertEqual(len(msa_t.sequences), 1)
        self.assertTrue(msa_t)
        self.assertEqual(msa_t.sequences[0].name, b"seq1")
        self.assertEqual(msa_t.sequences[0].sequence, "ACGT")
