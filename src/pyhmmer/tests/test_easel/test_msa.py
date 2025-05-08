import copy
import functools
import gc
import io
import os
import unittest
import tempfile

from pyhmmer import easel

from .. import __name__ as __package__
from .utils import EASEL_FOLDER, resource_files


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

    def test_indexed_key_error(self):
        msa = self.MSA()
        with self.assertRaises(KeyError):
            msa.indexed[b"xxx"]

    def test_indexed_type_error(self):
        msa = self.MSA()
        with self.assertRaises(TypeError):
            msa.indexed[1]


class TestTextMSA(TestMSA, unittest.TestCase):

    MSA = easel.TextMSA

    @staticmethod
    def read_msa(sto):
        with easel.MSAFile(sto, "stockholm") as msa_file:
            return msa_file.read()

    def test_indexed(self):
        s1 = easel.TextSequence(name=b"seq1", sequence="ATGC")
        s2 = easel.TextSequence(name=b"seq2", sequence="ATGG")
        msa = easel.TextMSA(sequences=[s1, s2])
        self.assertEqual(len(msa.indexed), 2)
        self.assertEqual(msa.indexed[b"seq1"], s1)
        self.assertEqual(msa.indexed[b"seq2"], s2)

    def test_rf(self):
        s1 = easel.TextSequence(name=b"seq1", sequence="ATGC")
        s2 = easel.TextSequence(name=b"seq2", sequence="ATGG")
        msa = easel.TextMSA(sequences=[s1, s2])
        self.assertIsNone(msa.reference)
        msa.reference = b"xxxx"
        self.assertEqual(msa.reference, b"xxxx")
        with self.assertRaises(ValueError):
            msa.reference = b"x"

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
    MSA = staticmethod(functools.partial(easel.DigitalMSA, alphabet))

    @staticmethod
    def read_msa(sto):
        with easel.MSAFile(sto, "stockholm", digital=True) as msa_file:
            return msa_file.read()

    def test_sample(self):
        rng = easel.Randomness(42)
        msa = easel.DigitalMSA.sample(self.alphabet, 10, 50, rng)
        self.assertLessEqual(len(msa.sequences), 10)
        self.assertLessEqual(len(msa), 50)

    def test_sample_seed(self):
        msa1 = easel.DigitalMSA.sample(self.alphabet, 10, 50, randomness=easel.Randomness(42))
        msa2 = easel.DigitalMSA.sample(self.alphabet, 10, 50, randomness=42)
        self.assertEqual(msa1, msa2)
        msa3 = easel.DigitalMSA.sample(self.alphabet, 10, 50, randomness=100)
        self.assertNotEqual(msa3, msa2)

        with self.assertRaises(TypeError):
            msa4 = easel.DigitalMSA.sample(self.alphabet, 10, 50, randomness="xxx")

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


    def test_identity_filter(self):
        # adapted from `utest_idfilter` in `esl_msaweight.c`
        dna = easel.Alphabet.dna()
        s1 = easel.DigitalSequence

        buffer = io.BytesIO(
            b"# STOCKHOLM 1.0\n\n"
            b"seq1  ..AAAAAAAA\n"
            b"seq2  AAAAAAAAAA\n"
            b"seq3  CCCCCCCCCC\n"
            b"seq4  GGGGGGGGGG\n"
            b"//\n"
        )

        with easel.MSAFile(buffer, format="stockholm", digital=True, alphabet=dna) as msa_file:
            msa = msa_file.read()

        msa2 = msa.identity_filter(1.0)
        self.assertEqual(len(msa2.sequences), 3)
        self.assertEqual(msa2.names[0], b"seq2")

        msa3 = msa.identity_filter(1.0, preference="origorder")
        self.assertEqual(len(msa3.sequences), 3)
        self.assertEqual(msa3.names[0], b"seq1")

        msa3 = msa.identity_filter(1.0, preference="random")
        self.assertEqual(len(msa3.sequences), 3)
        self.assertIn(msa3.names[0], [b"seq1", b"seq2"])

    @unittest.skipUnless(resource_files, "importlib.resources.files not available")
    def test_identity_filter_luxc(self):
        luxc = resource_files(__package__).joinpath("data", "msa", "LuxC.sto")
        msa = self.read_msa(luxc)
        filtered = msa.identity_filter()

        weighted = resource_files(__package__).joinpath("data", "msa", "LuxC.weighted.sto")
        expected = self.read_msa(weighted)

        self.assertEqual(len(filtered.sequences), len(expected.sequences))
        self.assertEqual(len(filtered.names), len(expected.names))
        for s1, s2 in zip(filtered.sequences, expected.sequences):
            self.assertEqual(s1, s2)
