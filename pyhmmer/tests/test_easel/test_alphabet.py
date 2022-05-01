import pickle
import unittest
import sys

from pyhmmer import easel


class TestAlphabet(unittest.TestCase):

    @unittest.skipIf(sys.implementation.name == "pypy", "`getsizeof` not supported on PyPY")
    def test_sizeof(self):
        alphabet = easel.Alphabet.dna()
        self.assertGreater(sys.getsizeof(alphabet), 0)

    def test_is_nucleotide(self):
        self.assertTrue(easel.Alphabet.dna().is_nucleotide())
        self.assertTrue(easel.Alphabet.rna().is_nucleotide())
        self.assertFalse(easel.Alphabet.amino().is_nucleotide())

    def test_is_dna(self):
        self.assertTrue(easel.Alphabet.dna().is_dna())
        self.assertFalse(easel.Alphabet.rna().is_dna())
        self.assertFalse(easel.Alphabet.amino().is_dna())

    def test_is_rna(self):
        self.assertTrue(easel.Alphabet.rna().is_rna())
        self.assertFalse(easel.Alphabet.dna().is_rna())
        self.assertFalse(easel.Alphabet.amino().is_rna())

    def test_is_amino(self):
        self.assertTrue(easel.Alphabet.amino().is_amino())
        self.assertFalse(easel.Alphabet.rna().is_amino())
        self.assertFalse(easel.Alphabet.dna().is_amino())

    def test_eq(self):
        self.assertEqual(easel.Alphabet.dna(), easel.Alphabet.dna())
        self.assertEqual(easel.Alphabet.rna(), easel.Alphabet.rna())
        self.assertEqual(easel.Alphabet.amino(), easel.Alphabet.amino())
        self.assertNotEqual(easel.Alphabet.amino(), easel.Alphabet.dna())
        self.assertNotEqual(easel.Alphabet.amino(), object())

    def test_repr(self):
        self.assertEqual(repr(easel.Alphabet.amino()), "pyhmmer.easel.Alphabet.amino()")
        self.assertEqual(repr(easel.Alphabet.rna()), "pyhmmer.easel.Alphabet.rna()")
        self.assertEqual(repr(easel.Alphabet.dna()), "pyhmmer.easel.Alphabet.dna()")

    def test_K(self):
        self.assertEqual(easel.Alphabet.dna().K, 4)
        self.assertEqual(easel.Alphabet.rna().K, 4)
        self.assertEqual(easel.Alphabet.amino().K, 20)

    def test_picklable(self):
        abc = easel.Alphabet.dna()
        p = pickle.dumps(abc)
        abc2 = pickle.loads(p)
        self.assertIsNot(abc2, abc)
        self.assertEqual(abc, abc)
