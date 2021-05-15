import pickle
import unittest

from pyhmmer import easel


class TestAlphabet(unittest.TestCase):

    def test_eq(self):
        self.assertEqual(easel.Alphabet.dna(), easel.Alphabet.dna())
        self.assertEqual(easel.Alphabet.rna(), easel.Alphabet.rna())
        self.assertEqual(easel.Alphabet.amino(), easel.Alphabet.amino())
        self.assertNotEqual(easel.Alphabet.amino(), easel.Alphabet.dna())
        self.assertNotEqual(easel.Alphabet.amino(), object())

    def test_repr(self):
        self.assertEqual(repr(easel.Alphabet.amino()), "Alphabet.amino()")
        self.assertEqual(repr(easel.Alphabet.rna()), "Alphabet.rna()")
        self.assertEqual(repr(easel.Alphabet.dna()), "Alphabet.dna()")

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
