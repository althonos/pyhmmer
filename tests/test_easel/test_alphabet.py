import unittest

from pyhmmer import easel


class TestAlphabet(unittest.TestCase):
    def test_eq(self):
        self.assertEqual(easel.Alphabet.dna(), easel.Alphabet.dna())
        self.assertEqual(easel.Alphabet.rna(), easel.Alphabet.rna())
        self.assertEqual(easel.Alphabet.amino(), easel.Alphabet.amino())
        self.assertNotEqual(easel.Alphabet.amino(), easel.Alphabet.dna())
        self.assertNotEqual(easel.Alphabet.amino(), object())
