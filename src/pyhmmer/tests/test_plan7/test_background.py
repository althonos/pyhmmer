import copy
import gc
import io
import os
import pickle
import tempfile
import unittest

from pyhmmer import easel, plan7


class TestBuilder(unittest.TestCase):

    def test_init(self):
        abc = easel.Alphabet.dna()
        bg = plan7.Background(abc)
        self.assertFalse(bg.uniform)

        bgu = plan7.Background(abc, uniform=True)
        self.assertTrue(bgu.uniform)

    def test_copy(self):
        abc = easel.Alphabet.dna()
        bg = plan7.Background(abc)

        bg_copy = copy.copy(bg)
        self.assertIs(bg.alphabet, bg_copy.alphabet)
        self.assertEqual(bg.uniform, bg_copy.uniform)
        self.assertEqual(bg.L, bg_copy.L)
        self.assertEqual(bg.transition_probability, bg_copy.transition_probability)
        self.assertEqual(bg.omega, bg_copy.omega)
        self.assertEqual(list(bg.residue_frequencies), list(bg_copy.residue_frequencies))

    def test_pickle(self):
        abc = easel.Alphabet.dna()
        bg = plan7.Background(abc)

        bg2 = pickle.loads(pickle.dumps(bg))
        self.assertEqual(bg.alphabet, bg2.alphabet)
        self.assertEqual(bg.transition_probability, bg2.transition_probability)
        self.assertEqual(bg.residue_frequencies, bg2.residue_frequencies)

    def test_repr(self):
        abc = easel.Alphabet.amino()
        bg = plan7.Background(abc)
        self.assertEqual(repr(bg), "Background(Alphabet.amino(), uniform=False)")
        bg2 = plan7.Background(abc, uniform=True)
        self.assertEqual(repr(bg2), "Background(Alphabet.amino(), uniform=True)")
