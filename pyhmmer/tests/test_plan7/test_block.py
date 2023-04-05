import abc
import unittest

from ...errors import AlphabetMismatch
from ...easel import Alphabet, DigitalSequence, Randomness
from ...plan7 import Background, HMM, Profile, OptimizedProfile, OptimizedProfileBlock


class TestOptimizedProfileBlock(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.rng = Randomness(seed=0)
        cls.alphabet = Alphabet.amino()
        cls.background = Background(cls.alphabet)

    @classmethod
    def _random_hmm(cls, name, M=100):
        hmm = HMM.sample(cls.alphabet, M, cls.rng)
        hmm.name = name
        return hmm

    @classmethod
    def _random_profile(cls, name, M=100):
        hmm = cls._random_hmm(name, M=M)
        profile = Profile(hmm.M, hmm.alphabet)
        profile.configure(hmm, cls.background, 200)
        return profile

    @classmethod
    def _random_optimized_profile(cls, name, M=100):
        return cls._random_profile(name, M=M).to_optimized()

    def test_alphabet_mismatch(self):
        om1 = self._random_optimized_profile(b"profile1")
        dna = Alphabet.dna()
        block = OptimizedProfileBlock(dna)
        self.assertRaises(AlphabetMismatch, block.append, om1)
        self.assertRaises(AlphabetMismatch, block.insert, 0, om1)

    def test_truth(self):
        om1 = self._random_optimized_profile(b"profile1")
        om2 = self._random_optimized_profile(b"profile2")

        self.assertTrue(not OptimizedProfileBlock(self.alphabet))
        block = OptimizedProfileBlock(self.alphabet, [om1, om2])
        self.assertTrue(block)

    def test_identity(self):
        self.assertIsNot(OptimizedProfileBlock(self.alphabet), OptimizedProfileBlock(self.alphabet))
    
    def test_len(self):
        om1 = self._random_optimized_profile(b"profile1")
        om2 = self._random_optimized_profile(b"profile2")
        om3 = self._random_optimized_profile(b"profile3")

        self.assertEqual(len(OptimizedProfileBlock(self.alphabet, [])), 0)
        self.assertEqual(len(OptimizedProfileBlock(self.alphabet, [om1])), 1)
        self.assertEqual(len(OptimizedProfileBlock(self.alphabet, [om1, om2, om3])), 3)

    def test_append(self):
        om1 = self._random_optimized_profile(b"profile1")
        om2 = self._random_optimized_profile(b"profile2")

        block = OptimizedProfileBlock(self.alphabet)
        self.assertEqual(len(block), 0)

        block.append(om1)
        self.assertEqual(len(block), 1)
        self.assertEqual(block[0], om1)

        block.append(om2)
        self.assertEqual(len(block), 2)
        self.assertEqual(block[1], om2)

    def test_iter(self):
        om1 = self._random_optimized_profile(b"profile1")
        om2 = self._random_optimized_profile(b"profile2")
        om3 = self._random_optimized_profile(b"profile3")

        block = OptimizedProfileBlock(self.alphabet, [om1, om2, om3])
        self.assertEqual(len(block), 3)

        iterator = iter(block)
        self.assertIs(next(iterator), om1)
        self.assertIs(next(iterator), om2)
        self.assertIs(next(iterator), om3)
        self.assertRaises(StopIteration, next, iterator)

    def test_clear(self):
        om1 = self._random_optimized_profile(b"profile1")
        om2 = self._random_optimized_profile(b"profile2")

        block = OptimizedProfileBlock(self.alphabet, [om1, om2])
        self.assertEqual(len(block), 2)
        block.clear()
        self.assertEqual(len(block), 0)
        block.clear()
        self.assertEqual(len(block), 0)
        
        block = OptimizedProfileBlock(self.alphabet)
        self.assertEqual(len(block), 0)
        block.clear()
        self.assertEqual(len(block), 0)

    def test_getitem(self):
        om1 = self._random_optimized_profile(b"profile1")
        om2 = self._random_optimized_profile(b"profile2")

        block = OptimizedProfileBlock(self.alphabet, [om1, om2])
        self.assertEqual(len(block), 2)

        self.assertIs(block[0], om1)
        self.assertIs(block[1], om2)
        self.assertIs(block[-1], om2)
        self.assertIs(block[-2], om1)

        with self.assertRaises(IndexError):
            block[3]
        with self.assertRaises(IndexError):
            block[-5]
        
    def test_setitem(self):
        om1 = self._random_optimized_profile(b"profile1")
        om2 = self._random_optimized_profile(b"profile2")
        om3 = self._random_optimized_profile(b"profile3")

        block = OptimizedProfileBlock(self.alphabet, [om1, om2])
        self.assertEqual(len(block), 2)

        block[0] = om3
        self.assertIs(block[0], om3)
        self.assertIs(block[1], om2)

    def test_remove(self):
        om1 = self._random_optimized_profile(b"profile1")
        om2 = self._random_optimized_profile(b"profile2")

        block = OptimizedProfileBlock(self.alphabet, [om1, om2])
        block.remove(om1)
        self.assertEqual(len(block), 1)
        self.assertIs(block[0], om2)

    def test_index(self):
        om1 = self._random_optimized_profile(b"profile1")
        om2 = self._random_optimized_profile(b"profile2")
        om3 = self._random_optimized_profile(b"profile3")

        block = OptimizedProfileBlock(self.alphabet, [om1, om2])
        self.assertEqual(len(block), 2)

        self.assertEqual(block.index(om1), 0)
        self.assertEqual(block.index(om2), 1)

        self.assertRaises(ValueError, block.index, om3)
        self.assertRaises(ValueError, block.index, om1, start=1)

    def test_pop(self):
        om1 = self._random_optimized_profile(b"profile1")
        om2 = self._random_optimized_profile(b"profile2")
        om3 = self._random_optimized_profile(b"profile3")

        block = OptimizedProfileBlock(self.alphabet, [om1, om2, om3])
        self.assertEqual(len(block), 3)

        self.assertIs(block.pop(), om3)
        self.assertIs(block.pop(0), om1)

        self.assertEqual(len(block), 1)
        self.assertIs(block[0], om2)

        self.assertIs(block.pop(-1), om2)
        self.assertEqual(len(block), 0)
        self.assertRaises(IndexError, block.pop)

        empty = OptimizedProfileBlock(self.alphabet)
        self.assertRaises(IndexError, block.pop)

    def test_setitem_slice(self):
        om1 = self._random_optimized_profile(b"profile1")
        om2 = self._random_optimized_profile(b"profile2")
        om3 = self._random_optimized_profile(b"profile3")

        block = OptimizedProfileBlock(self.alphabet, [om1, om2])
        self.assertEqual(len(block), 2)
        self.assertIs(block[0], om1)
        self.assertIs(block[1], om2)

        block[:] = [om1, om2, om3]
        self.assertEqual(len(block), 3)
        self.assertIs(block[0], om1)
        self.assertIs(block[1], om2)
        self.assertIs(block[2], om3)

        block[2:5] = [om3, om3, om3]
        self.assertEqual(len(block), 5)
        self.assertIs(block[3], om3)
        self.assertIs(block[4], om3)

    def test_delitem_slice(self):
        om1 = self._random_optimized_profile(b"profile1")
        om2 = self._random_optimized_profile(b"profile2")
        om3 = self._random_optimized_profile(b"profile3")

        block = OptimizedProfileBlock(self.alphabet, [om1, om2, om3])
        self.assertEqual(len(block), 3)
        del block[:]
        self.assertEqual(len(block), 0)

        block = OptimizedProfileBlock(self.alphabet, [om1, om2, om3])
        self.assertEqual(len(block), 3)
        del block[::2]
        self.assertEqual(len(block), 1)
        self.assertIs(block[0], om2)

    def test_contains(self):
        om1 = self._random_optimized_profile(b"profile1")
        om2 = self._random_optimized_profile(b"profile2")
        om3 = self._random_optimized_profile(b"profile3")

        block = OptimizedProfileBlock(self.alphabet, [om1, om2])
        self.assertEqual(len(block), 2)
        self.assertIn(om1, block)
        self.assertIn(om2, block)
        self.assertNotIn(om3, block)
        
        self.assertNotIn(42, block)
        self.assertNotIn(object(), block)

    def test_copy(self):
        om1 = self._random_optimized_profile(b"profile1")
        om2 = self._random_optimized_profile(b"profile2")
        om3 = self._random_optimized_profile(b"profile3")

        block = OptimizedProfileBlock(self.alphabet, [om1, om2, om3])
        block_copy = block.copy()

        self.assertEqual(len(block), len(block_copy))
        self.assertEqual(block, block_copy)

        self.assertIs(block[0], block_copy[0])
        self.assertIs(block[1], block_copy[1])
        self.assertIs(block[2], block_copy[2])

        block.remove(om3)
        self.assertEqual(len(block), 2)
        self.assertEqual(len(block_copy), 3)