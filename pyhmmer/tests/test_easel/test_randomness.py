import copy
import pickle
import unittest
import sys

from pyhmmer.easel import Randomness


class TestRandomness(unittest.TestCase):

    @unittest.skipIf(sys.implementation.name == "pypy", "`getsizeof` not supported on PyPY")
    def test_sizeof(self):
        rng = Randomness(42, fast=True)
        self.assertGreater(sys.getsizeof(rng), 0)

    def test_init_fast(self):
        rng = Randomness(42, fast=True)
        self.assertTrue(rng.is_fast())

    def test_init_mersenne(self):
        rng = Randomness(42)
        self.assertFalse(rng.is_fast())

    def test_init_error(self):
        self.assertRaises(TypeError, Randomness, "ok")
        self.assertRaises(OverflowError, Randomness, -1)

    def test_copy(self):
        rng = Randomness(42)
        rng.random() # advance the rng
        new = copy.copy(rng)
        self.assertEqual(new.random(), rng.random())

    def test_pickle(self):
        rng = Randomness(42)
        rng.random() # advance the rng
        new = pickle.loads(pickle.dumps(rng))
        self.assertEqual(new.random(), rng.random())

    def test_normalvariate(self):
        rng1 = Randomness(42)
        n1 = rng1.normalvariate(0, 1)
        rng2 = Randomness(42)
        n2 = rng2.normalvariate(0, 1)
        self.assertEqual(n1, n2)

    def test_random(self):
        rng1 = Randomness(42)
        n1 = rng1.random()
        rng2 = Randomness(42)
        n2 = rng2.random()
        self.assertEqual(n1, n2)
