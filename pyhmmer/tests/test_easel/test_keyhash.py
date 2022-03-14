import copy
import pickle
import unittest
import sys

from pyhmmer import easel


class TestKeyHash(unittest.TestCase):

    @unittest.skipIf(sys.implementation.name == "pypy", "`getsizeof` not supported on PyPY")
    def test_getsizeof(self):
        kh = easel.KeyHash()
        self.assertGreater(sys.getsizeof(kh), 0)

    def test_len(self):
        kh = easel.KeyHash()
        self.assertEqual(len(kh), 0)
        kh.add(b"first")
        self.assertEqual(len(kh), 1)
        kh.add(b"second")
        self.assertEqual(len(kh), 2)
        kh.add(b"first")
        self.assertEqual(len(kh), 2)
        kh.clear()
        self.assertEqual(len(kh), 0)

    def test_iter(self):
        kh = easel.KeyHash()
        self.assertRaises(StopIteration, next, iter(kh))
        kh.add(b"first")
        kh.add(b"second")
        self.assertEqual(list(iter(kh)), [b"first", b"second"])

    def test_contains(self):
        kh = easel.KeyHash()
        kh.add(b"first")
        self.assertIn(b"first", kh)
        self.assertNotIn(b"second", kh)
        self.assertNotIn(1, kh)

    def test_clear(self):
        kh = easel.KeyHash()
        kh.add(b"first")
        self.assertIn(b"first", kh)
        self.assertEqual(len(kh), 1)
        kh.clear()
        self.assertNotIn(b"first", kh)
        self.assertEqual(len(kh), 0)

    def test_copy(self):
        kh = easel.KeyHash()
        kh.add(b"first")
        k2 = copy.copy(kh)
        self.assertEqual(list(kh), list(k2))
        kh.clear()
        self.assertNotEqual(list(kh), list(k2))

    def test_pickle(self):
        kh = easel.KeyHash()
        kh.add(b"first")
        kh.add(b"second")

        kh2 = pickle.loads(pickle.dumps(kh))
        self.assertEqual(len(kh2), 2)
        self.assertEqual(list(iter(kh2)), [b"first", b"second"])
        self.assertEqual(kh.__getstate__(), kh2.__getstate__())
