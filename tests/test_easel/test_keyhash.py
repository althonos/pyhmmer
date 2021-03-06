import copy
import unittest

from pyhmmer import easel


class TestKeyHash(unittest.TestCase):

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
