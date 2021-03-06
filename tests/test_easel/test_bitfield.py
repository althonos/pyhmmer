import unittest

from pyhmmer import easel


class TestBitfield(unittest.TestCase):
    def test_index_error(self):
        bitfield = easel.Bitfield(8)
        with self.assertRaises(IndexError):
            bitfield[9] = 1
        with self.assertRaises(IndexError):
            bitfield[-9] = 1
