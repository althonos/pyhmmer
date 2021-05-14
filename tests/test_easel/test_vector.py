import unittest

from pyhmmer.easel import VectorF


class TestVectorF(unittest.TestCase):

    def test_init(self):
        vec = VectorF([1, 2, 3])
        self.assertEqual(vec[0], 1)
        self.assertEqual(vec[1], 2)
        self.assertEqual(vec[2], 3)

    def test_init_error(self):
        self.assertRaises(ValueError, VectorF, [])
        self.assertRaises(TypeError, VectorF, 1)

    def test_len(self):
        vec = VectorF([1, 2, 3])
        self.assertEqual(len(vec), 3)
        vec = VectorF.zeros(100)
        self.assertEqual(len(vec), 100)

    def test_min(self):
        vec = VectorF([1, 2, 3])
        self.assertEqual(vec.min(), 1)

    def test_max(self):
        vec = VectorF([1, 2, 3])
        self.assertEqual(vec.max(), 3)

    def test_argmin(self):
        vec = VectorF([4, 2, 8])
        self.assertEqual(vec.argmin(), 1)

    def test_argmax(self):
        vec = VectorF([2, 8, 4])
        self.assertEqual(vec.argmax(), 1)

    def test_copy(self):
        vec = VectorF([1, 2, 3])
        vec2 = vec.copy()
        del vec
        self.assertEqual(vec2[0], 1)
        self.assertEqual(vec2[1], 2)
        self.assertEqual(vec2[2], 3)

    def test_normalize(self):
        vec = VectorF([1, 3])
        vec.normalize()
        self.assertEqual(vec[0], 1/4)
        self.assertEqual(vec[1], 3/4)

    def test_reverse(self):
        vec = VectorF([1, 2, 3])
        vec.reverse()
        self.assertEqual(vec[0], 3)
        self.assertEqual(vec[1], 2)
        self.assertEqual(vec[2], 1)

    def test_add(self):
        vec = VectorF([1, 2, 3])
        vec2 = vec + 1
        self.assertEqual(vec2[0], 2)
        self.assertEqual(vec2[1], 3)
        self.assertEqual(vec2[2], 4)

    def test_iadd(self):
        vec = VectorF([1, 2, 3])
        vec += 3
        self.assertEqual(vec[0], 4)
        self.assertEqual(vec[1], 5)
        self.assertEqual(vec[2], 6)
        vec += VectorF([10, 11, 12])
        self.assertEqual(vec[0], 14)
        self.assertEqual(vec[1], 16)
        self.assertEqual(vec[2], 18)

    def test_imul(self):
        vec = VectorF([1, 2, 3])
        vec *= 3
        self.assertEqual(vec[0], 3)
        self.assertEqual(vec[1], 6)
        self.assertEqual(vec[2], 9)

    def test_matmul(self):
        u = VectorF([4, 5, 6])
        v = VectorF([1, 2, 3])
        self.assertEqual(u @ v, 1*4 + 2*5 + 3*6)

    def test_sum(self):
        vec = VectorF([1, 2, 3])
        self.assertEqual(vec.sum(), 1 + 2 + 3)

    def test_memoryview_tolist(self):
        vec = VectorF([1, 2, 3])
        mem = memoryview(vec)
        self.assertEqual(mem.tolist(), [1.0, 2.0, 3.0])
