import array
import unittest
import pickle
import struct

from pyhmmer.easel import Vector, VectorF, VectorU8


class _TestVectorBase(object):

    Vector = NotImplemented

    def test_pickle(self):
        v1 = self.Vector([1, 2, 3, 4, 5, 6])
        v2 = pickle.loads(pickle.dumps(v1))

        for i in range(len(v1)):
            self.assertEqual(v1[i], v2[i])

    def test_empty_vector(self):
        v1 = self.Vector([])
        v2 = self.Vector.zeros(0)
        self.assertEqual(len(v1), 0)
        self.assertEqual(len(v2), 0)
        self.assertEqual(v1, v2)
        self.assertFalse(v2)
        self.assertFalse(v1)

    def test_init(self):
        vec = self.Vector([1, 2, 3])
        self.assertEqual(vec[0], 1)
        self.assertEqual(vec[1], 2)
        self.assertEqual(vec[2], 3)

    def test_init_error(self):
        self.assertRaises(TypeError, self.Vector, 1)
        self.assertRaises(TypeError, self.Vector.zeros, [1, 2, 3])
        self.assertRaises(TypeError, self.Vector.zeros, "1")

    def test_shape(self):
        vec = self.Vector([1, 2, 3])
        self.assertEqual(vec.shape, (3,))

    def test_len(self):
        vec = self.Vector([1, 2, 3])
        self.assertEqual(len(vec), 3)
        vec = self.Vector.zeros(100)
        self.assertEqual(len(vec), 100)

    def test_copy(self):
        vec = self.Vector([1, 2, 3])
        vec2 = vec.copy()
        del vec
        self.assertIsInstance(vec2, self.Vector)
        self.assertEqual(vec2[0], 1)
        self.assertEqual(vec2[1], 2)
        self.assertEqual(vec2[2], 3)

    def test_reverse(self):
        vec = self.Vector([1, 2, 3])
        vec.reverse()
        self.assertEqual(vec[0], 3)
        self.assertEqual(vec[1], 2)
        self.assertEqual(vec[2], 1)
        vec2 = self.Vector([1, 2, 3, 4])
        vec2.reverse()
        self.assertEqual(vec2[0], 4)
        self.assertEqual(vec2[1], 3)
        self.assertEqual(vec2[2], 2)
        self.assertEqual(vec2[3], 1)

    def test_add(self):
        vec = self.Vector([1, 2, 3])
        vec2 = vec + 1
        self.assertEqual(vec2[0], 2)
        self.assertEqual(vec2[1], 3)
        self.assertEqual(vec2[2], 4)

        with self.assertRaises(ValueError):
            vec + self.Vector([1])

    def test_iadd_scalar(self):
        vec = self.Vector([1, 2, 3])
        vec += 3
        self.assertEqual(vec[0], 4)
        self.assertEqual(vec[1], 5)
        self.assertEqual(vec[2], 6)

    def test_iadd_vector(self):
        vec = self.Vector([4, 5, 6])
        vec += self.Vector([10, 11, 12])
        self.assertEqual(vec[0], 14)
        self.assertEqual(vec[1], 16)
        self.assertEqual(vec[2], 18)

    def test_mul_scalar(self):
        vec = self.Vector([1, 2, 3])
        v2 = vec * 3
        self.assertEqual(v2[0], 3)
        self.assertEqual(v2[1], 6)
        self.assertEqual(v2[2], 9)

    def test_mul_vector(self):
        vec = self.Vector([1, 2, 3])
        v2 = self.Vector([3, 6, 9])
        v3 = vec * v2
        self.assertEqual(v3[0], 3)
        self.assertEqual(v3[1], 12)
        self.assertEqual(v3[2], 27)

    def test_imul_scalar(self):
        vec = self.Vector([1, 2, 3])
        vec *= 3
        self.assertEqual(vec[0], 3)
        self.assertEqual(vec[1], 6)
        self.assertEqual(vec[2], 9)

    def test_matmul_vector(self):
        u = self.Vector([4, 5, 6])
        v = self.Vector([1, 2, 3])
        self.assertEqual(u @ v, 1*4 + 2*5 + 3*6)

    def test_sum(self):
        vec = self.Vector([1, 2, 3])
        self.assertEqual(vec.sum(), 1 + 2 + 3)

    def test_slice(self):
        vec = self.Vector([1, 2, 3, 4])

        v1 = vec[:]
        self.assertEqual(len(v1), 4)
        self.assertEqual(v1[0], 1)
        self.assertEqual(v1[-1], 4)

        v2 = vec[1:3]
        self.assertEqual(len(v2), 2)
        self.assertEqual(v2[0], 2)
        self.assertEqual(v2[1], 3)

        v3 = vec[:-1]
        self.assertEqual(len(v3), 3)
        self.assertEqual(v3[-1], 3)

        v4 = vec[0:10]
        self.assertEqual(len(v4), 4)
        self.assertEqual(v4[-1], 4)

        with self.assertRaises(ValueError):
            vec[::-1]

    def test_min(self):
        vec = self.Vector([1, 2, 3])
        self.assertEqual(vec.min(), 1)

    def test_max(self):
        vec = self.Vector([1, 2, 3])
        self.assertEqual(vec.max(), 3)

    def test_argmin(self):
        vec = self.Vector([4, 2, 8])
        self.assertEqual(vec.argmin(), 1)

    def test_argmax(self):
        vec = self.Vector([2, 8, 4])
        self.assertEqual(vec.argmax(), 1)


class TestVector(unittest.TestCase):

    def test_abstract(self):
        self.assertRaises(TypeError, Vector, [1, 2, 3])
        self.assertRaises(TypeError, Vector.zeros, 1)


class TestVectorF(_TestVectorBase, unittest.TestCase):

    Vector = VectorF

    def test_strides(self):
        vec = self.Vector([1, 2, 3])
        sizeof_float = len(struct.pack('f', 1.0))
        self.assertEqual(vec.strides, (sizeof_float,))

    def test_normalize(self):
        vec = self.Vector([1, 3])
        vec.normalize()
        self.assertEqual(vec[0], 1/4)
        self.assertEqual(vec[1], 3/4)

    def test_memoryview_tolist(self):
        vec = self.Vector([1, 2, 3])
        mem = memoryview(vec)
        self.assertEqual(mem.tolist(), [1.0, 2.0, 3.0])


class TestVectorU8(_TestVectorBase, unittest.TestCase):

    Vector = VectorU8

    def test_strides(self):
        vec = self.Vector([1, 2, 3])
        sizeof_u8 = len(struct.pack('B', 1))
        self.assertEqual(vec.strides, (sizeof_u8,))

    def test_sum_wrapping(self):
        vec = self.Vector([124, 72, 116])
        self.assertEqual(vec.sum(), (124 + 72 + 116) % 256)

    def test_memoryview_tolist(self):
        vec = self.Vector([1, 2, 3])
        mem = memoryview(vec)
        self.assertEqual(mem.tolist(), [1, 2, 3])

    def test_eq_bytebuffer(self):
        vec = self.Vector([1, 2, 3])
        b1 = bytearray([1, 2, 3])
        self.assertEqual(vec, b1)

        b2 = array.array('B', [1, 2, 3])
        self.assertEqual(vec, b2)

        b3 = array.array('B', [1, 2, 3, 4])
        self.assertNotEqual(vec, b3)

        b4 = array.array('L', [1, 2, 3])
        self.assertNotEqual(vec, b4)
