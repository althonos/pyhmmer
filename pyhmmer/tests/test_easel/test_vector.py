import array
import unittest
import pickle
import struct
import sys

from pyhmmer.easel import Vector, VectorF, VectorU8


class _TestVectorBase(object):

    Vector = NotImplemented

    def test_pickle(self):
        v1 = self.Vector(range(6))
        v2 = pickle.loads(pickle.dumps(v1))
        self.assertSequenceEqual(v1, v2)

    def test_pickle_protocol4(self):
        v1 = self.Vector(range(6))
        v2 = pickle.loads(pickle.dumps(v1, protocol=4))
        self.assertEqual(v1.shape, v2.shape)
        self.assertSequenceEqual(v1, v2)
        self.assertSequenceEqual(memoryview(v1), memoryview(v2))

    @unittest.skipUnless(sys.version_info >= (3, 8), "pickle protocol 5 requires Python 3.8+")
    def test_pickle_protocol5(self):
        v1 = self.Vector(range(6))
        v2 = pickle.loads(pickle.dumps(v1, protocol=5))
        self.assertEqual(v1.shape, v2.shape)
        self.assertSequenceEqual(v1, v2)
        self.assertSequenceEqual(memoryview(v1), memoryview(v2))

    def test_empty_vector(self):
        v1 = self.Vector([])
        v2 = self.Vector.zeros(0)
        v3 = self.Vector()
        self.assertEqual(len(v1), 0)
        self.assertEqual(len(v2), 0)
        self.assertEqual(len(v3), 0)
        self.assertEqual(v1, v2)
        self.assertEqual(v1, v3)
        self.assertFalse(v1)
        self.assertFalse(v2)
        self.assertFalse(v3)

        if sys.implementation.name != "pypy":
            v3 = self.Vector.zeros(3)
            self.assertLess(sys.getsizeof(v1), sys.getsizeof(v3))

    def test_init(self):
        vec = self.Vector([1, 2, 3])
        self.assertEqual(vec[0], 1)
        self.assertEqual(vec[1], 2)
        self.assertEqual(vec[2], 3)

    def test_init_memcpy(self):
        v1 = self.Vector([1, 2, 3])
        a  = array.array(v1.format, v1)
        v2 = self.Vector(a)
        self.assertEqual(v1, v2)

    def test_init_error(self):
        self.assertRaises(TypeError, self.Vector, 1)
        self.assertRaises(TypeError, self.Vector.zeros, [1, 2, 3])
        self.assertRaises(TypeError, self.Vector.zeros, "1")

    def test_shape(self):
        vec = self.Vector([1, 2, 3])
        self.assertEqual(vec.shape, (3,))
        vec2 = self.Vector.zeros(100)
        self.assertEqual(vec2.shape, (100,))
        vec3 = self.Vector.zeros(0)
        self.assertEqual(vec3.shape, (0,))

    def test_len(self):
        vec = self.Vector([1, 2, 3])
        self.assertEqual(len(vec), 3)
        vec2 = self.Vector.zeros(100)
        self.assertEqual(len(vec2), 100)
        vec3 = self.Vector([])
        self.assertEqual(len(vec3), 0)

    def test_copy(self):
        vec = self.Vector([1, 2, 3])
        vec2 = vec.copy()
        del vec
        self.assertIsInstance(vec2, self.Vector)
        self.assertEqual(vec2[0], 1)
        self.assertEqual(vec2[1], 2)
        self.assertEqual(vec2[2], 3)

        vec3 = self.Vector([])
        vec4 = vec3.copy()
        self.assertEqual(vec3, vec4)
        self.assertEqual(len(vec4), 0)

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

        vec3 = self.Vector([])
        vec3.reverse()
        self.assertEqual(vec3, self.Vector([]))
        self.assertEqual(len(vec3), 0)

    def test_add(self):
        vec = self.Vector([1, 2, 3])
        vec2 = vec + 1
        self.assertEqual(vec2[0], 2)
        self.assertEqual(vec2[1], 3)
        self.assertEqual(vec2[2], 4)

        with self.assertRaises(ValueError):
            vec + self.Vector([1])

        v2 = self.Vector([])
        v3 = v2 + self.Vector([])
        self.assertEqual(v3, self.Vector([]))

    def test_iadd_scalar(self):
        vec = self.Vector([1, 2, 3])
        vec += 3
        self.assertEqual(vec[0], 4)
        self.assertEqual(vec[1], 5)
        self.assertEqual(vec[2], 6)

        v2 = self.Vector([])
        v2 += 3
        self.assertEqual(v2, self.Vector([]))

    def test_iadd_vector(self):
        vec = self.Vector([4, 5, 6])
        vec += self.Vector([10, 11, 12])
        self.assertEqual(vec[0], 14)
        self.assertEqual(vec[1], 16)
        self.assertEqual(vec[2], 18)

        with self.assertRaises(ValueError):
            vec += self.Vector([1])

        v2 = self.Vector([])
        v2 += self.Vector([])
        self.assertEqual(v2, self.Vector([]))

    def test_sub(self):
        vec = self.Vector([1, 2, 3])
        v2 = vec - 1
        self.assertEqual(v2[0], 0)
        self.assertEqual(v2[1], 1)
        self.assertEqual(v2[2], 2)

        v3 = self.Vector([8, 10, 12])
        v4 = self.Vector([1, 2, 3])
        v5 = v3 - v4
        self.assertEqual(v5[0], 7)
        self.assertEqual(v5[1], 8)
        self.assertEqual(v5[2], 9)

    def test_isub_scalar(self):
        vec = self.Vector([4, 5, 6])
        vec -= 2
        self.assertEqual(vec[0], 2)
        self.assertEqual(vec[1], 3)
        self.assertEqual(vec[2], 4)

    def test_isub_vector(self):
        vec = self.Vector([4, 5, 6])
        vec -= self.Vector([2, 3, 2])
        self.assertEqual(vec[0], 2)
        self.assertEqual(vec[1], 2)
        self.assertEqual(vec[2], 4)

    def test_mul_scalar(self):
        vec = self.Vector([1, 2, 3])
        v2 = vec * 3
        self.assertEqual(v2[0], 3)
        self.assertEqual(v2[1], 6)
        self.assertEqual(v2[2], 9)

        v2 = self.Vector([])
        v3 = v2 * 3
        self.assertEqual(v3, self.Vector([]))

    def test_mul_vector(self):
        vec = self.Vector([1, 2, 3])
        v2 = self.Vector([3, 6, 9])
        v3 = vec * v2
        self.assertEqual(v3[0], 3)
        self.assertEqual(v3[1], 12)
        self.assertEqual(v3[2], 27)

        v2 = self.Vector([])
        v3 = v2 * self.Vector([])
        self.assertEqual(v3, self.Vector([]))

    def test_imul_scalar(self):
        vec = self.Vector([1, 2, 3])
        vec *= 3
        self.assertEqual(vec[0], 3)
        self.assertEqual(vec[1], 6)
        self.assertEqual(vec[2], 9)

        v2 = self.Vector([])
        v2 *= 3
        self.assertEqual(v2, self.Vector([]))

    def test_matmul_vector(self):
        u = self.Vector([4, 5, 6])
        v = self.Vector([1, 2, 3])
        self.assertEqual(u @ v, 1*4 + 2*5 + 3*6)

        x = self.Vector([])
        y = self.Vector([])
        self.assertEqual(x @ y, 0)

    def test_sum(self):
        vec = self.Vector([1, 2, 3])
        self.assertEqual(vec.sum(), 1 + 2 + 3)

        vec2 = self.Vector([])
        self.assertEqual(vec2.sum(), 0)

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

        v2 = self.Vector([])
        self.assertRaises(ValueError, v2.min)

    def test_max(self):
        vec = self.Vector([1, 2, 3])
        self.assertEqual(vec.max(), 3)

        v2 = self.Vector([])
        self.assertRaises(ValueError, v2.max)

    def test_argmin(self):
        vec = self.Vector([4, 2, 8])
        self.assertEqual(vec.argmin(), 1)

        v2 = self.Vector([])
        self.assertRaises(ValueError, v2.argmin)

    def test_argmax(self):
        vec = self.Vector([2, 8, 4])
        self.assertEqual(vec.argmax(), 1)

        v2 = self.Vector([])
        self.assertRaises(ValueError, v2.argmax)


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

        vec2 = self.Vector([])
        vec2.normalize()

    def test_memoryview_tolist(self):
        vec = self.Vector([1, 2, 3])
        mem = memoryview(vec)
        self.assertEqual(mem.tolist(), [1.0, 2.0, 3.0])

    def test_neg(self):
        vec = self.Vector([1, 2, 3])
        v2 = -vec
        self.assertEqual(v2[0], -1)
        self.assertEqual(v2[1], -2)
        self.assertEqual(v2[2], -3)

    def test_div_scalar(self):
        vec = self.Vector([1, 2, 3])
        v2 = vec / 2
        self.assertEqual(v2[0], 0.5)
        self.assertEqual(v2[1], 1.0)
        self.assertEqual(v2[2], 1.5)

        v2 = self.Vector([])
        v3 = v2 / 3
        self.assertEqual(v3, self.Vector([]))

    def test_div_vector(self):
        vec = self.Vector([1, 2, 3])
        v2 = self.Vector([2, 4, 6])
        v3 = vec / v2
        self.assertEqual(v3[0], 0.5)
        self.assertEqual(v3[1], 0.5)
        self.assertEqual(v3[2], 0.5)

        v2 = self.Vector([])
        v3 = v2 / self.Vector([])
        self.assertEqual(v3, self.Vector([]))

    def test_idiv_scalar(self):
        vec = self.Vector([1, 2, 3])
        vec /= 2
        self.assertEqual(vec[0], 0.5)
        self.assertEqual(vec[1], 1.0)
        self.assertEqual(vec[2], 1.5)

        vec = self.Vector([])
        vec /= 3
        self.assertEqual(vec, self.Vector([]))

    def test_idiv_vector(self):
        vec = self.Vector([1, 2, 3])
        vec /= self.Vector([2, 4, 6])
        self.assertEqual(vec[0], 0.5)
        self.assertEqual(vec[1], 0.5)
        self.assertEqual(vec[2], 0.5)

        vec = self.Vector([])
        vec /= self.Vector([])
        self.assertEqual(vec, self.Vector([]))

    def test_from_raw_bytes_littleendian(self):
        vec = self.Vector._from_raw_bytes(b'\x00\x00\x80>', 1, byteorder="little")
        self.assertEqual(len(vec), 1)
        self.assertEqual(vec[0], 0.25)

    def test_from_raw_bytes_bigendian(self):
        vec = self.Vector._from_raw_bytes(b'>\x80\x00\x00', 1, byteorder="big")
        self.assertEqual(len(vec), 1)
        self.assertEqual(vec[0], 0.25)


class TestVectorU8(_TestVectorBase, unittest.TestCase):

    Vector = VectorU8

    def test_strides(self):
        vec = self.Vector([1, 2, 3])
        sizeof_u8 = len(struct.pack('B', 1))
        self.assertEqual(vec.strides, (sizeof_u8,))

    def test_isub_wrapping(self):
        vec = self.Vector([0, 1, 2])
        vec -= 1
        self.assertEqual(vec[0], 255)
        self.assertEqual(vec[1], 0)
        self.assertEqual(vec[2], 1)

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

    def test_floordiv_scalar(self):
        vec = self.Vector([1, 2, 3])
        v2 = vec // 2
        self.assertEqual(v2[0], 0)
        self.assertEqual(v2[1], 1)
        self.assertEqual(v2[2], 1)

        v2 = self.Vector([])
        v3 = v2 // 3
        self.assertEqual(v3, self.Vector([]))

    def test_floordiv_vector(self):
        vec = self.Vector([1, 2, 3])
        v2 = self.Vector([2, 4, 1])
        v3 = vec // v2
        self.assertEqual(v3[0], 0)
        self.assertEqual(v3[1], 0)
        self.assertEqual(v3[2], 3)

        v2 = self.Vector([])
        v3 = v2 // self.Vector([])
        self.assertEqual(v3, self.Vector([]))

    def test_ifloordiv_scalar(self):
        vec = self.Vector([1, 2, 3])
        vec //= 2
        self.assertEqual(vec[0], 0)
        self.assertEqual(vec[1], 1)
        self.assertEqual(vec[2], 1)

        vec = self.Vector([])
        vec //= 3
        self.assertEqual(vec, self.Vector([]))

    def test_ifloordiv_vector(self):
        vec = self.Vector([1, 2, 3])
        vec //= self.Vector([2, 4, 6])
        self.assertEqual(vec[0], 0)
        self.assertEqual(vec[1], 0)
        self.assertEqual(vec[2], 0)

        vec = self.Vector([])
        vec //= self.Vector([])
        self.assertEqual(vec, self.Vector([]))
