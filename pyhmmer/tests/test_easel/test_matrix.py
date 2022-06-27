import pickle
import unittest
import sys

from pyhmmer.easel import Matrix, MatrixF, MatrixU8, Vector, VectorF, VectorU8


class _TestMatrixBase(object):

    Matrix = NotImplemented

    @unittest.skipIf(sys.implementation.name == "pypy", "`getsizeof` not supported on PyPY")
    def test_sizeof(self):
        m1 = self.Matrix([[1, 2], [3, 4]])
        m2 = self.Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        self.assertGreater(sys.getsizeof(m2), sys.getsizeof(m1))

    def test_pickle(self):
        m1 = self.Matrix([[1, 2], [3, 4]])
        m2 = pickle.loads(pickle.dumps(m1))
        self.assertEqual(m1.shape, m2.shape)
        self.assertSequenceEqual(memoryview(m1), memoryview(m2))
        # self.assertSequenceEqual(m1, m2)
        for i in range(m1.shape[0]):
            self.assertEqual(m1[i], m2[i])
            for j in range(m1.shape[1]):
                self.assertEqual(m1[i,j], m2[i,j])

    def test_pickle_protocol4(self):
        m1 = self.Matrix([[1, 2], [3, 4]])
        m2 = pickle.loads(pickle.dumps(m1, protocol=4))
        self.assertEqual(m1.shape, m2.shape)
        self.assertSequenceEqual(memoryview(m1), memoryview(m2))
        # self.assertSequenceEqual(m1, m2)
        for i in range(m1.shape[0]):
            self.assertEqual(m1[i], m2[i])
            for j in range(m1.shape[1]):
                self.assertEqual(m1[i,j], m2[i,j])

    @unittest.skipUnless(sys.version_info >= (3, 8), "pickle protocol 5 requires Python 3.8+")
    def test_pickle_protocol5(self):
        m1 = self.Matrix([[1, 2], [3, 4]])
        m2 = pickle.loads(pickle.dumps(m1, protocol=5))
        self.assertEqual(m1.shape, m2.shape)
        self.assertSequenceEqual(memoryview(m1), memoryview(m2))
        # self.assertSequenceEqual(m1, m2)
        for i in range(m1.shape[0]):
            self.assertEqual(m1[i], m2[i])
            for j in range(m1.shape[1]):
                self.assertEqual(m1[i,j], m2[i,j])

    def test_empty_matrix(self):
        m1 = self.Matrix([])
        m2 = self.Matrix.zeros(0, 0)
        self.assertEqual(m1, m2)
        self.assertFalse(m2)
        self.assertFalse(m1)

    def test_init(self):
        mat = self.Matrix([ [1, 2], [3, 4] ])
        self.assertEqual(mat[0, 0], 1)
        self.assertEqual(mat[0, 1], 2)
        self.assertEqual(mat[1, 0], 3)
        self.assertEqual(mat[1, 1], 4)

    def test_init_error(self):
        self.assertRaises(ValueError, self.Matrix, [ [1.0, 2.0], [1.0] ])
        self.assertRaises(TypeError, self.Matrix, 1)
        self.assertRaises(TypeError, self.Matrix.zeros, [])

    def test_zeros(self):
        m = self.Matrix.zeros(2, 2)
        self.assertEqual(m.shape, (2, 2))
        self.assertEqual(m, self.Matrix([ [0, 0], [0, 0] ]))

        m2 = self.Matrix.zeros(4, 3)
        self.assertEqual(m2.shape, (4, 3))

        m3 = self.Matrix.zeros(0, 0)
        self.assertEqual(m3.shape, (0, 0))

        m4 = self.Matrix.zeros(0, 2)
        self.assertEqual(m4.shape, (0, 2))

        self.assertRaises(ValueError, self.Matrix.zeros, -2, 2)
        self.assertRaises(ValueError, self.Matrix.zeros, 2, -2)
        self.assertRaises(ValueError, self.Matrix.zeros, 0, -2)

    def test_eq(self):
        m1 = self.Matrix([ [1,2], [3,4] ])
        m2 = self.Matrix([ [1,2], [3,4] ])
        self.assertEqual(m1, m2)
        m3 = self.Matrix([ [2,2], [3,4] ])
        self.assertNotEqual(m1, m3)
        m4 = self.Matrix([ [1,2], [3,4], [5,6] ])
        self.assertNotEqual(m1, m4)

    def test_len(self):
        mat = self.Matrix([[1, 2], [3, 4]])
        self.assertEqual(len(mat), 2)
        mat = self.Matrix.zeros(10, 10)
        self.assertEqual(len(mat), 10)

    def test_min(self):
        mat = self.Matrix([ [1, 2], [3, 4] ])
        self.assertEqual(mat.min(), 1)

    def test_max(self):
        mat = self.Matrix([ [1, 2], [3, 4] ])
        self.assertEqual(mat.max(), 4)

    def test_argmin(self):
        mat = self.Matrix([ [1, 2], [3, 4] ])
        self.assertEqual(mat.argmin(), (0, 0))
        mat = self.Matrix([ [2, 1], [4, 3] ])
        self.assertEqual(mat.argmin(), (0, 1))

    def test_argmax(self):
        mat = self.Matrix([ [1, 2], [3, 4] ])
        self.assertEqual(mat.argmax(), (1, 1))
        mat = self.Matrix([ [2, 1], [4, 3] ])
        self.assertEqual(mat.argmax(), (1, 0))

    def test_copy(self):
        mat = self.Matrix([ [1, 2], [3, 4] ])
        mat2 = mat.copy()
        del mat
        self.assertIsInstance(mat2, self.Matrix)
        self.assertEqual(mat2[0, 0], 1)
        self.assertEqual(mat2[0, 1], 2)
        self.assertEqual(mat2[1, 0], 3)
        self.assertEqual(mat2[1, 1], 4)

    def test_add_matrix(self):
        m1 = self.Matrix([ [1, 2], [3, 4] ])
        m2 = self.Matrix([ [2, 2], [3, 3] ])
        m3 = m1 + m2
        self.assertEqual(m3[0, 0], 3)
        self.assertEqual(m3[0, 1], 4)
        self.assertEqual(m3[1, 0], 6)
        self.assertEqual(m3[1, 1], 7)

    def test_add_scalar(self):
        m1 = self.Matrix([ [2, 2], [3, 3] ])
        m2 = m1 + 2.0
        self.assertEqual(m2[0, 0], 4)
        self.assertEqual(m2[0, 1], 4)
        self.assertEqual(m2[1, 0], 5)
        self.assertEqual(m2[1, 1], 5)

    def test_iadd_scalar(self):
        mat = self.Matrix([ [1, 2], [3, 4] ])
        mat += 3
        self.assertEqual(mat[0, 0], 4)
        self.assertEqual(mat[0, 1], 5)
        self.assertEqual(mat[1, 0], 6)
        self.assertEqual(mat[1, 1], 7)

    def test_iadd_matrix(self):
        mat = self.Matrix([ [4, 5], [6, 7] ])
        mat += self.Matrix([ [2, 2], [3, 3] ])
        self.assertEqual(mat[0, 0], 6)
        self.assertEqual(mat[0, 1], 7)
        self.assertEqual(mat[1, 0], 9)
        self.assertEqual(mat[1, 1], 10)

    def test_mul_scalar(self):
        mat = self.Matrix([ [1, 2], [3, 4] ])
        m2 = mat * 3
        self.assertEqual(m2[0, 0], 3)
        self.assertEqual(m2[0, 1], 6)
        self.assertEqual(m2[1, 0], 9)
        self.assertEqual(m2[1, 1], 12)

    def test_mul_matrix(self):
        mat = self.Matrix([ [3, 6], [9, 12] ])
        m2 = mat * self.Matrix([ [2, 2], [3, 3] ])
        self.assertEqual(m2[0, 0], 6)
        self.assertEqual(m2[0, 1], 12)
        self.assertEqual(m2[1, 0], 27)
        self.assertEqual(m2[1, 1], 36)

    def test_imul_scalar(self):
        mat = self.Matrix([ [1, 2], [3, 4] ])
        mat *= 3
        self.assertEqual(mat[0, 0], 3)
        self.assertEqual(mat[0, 1], 6)
        self.assertEqual(mat[1, 0], 9)
        self.assertEqual(mat[1, 1], 12)

    def test_imul_matrix(self):
        mat = self.Matrix([ [3, 6], [9, 12] ])
        mat *= self.Matrix([ [2, 2], [3, 3] ])
        self.assertEqual(mat[0, 0], 6)
        self.assertEqual(mat[0, 1], 12)
        self.assertEqual(mat[1, 0], 27)
        self.assertEqual(mat[1, 1], 36)

    def test_sum(self):
        mat = self.Matrix([ [1, 2], [3, 4] ])
        self.assertEqual(mat.sum(), 1 + 2 + 3 + 4)

    def test_getitem_slice(self):
        mat = self.Matrix([ [1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12], [13, 14, 15, 16] ])

        m2 = mat[:2]
        self.assertEqual(m2, self.Matrix([[1, 2, 3, 4], [5, 6, 7, 8]]))

        m3 = mat[1:3]
        self.assertEqual(m3, self.Matrix([[5, 6, 7, 8], [9, 10, 11, 12]]))

        m4 = mat[-3:-1]
        self.assertEqual(m4, self.Matrix([[5, 6, 7, 8], [9, 10, 11, 12]]))


class TestMatrix(unittest.TestCase):

    def test_abstract(self):
        self.assertRaises(TypeError, Matrix, [ [1, 2], [3, 4] ])
        self.assertRaises(TypeError, Matrix.zeros, 2, 2)


class TestMatrixF(_TestMatrixBase, unittest.TestCase):

    Matrix = MatrixF

    def test_memoryview_tolist(self):
        mat = MatrixF([ [1, 2], [3, 4] ])
        mem = memoryview(mat)
        self.assertEqual(mem.tolist(), [ [1.0, 2.0], [3.0, 4.0] ])

    def test_getitem_row(self):
        mat = MatrixF([ [1, 2], [3, 4] ])

        self.assertIsInstance(mat[0], VectorF)
        self.assertEqual(list(mat[0]), [1.0, 2.0])
        self.assertEqual(list(mat[-1]), [3.0, 4.0])

        with self.assertRaises(IndexError):
            row = mat[-10]

    def test_getitem_element(self):
        mat = MatrixF([ [1, 2], [3, 4] ])

        self.assertEqual(mat[0, 0], 1.0)
        self.assertEqual(mat[0, 1], 2.0)
        self.assertEqual(mat[1, 0], 3.0)
        self.assertEqual(mat[1, 1], 4.0)

        self.assertEqual(mat[0, -1], 2.0)
        self.assertEqual(mat[-1, 0], 3.0)
        self.assertEqual(mat[-1, -1], 4.0)

        with self.assertRaises(IndexError):
            x = mat[-10, 0]
        with self.assertRaises(IndexError):
            x = mat[0, -10]
        with self.assertRaises(IndexError):
            x = mat[0, 10]
        with self.assertRaises(IndexError):
            x = mat[10, 0]

    def test_from_raw_bytes_littleendian(self):
        mat = self.Matrix._from_raw_bytes(b'\x00\x00\x80>', 1, 1, byteorder="little")
        self.assertEqual(len(mat), 1)
        self.assertEqual(mat[0][0], 0.25)

    def test_from_raw_bytes_bigendian(self):
        mat = self.Matrix._from_raw_bytes(b'>\x80\x00\x00', 1, 1, byteorder="big")
        self.assertEqual(len(mat), 1)
        self.assertEqual(mat[0][0], 0.25)


class TestMatrixU8(_TestMatrixBase, unittest.TestCase):

    Matrix = MatrixU8

    def test_memoryview_tolist(self):
        mat = MatrixU8([ [1, 2], [3, 4] ])
        mem = memoryview(mat)
        self.assertEqual(mem.tolist(), [ [1, 2], [3, 4] ])

    def test_getitem_row(self):
        mat = MatrixU8([ [1, 2], [3, 4] ])

        self.assertIsInstance(mat[0], VectorU8)
        self.assertEqual(list(mat[0]), [1, 2])
        self.assertEqual(list(mat[-1]), [3, 4])

        with self.assertRaises(IndexError):
            row = mat[-10]

    def test_getitem_element(self):
        mat = MatrixU8([ [1, 2], [3, 4] ])

        self.assertEqual(mat[0, 0], 1)
        self.assertEqual(mat[0, 1], 2)
        self.assertEqual(mat[1, 0], 3)
        self.assertEqual(mat[1, 1], 4)

        self.assertEqual(mat[0, -1], 2)
        self.assertEqual(mat[-1, 0], 3)
        self.assertEqual(mat[-1, -1], 4)

        with self.assertRaises(IndexError):
            x = mat[-10, 0]
        with self.assertRaises(IndexError):
            x = mat[0, -10]
        with self.assertRaises(IndexError):
            x = mat[0, 10]
        with self.assertRaises(IndexError):
            x = mat[10, 0]
