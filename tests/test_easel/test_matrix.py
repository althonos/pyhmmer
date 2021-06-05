import unittest

from pyhmmer.easel import Matrix, MatrixF, MatrixU8, Vector, VectorF, VectorU8


class TestMatrix(unittest.TestCase):

    def test_abstract(self):
        self.assertRaises(TypeError, Matrix, [ [1, 2], [3, 4] ])
        self.assertRaises(TypeError, Matrix.zeros, 2, 2)


class TestMatrixF(unittest.TestCase):

    def test_init(self):
        mat = MatrixF([ [1, 2], [3, 4] ])
        self.assertEqual(mat[0, 0], 1)
        self.assertEqual(mat[0, 1], 2)
        self.assertEqual(mat[1, 0], 3)
        self.assertEqual(mat[1, 1], 4)

    def test_init_error(self):
        self.assertRaises(ValueError, MatrixF, [])
        self.assertRaises(ValueError, MatrixF, [ [], [] ])
        self.assertRaises(ValueError, MatrixF, [ [1.0, 2.0], [1.0] ])
        self.assertRaises(TypeError, MatrixF, 1)

    def test_len(self):
        mat = MatrixF([ [1, 2], [3, 4] ])
        self.assertEqual(len(mat), 2)
        mat = MatrixF.zeros(10, 10)
        self.assertEqual(len(mat), 10)

    def test_min(self):
        mat = MatrixF([ [1, 2], [3, 4] ])
        self.assertEqual(mat.min(), 1)

    def test_max(self):
        mat = MatrixF([ [1, 2], [3, 4] ])
        self.assertEqual(mat.max(), 4)

    def test_argmin(self):
        mat = MatrixF([ [1, 2], [3, 4] ])
        self.assertEqual(mat.argmin(), (0, 0))
        mat = MatrixF([ [2, 1], [4, 3] ])
        self.assertEqual(mat.argmin(), (0, 1))

    def test_argmax(self):
        mat = MatrixF([ [1, 2], [3, 4] ])
        self.assertEqual(mat.argmax(), (1, 1))
        mat = MatrixF([ [2, 1], [4, 3] ])
        self.assertEqual(mat.argmax(), (1, 0))

    def test_copy(self):
        mat = MatrixF([ [1, 2], [3, 4] ])
        mat2 = mat.copy()
        del mat
        self.assertEqual(mat2[0, 0], 1)
        self.assertEqual(mat2[0, 1], 2)
        self.assertEqual(mat2[1, 0], 3)
        self.assertEqual(mat2[1, 1], 4)

    def test_add(self):
        m1 = MatrixF([ [1, 2], [3, 4] ])
        m2 = MatrixF([ [2, 2], [3, 3] ])

        m3 = m1 + m2
        self.assertEqual(m3[0, 0], 3)
        self.assertEqual(m3[0, 1], 4)
        self.assertEqual(m3[1, 0], 6)
        self.assertEqual(m3[1, 1], 7)

        m4 = m2 + 2.0
        self.assertEqual(m4[0, 0], 4)
        self.assertEqual(m4[0, 1], 4)
        self.assertEqual(m4[1, 0], 5)
        self.assertEqual(m4[1, 1], 5)

    def test_iadd(self):
        mat = MatrixF([ [1, 2], [3, 4] ])
        mat += 3
        self.assertEqual(mat[0, 0], 4)
        self.assertEqual(mat[0, 1], 5)
        self.assertEqual(mat[1, 0], 6)
        self.assertEqual(mat[1, 1], 7)

        mat += MatrixF([ [2, 2], [3, 3] ])
        self.assertEqual(mat[0, 0], 6)
        self.assertEqual(mat[0, 1], 7)
        self.assertEqual(mat[1, 0], 9)
        self.assertEqual(mat[1, 1], 10)

    def test_mul(self):
        mat = MatrixF([ [1, 2], [3, 4] ])

        m2 = mat * 3
        self.assertEqual(m2[0, 0], 3)
        self.assertEqual(m2[0, 1], 6)
        self.assertEqual(m2[1, 0], 9)
        self.assertEqual(m2[1, 1], 12)

        m3 = m2 * MatrixF([ [2, 2], [3, 3] ])
        self.assertEqual(m3[0, 0], 6)
        self.assertEqual(m3[0, 1], 12)
        self.assertEqual(m3[1, 0], 27)
        self.assertEqual(m3[1, 1], 36)

    def test_imul(self):
        mat = MatrixF([ [1, 2], [3, 4] ])

        mat *= 3
        self.assertEqual(mat[0, 0], 3)
        self.assertEqual(mat[0, 1], 6)
        self.assertEqual(mat[1, 0], 9)
        self.assertEqual(mat[1, 1], 12)

        mat *= MatrixF([ [2, 2], [3, 3] ])
        self.assertEqual(mat[0, 0], 6)
        self.assertEqual(mat[0, 1], 12)
        self.assertEqual(mat[1, 0], 27)
        self.assertEqual(mat[1, 1], 36)

    def test_sum(self):
        mat = MatrixF([ [1, 2], [3, 4] ])
        self.assertEqual(mat.sum(), 1 + 2 + 3 + 4)

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


class TestMatrixU8(unittest.TestCase):

    def test_init(self):
        mat = MatrixU8([ [1, 2], [3, 4] ])
        self.assertEqual(mat[0, 0], 1)
        self.assertEqual(mat[0, 1], 2)
        self.assertEqual(mat[1, 0], 3)
        self.assertEqual(mat[1, 1], 4)

    def test_init_error(self):
        self.assertRaises(ValueError, MatrixU8, [])
        self.assertRaises(ValueError, MatrixU8, [ [], [] ])
        self.assertRaises(ValueError, MatrixU8, [ [1.0, 2.0], [1.0] ])
        self.assertRaises(TypeError, MatrixU8, 1)

    def test_len(self):
        mat = MatrixU8([ [1, 2], [3, 4] ])
        self.assertEqual(len(mat), 2)
        mat = MatrixU8.zeros(10, 10)
        self.assertEqual(len(mat), 10)

    def test_min(self):
        mat = MatrixU8([ [1, 2], [3, 4] ])
        self.assertEqual(mat.min(), 1)

    def test_max(self):
        mat = MatrixF([ [1, 2], [3, 4] ])
        self.assertEqual(mat.max(), 4)

    def test_argmin(self):
        mat = MatrixU8([ [1, 2], [3, 4] ])
        self.assertEqual(mat.argmin(), (0, 0))
        mat = MatrixU8([ [2, 1], [4, 3] ])
        self.assertEqual(mat.argmin(), (0, 1))

    def test_argmax(self):
        mat = MatrixU8([ [1, 2], [3, 4] ])
        self.assertEqual(mat.argmax(), (1, 1))
        mat = MatrixU8([ [2, 1], [4, 3] ])
        self.assertEqual(mat.argmax(), (1, 0))

    def test_copy(self):
        mat = MatrixU8([ [1, 2], [3, 4] ])
        mat2 = mat.copy()
        del mat
        self.assertEqual(mat2[0, 0], 1)
        self.assertEqual(mat2[0, 1], 2)
        self.assertEqual(mat2[1, 0], 3)
        self.assertEqual(mat2[1, 1], 4)

    def test_add(self):
        m1 = MatrixU8([ [1, 2], [3, 4] ])
        m2 = MatrixU8([ [2, 2], [3, 3] ])

        m3 = m1 + m2
        self.assertEqual(m3[0, 0], 3)
        self.assertEqual(m3[0, 1], 4)
        self.assertEqual(m3[1, 0], 6)
        self.assertEqual(m3[1, 1], 7)

        m4 = m2 + 2.0
        self.assertEqual(m4[0, 0], 4)
        self.assertEqual(m4[0, 1], 4)
        self.assertEqual(m4[1, 0], 5)
        self.assertEqual(m4[1, 1], 5)

    def test_iadd(self):
        mat = MatrixU8([ [1, 2], [3, 4] ])
        mat += 3
        self.assertEqual(mat[0, 0], 4)
        self.assertEqual(mat[0, 1], 5)
        self.assertEqual(mat[1, 0], 6)
        self.assertEqual(mat[1, 1], 7)

        mat += MatrixU8([ [2, 2], [3, 3] ])
        self.assertEqual(mat[0, 0], 6)
        self.assertEqual(mat[0, 1], 7)
        self.assertEqual(mat[1, 0], 9)
        self.assertEqual(mat[1, 1], 10)

    def test_mul(self):
        mat = MatrixU8([ [1, 2], [3, 4] ])

        m2 = mat * 3
        self.assertEqual(m2[0, 0], 3)
        self.assertEqual(m2[0, 1], 6)
        self.assertEqual(m2[1, 0], 9)
        self.assertEqual(m2[1, 1], 12)

        m3 = m2 * MatrixU8([ [2, 2], [3, 3] ])
        self.assertEqual(m3[0, 0], 6)
        self.assertEqual(m3[0, 1], 12)
        self.assertEqual(m3[1, 0], 27)
        self.assertEqual(m3[1, 1], 36)

    def test_imul(self):
        mat = MatrixU8([ [1, 2], [3, 4] ])

        mat *= 3
        self.assertEqual(mat[0, 0], 3)
        self.assertEqual(mat[0, 1], 6)
        self.assertEqual(mat[1, 0], 9)
        self.assertEqual(mat[1, 1], 12)

        mat *= MatrixU8([ [2, 2], [3, 3] ])
        self.assertEqual(mat[0, 0], 6)
        self.assertEqual(mat[0, 1], 12)
        self.assertEqual(mat[1, 0], 27)
        self.assertEqual(mat[1, 1], 36)

    def test_sum(self):
        mat = MatrixU8([ [1, 2], [3, 4] ])
        self.assertEqual(mat.sum(), 1 + 2 + 3 + 4)

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
