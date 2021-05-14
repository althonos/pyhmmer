import unittest

from pyhmmer.easel import MatrixF


class TestVectorF(unittest.TestCase):

    def test_init(self):
        mat = MatrixF([ [1, 2], [3, 4] ])
        self.assertEqual(mat[0, 0], 1)
        self.assertEqual(mat[0, 1], 2)
        self.assertEqual(mat[1, 0], 3)
        self.assertEqual(mat[1, 1], 4)

    def test_init_error(self):
        self.assertRaises(ValueError, MatrixF, [])
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

    def test_iadd(self):
        mat = MatrixF([ [1, 2], [3, 4] ])
        mat += 3
        self.assertEqual(mat[0, 0], 4)
        self.assertEqual(mat[0, 1], 5)
        self.assertEqual(mat[1, 0], 6)
        self.assertEqual(mat[1, 1], 7)

    def test_imul(self):
        mat = MatrixF([ [1, 2], [3, 4] ])
        mat *= 3
        self.assertEqual(mat[0, 0], 3)
        self.assertEqual(mat[0, 1], 6)
        self.assertEqual(mat[1, 0], 9)
        self.assertEqual(mat[1, 1], 12)

    def test_sum(self):
        mat = MatrixF([ [1, 2], [3, 4] ])
        self.assertEqual(mat.sum(), 1 + 2 + 3 + 4)

    def test_memoryview_tolist(self):
        mat = MatrixF([ [1, 2], [3, 4] ])
        mem = memoryview(mat)
        self.assertEqual(mem.tolist(), [ [1.0, 2.0], [3.0, 4.0] ])
