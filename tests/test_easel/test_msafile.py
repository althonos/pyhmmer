import copy
import gc
import io
import os
import unittest
import tempfile

from pyhmmer import easel


class TestMSAFile(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.easel_folder = os.path.realpath(
            os.path.join(__file__, os.pardir, os.pardir, os.pardir, "vendor", "easel")
        )

    def test_author(self):
        trna_5 = os.path.join(self.easel_folder, "testsuite", "trna-5.stk")
        with easel.MSAFile(trna_5) as f:
            msa = f.read()
        self.assertEqual(msa.author, b"Infernal 0.1")

    def test_readformat_stockholm(self):
        stockholm = os.path.join(self.easel_folder, "formats", "stockholm.1")

        # check reading with specified format works
        with easel.MSAFile(stockholm, "stockholm") as f:
            msa = f.read()
            self.assertNotEqual(msa.name, b"")

        # check reading without specifying the format works too
        with easel.MSAFile(stockholm) as f:
            msa2 = f.read()
            self.assertNotEqual(msa2.name, b"")

        # check reading while giving another file format fails
        with easel.MSAFile(stockholm, "clustal") as f:
            self.assertRaises(ValueError, f.read)
