import copy
import gc
import io
import os
import unittest
import tempfile

from pyhmmer import easel


class TestMSA(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.formats_folder = os.path.realpath(
            os.path.join(
                __file__, os.pardir, os.pardir, os.pardir, "vendor", "easel", "formats"
            )
        )

    def test_write_roundtrip_stockholm(self):
        sto = os.path.join(self.formats_folder, "stockholm.1")
        with easel.MSAFile(sto, "stockholm") as msa_file:
            msa = msa_file.read()
        with io.BytesIO() as buffer:
            msa.write(buffer, "stockholm")
            actual = buffer.getvalue().decode()
        with open(sto) as f:
            expected = f.read()
        self.assertMultiLineEqual(actual, expected)

    def test_write_invalid_format(self):
        sto = os.path.join(self.formats_folder, "stockholm.1")
        with easel.MSAFile(sto, "stockholm") as msa_file:
            msa = msa_file.read()

        self.assertRaises(ValueError, msa.write, io.BytesIO(), "invalidformat")
