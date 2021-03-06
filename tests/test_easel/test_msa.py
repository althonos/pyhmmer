import copy
import gc
import io
import os
import unittest
import tempfile

from pyhmmer import easel


class TestMSA(object):
    @classmethod
    def setUpClass(cls):
        cls.formats_folder = os.path.realpath(
            os.path.join(
                __file__, os.pardir, os.pardir, os.pardir, "vendor", "easel", "formats"
            )
        )

    def test_eq(self):
        sto = os.path.join(self.formats_folder, "stockholm.1")
        msa = self.read_msa(sto)
        msa2 = copy.copy(msa)
        self.assertEqual(msa, msa2)
        msa2.name = b"name"
        self.assertNotEqual(msa, msa2)

    def test_write_roundtrip_stockholm(self):
        sto = os.path.join(self.formats_folder, "stockholm.1")
        msa = self.read_msa(sto)
        with io.BytesIO() as buffer:
            msa.write(buffer, "stockholm")
            actual = buffer.getvalue().decode()
        with open(sto) as f:
            expected = f.read()
        self.assertMultiLineEqual(actual, expected)

    def test_write_invalid_format(self):
        sto = os.path.join(self.formats_folder, "stockholm.1")
        msa = self.read_msa(sto)
        self.assertRaises(ValueError, msa.write, io.BytesIO(), "invalidformat")

    def test_sequences(self):
        sto = os.path.join(self.formats_folder, "stockholm.1")
        msa = self.read_msa(sto)
        self.assertEqual(len(msa.sequences), 2)
        self.assertEqual(msa.sequences[0].name, b"seq1")
        self.assertEqual(msa.sequences[1].name, b"seq2")


class TestTextMSA(TestMSA, unittest.TestCase):

    @staticmethod
    def read_msa(sto):
        with easel.MSAFile(sto, "stockholm") as msa_file:
            return msa_file.read()


class TestDigitalMSA(TestMSA, unittest.TestCase):

    @staticmethod
    def read_msa(sto):
        with easel.MSAFile(sto, "stockholm") as msa_file:
            msa_file.set_digital(msa_file.guess_alphabet())
            return msa_file.read()
