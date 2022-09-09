import abc
import copy
import gc
import io
import os
import unittest
import tempfile
import warnings

from pyhmmer import easel


class _TestSequenceBase(abc.ABC):

    @classmethod
    @abc.abstractmethod
    def Sequence(cls, **kwargs):
        pass

    def test_init_empty(self):
        seq = self.Sequence()
        self.assertEqual(seq.name, b"")
        self.assertEqual(len(seq), 0)
        self.assertEqual(len(seq.sequence), 0)

    def test_setter_name(self):
        seq = self.Sequence()
        self.assertEqual(seq.name, b"")
        seq.name = b"OTHER"
        self.assertEqual(seq.name, b"OTHER")

    def test_setter_description(self):
        seq = self.Sequence()
        self.assertIs(seq.description, b"")
        seq.description = b"a test sequence"
        self.assertEqual(seq.description, b"a test sequence")


class TestSequence(unittest.TestCase):

    def test_init_abstract(self):
        self.assertRaises(TypeError, easel.Sequence, name=b"ecoRI", sequence="GAATTC")

    def test_write_roundtrip(self):
        seq = easel.TextSequence(name=b"ecoRI", sequence="GAATTC")
        with io.BytesIO() as buffer:
            seq.write(buffer)
            seq2 = easel.SequenceFile.parse(buffer.getvalue(), "fasta")
        self.assertEqual(seq, seq2)


class TestDigitalSequence(_TestSequenceBase, unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.abc = easel.Alphabet.dna()

    @classmethod
    def Sequence(cls, **kwargs):
        return easel.DigitalSequence(cls.abc, **kwargs)

    def test_init_kwargs(self):
        arr = bytearray([0, 1, 2, 3])
        seq = easel.DigitalSequence(self.abc, name=b"TEST", sequence=arr)
        self.assertEqual(seq.name, b"TEST")
        self.assertEqual(bytearray(seq.sequence), arr)
        self.assertEqual(len(seq), 4)

    def test_copy(self):
        arr = bytearray([0, 1, 2, 3])
        seq = easel.DigitalSequence(self.abc, name=b"TEST", sequence=arr)

        # copy the sequence
        cpy = seq.copy()

        # check the copy is equal and attributes are equal
        self.assertEqual(seq, cpy)
        self.assertEqual(seq.name, cpy.name)
        self.assertEqual(seq.checksum(), cpy.checksum())

        # check attributes can be read even after the original object
        # is (hopefully) deallocated, to make sure internal strings were copied
        del seq
        gc.collect()
        self.assertEqual(cpy.name, b"TEST")

        # check `copy.copy` works too
        cpy2 = copy.copy(cpy)
        self.assertEqual(cpy, cpy2)
        self.assertEqual(cpy.name, cpy2.name)

    def test_eq(self):
        arr = bytearray([0, 1, 2, 3])
        seq1 = easel.DigitalSequence(self.abc, name=b"TEST", sequence=arr)
        self.assertEqual(seq1, seq1)
        self.assertNotEqual(seq1, object())

        seq2 = easel.DigitalSequence(self.abc, name=b"TEST", sequence=arr)
        self.assertEqual(seq1, seq2)

        seq3 = easel.DigitalSequence(self.abc, name=b"OTHER", sequence=arr)
        self.assertNotEqual(seq1, seq3)

        arr2 = bytearray([0, 1, 2, 0])
        seq4 = easel.DigitalSequence(self.abc, name=b"TEST", sequence=arr2)
        self.assertNotEqual(seq1, seq4)

    def test_eq_not_seq(self):
        arr = bytearray([0, 1, 2, 3])
        seq1 = easel.DigitalSequence(self.abc, name=b"TEST", sequence=arr)
        self.assertNotEqual(seq1, 1)
        self.assertNotEqual(seq1, b"hello")

    def test_textize_roundtrip(self):
        arr = bytearray([0, 1, 2, 3])
        dsq1 = easel.DigitalSequence(self.abc, name=b"TEST", sequence=arr)
        tsq1 = dsq1.textize()
        self.assertEqual(tsq1.name, dsq1.name)

        dsq2 = tsq1.digitize(self.abc)
        self.assertEqual(dsq1.name, dsq2.name)
        self.assertEqual(list(dsq1.sequence), list(dsq2.sequence))
        self.assertEqual(dsq1, dsq2)

    def test_reverse_complement(self):
        seq = easel.TextSequence(sequence="ATGC").digitize(self.abc)
        rc = seq.reverse_complement()
        self.assertEqual(rc.textize().sequence, "GCAT")

    def test_reverse_complement_inplace(self):
        seq = easel.TextSequence(sequence="ATGC").digitize(self.abc)
        seq.reverse_complement(inplace=True)
        self.assertEqual(seq.textize().sequence, "GCAT")

    def test_reverse_complement_protein(self):
        abc = easel.Alphabet.amino()
        seq = easel.TextSequence(sequence="MEMLP").digitize(abc)
        self.assertRaises(ValueError, seq.reverse_complement)
        self.assertRaises(ValueError, seq.reverse_complement, inplace=True)

    def test_invalid_characters(self):
        self.assertRaises(ValueError, easel.DigitalSequence, self.abc, name=b"TEST", sequence=b"test")


class TestTextSequence(_TestSequenceBase, unittest.TestCase):

    @classmethod
    def Sequence(cls, **kwargs):
        return easel.TextSequence(**kwargs)

    def test_init_kwargs(self):
        seq = easel.TextSequence(name=b"TEST", sequence="ATGC")
        self.assertEqual(seq.name, b"TEST")
        self.assertEqual(seq.sequence, "ATGC")
        self.assertEqual(len(seq), 4)

    def test_copy(self):
        seq = easel.TextSequence(name=b"TEST", sequence="ATGC")

        # copy the sequence
        cpy = seq.copy()

        # check the copy is equal and attributes are equal
        self.assertEqual(seq, cpy)
        self.assertEqual(seq.name, cpy.name)
        self.assertEqual(seq.checksum(), cpy.checksum())

        # check attributes can be read even after the original object
        # is (hopefully) deallocated, to make sure internal strings were copied
        del seq
        gc.collect()
        self.assertEqual(cpy.name, b"TEST")

        # check `copy.copy` works too
        cpy2 = copy.copy(cpy)
        self.assertEqual(cpy, cpy2)
        self.assertEqual(cpy.name, cpy2.name)

    def test_eq(self):
        seq1 = easel.TextSequence(name=b"TEST", sequence="ATGC")
        self.assertEqual(seq1, seq1)
        self.assertNotEqual(seq1, object())

        seq2 = easel.TextSequence(name=b"TEST", sequence="ATGC")
        self.assertEqual(seq1, seq2)

        seq3 = easel.TextSequence(name=b"OTHER", sequence="ATGC")
        self.assertNotEqual(seq1, seq3)

        seq4 = easel.TextSequence(name=b"TEST", sequence="ATGT")
        self.assertNotEqual(seq1, seq4)

    def test_eq_not_seq(self):
        seq1 = easel.TextSequence(name=b"TEST", sequence="ATGC")
        self.assertNotEqual(seq1, 1)
        self.assertNotEqual(seq1, b"hello")

    def test_digitize_roundtrip(self):
        seq1 = easel.TextSequence(name=b"TEST", sequence="ATGC")
        dsq1 = seq1.digitize(easel.Alphabet.dna())
        self.assertEqual(seq1.name, dsq1.name)

        seq2 = dsq1.textize()
        self.assertEqual(seq1.name, seq2.name)
        self.assertEqual(seq1, seq2)

    def test_digitize_invalid(self):
        seq1 = easel.TextSequence(name=b"TEST", sequence=">v3ry 1nv4l1d")
        self.assertRaises(ValueError, seq1.digitize, easel.Alphabet.dna())
        self.assertRaises(ValueError, seq1.digitize, easel.Alphabet.amino())

    def test_reverse_complement(self):
        seq = easel.TextSequence(sequence="ATGC")
        rc = seq.reverse_complement()
        self.assertEqual(rc.sequence, "GCAT")

    def test_reverse_complement_inplace(self):
        seq = easel.TextSequence(sequence="ATGC")
        seq.reverse_complement(inplace=True)
        self.assertEqual(seq.sequence, "GCAT")

    def test_reverse_complement_protein(self):
        seq = easel.TextSequence(sequence="MEMLP")

        with warnings.catch_warnings():
            warnings.simplefilter('error')
            self.assertRaises(Warning, seq.reverse_complement)

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            rc = seq.reverse_complement()
            self.assertEqual(rc.sequence, "NNKNK")

    def test_setter_residue_markup(self):
        seq = easel.TextSequence(sequence="MEMLP")
        self.assertEqual(seq.residue_markups, {})
        seq.residue_markups = {b"kind": b"DADDC"}
        self.assertEqual(seq.residue_markups[b"kind"], b"DADDC")
        with self.assertRaises(ValueError):
            seq.residue_markups = {b"kind": b"D"}
