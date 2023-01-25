import copy
import io
import itertools
import os
import shutil
import unittest
import tempfile
import pkg_resources

import pyhmmer
from pyhmmer.errors import EaselError, AlphabetMismatch
from pyhmmer.easel import (
    Alphabet, 
    GeneticCode,
    TextSequence, 
    DigitalSequence, 
    DigitalSequenceBlock, 
    VectorU8
)

class TestGeneticCode(unittest.TestCase):

    def test_invalid_table(self):
        with self.assertRaises(ValueError):
            gencode = GeneticCode(translation_table=42)


class TestTranslate(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.nt = Alphabet.dna()
        cls.aa = Alphabet.amino()
        cls.gencode = GeneticCode(cls.nt, cls.aa)

    def test_empty(self):
        seq = DigitalSequence(self.nt)
        prot = self.gencode.translate(seq)
        self.assertEqual(len(prot), 0)
        self.assertEqual(prot.sequence, VectorU8())

    def test_alphabet_mismatch(self):
        seq = TextSequence(sequence="MRRKYL").digitize(self.aa)
        with self.assertRaises(AlphabetMismatch):
            self.gencode.translate(seq)

    def test_length_error(self):
        seq = TextSequence(sequence="ATGC").digitize(self.nt)
        with self.assertRaises(ValueError):
            self.gencode.translate(seq)


class TestTranslateBlock(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.nt = Alphabet.dna()
        cls.aa = Alphabet.amino()
        cls.gencode = GeneticCode(cls.nt, cls.aa)

    def test_empty(self):
        block = DigitalSequenceBlock(self.nt)
        prots = self.gencode.translate_block(block)
        self.assertEqual(len(prots), 0)

    def test_alphabet_mismatch(self):
        block = DigitalSequenceBlock(self.aa)
        with self.assertRaises(AlphabetMismatch):
            self.gencode.translate_block(block)

    def test_length_error(self):
        seq1 = TextSequence(sequence="ATG").digitize(self.nt)
        seq2 = TextSequence(sequence="ATGC").digitize(self.nt)
        block = DigitalSequenceBlock(self.nt, (seq1, seq2))
        with self.assertRaises(ValueError):
            self.gencode.translate_block(block)
