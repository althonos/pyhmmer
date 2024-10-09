import copy
import io
import itertools
import os
import shutil
import unittest
import tempfile
import pickle

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
            gencode = GeneticCode(42)

    def test_translate_empty(self):
        gencode = GeneticCode()
        prot = gencode.translate(bytearray())
        self.assertEqual(len(prot), 0)
        self.assertEqual(prot, VectorU8())

    def test_length_error(self):
        gencode = GeneticCode()
        with self.assertRaises(ValueError):
            gencode.translate(bytearray([0, 0, 0, 1, 1]))

    def test_pickle(self):
        gencode = GeneticCode(11)
        gencode2 = pickle.loads(pickle.dumps(gencode))
        self.assertEqual(gencode.translation_table, gencode2.translation_table)
        self.assertEqual(gencode.nucleotide_alphabet, gencode2.nucleotide_alphabet)
        self.assertEqual(gencode.amino_alphabet, gencode2.amino_alphabet)