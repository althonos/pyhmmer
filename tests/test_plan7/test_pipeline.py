import abc
import io
import itertools
import os
import unittest
import tempfile
import threading
import pkg_resources

import pyhmmer
from pyhmmer.plan7 import Background, Builder, Pipeline, HMMFile, TopHits
from pyhmmer.easel import Alphabet, SequenceFile, TextSequence
from pyhmmer.errors import AlphabetMismatch


class TestSearchPipeline(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.alphabet = Alphabet.amino()
        with SequenceFile(pkg_resources.resource_filename("tests", "data/seqs/938293.PRJEB85.HG003687.faa")) as f:
            f.set_digital(cls.alphabet)
            cls.references = list(f)

    def test_alphabet_mismatch(self):
        pipeline = Pipeline(alphabet=Alphabet.dna())

        # mismatch between pipeline alphabet and query alphabet
        dsq = TextSequence(sequence="IRGIY").digitize(self.alphabet)
        self.assertRaises(AlphabetMismatch, pipeline.search_seq, dsq, self.references)

        # mismatch between pipeline alphabet and database alphabet
        dsq = TextSequence(sequence="ATGC").digitize(pipeline.alphabet)
        self.assertRaises(AlphabetMismatch, pipeline.search_seq, dsq, self.references)

    def test_search_hmm(self):
        seq = TextSequence(sequence="IRGIYNIIKSVAEDIEIGIIPPSKDHVTISSFKSPRIADT")
        bg = Background(self.alphabet)
        hmm, _, _ = Builder(self.alphabet).build(seq.digitize(self.alphabet), bg)
        pipeline = Pipeline(alphabet=self.alphabet)
        hits = pipeline.search_hmm(hmm, self.references)
        self.assertEqual(len(hits), 1)

    def test_search_seq(self):
        seq = TextSequence(sequence="IRGIYNIIKSVAEDIEIGIIPPSKDHVTISSFKSPRIADT")
        pipeline = Pipeline(alphabet=self.alphabet)
        hits = pipeline.search_seq(seq.digitize(self.alphabet), self.references)
        self.assertEqual(len(hits), 1)


class TestScanPipeline(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.alphabet = Alphabet.amino()
        with SequenceFile(pkg_resources.resource_filename("tests", "data/seqs/938293.PRJEB85.HG003687.faa")) as f:
            f.set_digital(cls.alphabet)
            cls.references = list(f)
        with HMMFile(pkg_resources.resource_filename("tests", "data/hmms/txt/t2pks.hmm")) as f:
            cls.hmms = list(f)

    def test_alphabet_mismatch(self):
        pipeline = Pipeline(alphabet=Alphabet.dna())

        # mismatch between pipeline alphabet and query alphabet
        dsq = TextSequence(sequence="IRGIY").digitize(self.alphabet)
        self.assertRaises(AlphabetMismatch, pipeline.scan_seq, dsq, self.hmms)

        # mismatch between pipeline alphabet and database alphabet
        dsq = TextSequence(sequence="ATGC").digitize(pipeline.alphabet)
        self.assertRaises(AlphabetMismatch, pipeline.scan_seq, dsq, self.hmms)

    def test_scan_seq(self):
        seq = next(x for x in self.references if x.name == b"938293.PRJEB85.HG003687_188")
        pipeline = Pipeline(alphabet=self.alphabet)
        hits = pipeline.scan_seq(seq, self.hmms)
        self.assertEqual(len(hits), 6)  # number found with `hmmscan`
