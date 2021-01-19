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


class TestPipelinesearch(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.alphabet = Alphabet.amino()
        with SequenceFile(pkg_resources.resource_filename("tests", "data/seqs/938293.PRJEB85.HG003687.faa")) as f:
            cls.references = [x.digitize(cls.alphabet) for x in f]

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
