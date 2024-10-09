import copy
import io
import itertools
import os
import shutil
import unittest
import tempfile

import pyhmmer
from pyhmmer.errors import EaselError, AlphabetMismatch
from pyhmmer.easel import Alphabet, SequenceFile, DigitalMSA, TextMSA
from pyhmmer.plan7 import HMM, HMMFile, TraceAligner, Traces

from .. import __name__ as __package__
from .utils import HMMER_FOLDER, resource_files


@unittest.skipUnless(resource_files, "importlib.resources.files not available")
class TestTraces(unittest.TestCase):

    def setUp(self):
        hmm_path = resource_files(__package__).joinpath("data", "hmms", "txt", "LuxC.hmm")
        seqs_path = resource_files(__package__).joinpath("data", "seqs", "LuxC.faa")

        with HMMFile(hmm_path) as hmm_file:
            self.hmm = next(hmm_file)
        with SequenceFile(seqs_path, digital=True, alphabet=self.hmm.alphabet) as seqs_file:
            self.seqs = seqs_file.read_block()

        self.aligner = TraceAligner()
        self.traces = self.aligner.compute_traces(self.hmm, self.seqs)

    def test_eq(self):
        traces2 = self.aligner.compute_traces(self.hmm, self.seqs)
        self.assertSequenceEqual(self.traces, traces2)
        self.assertEqual(self.traces, traces2)

    def test_bool(self):
        self.assertTrue(self.traces)
        self.assertFalse(Traces())

    def test_index_error(self):
        with self.assertRaises(IndexError):
            self.traces[len(self.traces) + 1]
        with self.assertRaises(IndexError):
            self.traces[-len(self.traces) - 1]

    def test_index_negative(self):
        t1 = self.traces[len(self.traces) - 1]
        t2 = self.traces[-1]
        self.assertEqual(t1, t2)

    def test_trace_M(self):
        for trace in self.traces:
            self.assertEqual(trace.M, self.hmm.M)

    def test_trace_L(self):
        for trace, seq in zip(self.traces, self.seqs):
            self.assertEqual(trace.L, len(seq))
