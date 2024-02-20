import copy
import io
import itertools
import os
import shutil
import unittest
import tempfile

import pyhmmer
from pyhmmer.errors import EaselError, AlphabetMismatch
from pyhmmer.easel import Alphabet, SequenceFile, DigitalMSA, TextMSA, DigitalSequenceBlock
from pyhmmer.plan7 import HMM, HMMFile, TraceAligner, Traces

from .. import __name__ as __package__
from .utils import HMMER_FOLDER, resource_files


@unittest.skipUnless(resource_files, "importlib.resources.files not available")
class TestTraceAligner(unittest.TestCase):

    def setUp(self):
        hmm_path = resource_files(__package__).joinpath("data", "hmms", "txt", "LuxC.hmm")
        seqs_path = resource_files(__package__).joinpath("data", "seqs", "LuxC.faa")

        with HMMFile(hmm_path) as hmm_file:
            self.hmm = next(hmm_file)
        with SequenceFile(seqs_path, digital=True, alphabet=self.hmm.alphabet) as seqs_file:
            self.seqs = seqs_file.read_block()

    def test_align_traces(self):
        # hmmalign tests/data/hmms/txt/LuxC.hmm tests/data/seqs/LuxC.faa
        aligner = TraceAligner()
        traces = aligner.compute_traces(self.hmm, self.seqs)
        self.assertEqual(len(traces), len(self.seqs))
        msa = aligner.align_traces(self.hmm, self.seqs, traces, all_consensus_cols=True)
        self.assertEqual(len(msa.sequences), len(self.seqs))
        self.assertEqual(len(msa), 567)

    def test_align_traces_trim(self):
        # hmmalign --trim tests/data/hmms/txt/LuxC.hmm tests/data/seqs/LuxC.faa
        aligner = TraceAligner()
        traces = aligner.compute_traces(self.hmm, self.seqs)
        self.assertEqual(len(traces), len(self.seqs))
        msa = aligner.align_traces(self.hmm, self.seqs, traces, all_consensus_cols=True, trim=True)
        self.assertEqual(len(msa.sequences), len(self.seqs))
        self.assertEqual(len(msa), 429)

    def test_align_traces_mismatch(self):
        aligner = TraceAligner()
        traces = Traces()
        self.assertRaises(ValueError, aligner.align_traces, self.hmm, self.seqs, traces)

    def test_align_traces_msa_type(self):
        aligner = TraceAligner()
        traces = Traces()
        seqs = DigitalSequenceBlock(self.hmm.alphabet)
        msa = aligner.align_traces(self.hmm, seqs, traces)
        self.assertIsInstance(msa, TextMSA)
        msa_d = aligner.align_traces(self.hmm, seqs, traces, digitize=True)
        self.assertIsInstance(msa_d, DigitalMSA)
