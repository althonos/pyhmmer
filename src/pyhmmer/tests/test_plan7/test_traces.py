import copy
import io
import itertools
import os
import shutil
import unittest
import tempfile

import pyhmmer
from pyhmmer.errors import EaselError, AlphabetMismatch
from pyhmmer.easel import Alphabet, SequenceFile, TextSequence, DigitalMSA, TextMSA
from pyhmmer.plan7 import HMM, HMMFile, TraceAligner, Traces, Trace

from .. import __name__ as __package__
from .utils import HMMER_FOLDER, resource_files


class TestTraces(unittest.TestCase):

    def _new_trace(self, n=3):
        seq = TextSequence(sequence="N"*n)
        trace = Trace.from_sequence(seq)
        return trace

    def _new_block(self, traces=()):
        return Traces(traces)

    @unittest.skipUnless(resource_files, "importlib.resources.files not available")
    def test_aligner_traces(self):
        hmm_path = resource_files(__package__).joinpath("data", "hmms", "txt", "LuxC.hmm")
        seqs_path = resource_files(__package__).joinpath("data", "seqs", "LuxC.faa")

        with HMMFile(hmm_path) as hmm_file:
            hmm = next(hmm_file)
        with SequenceFile(seqs_path, digital=True, alphabet=hmm.alphabet) as seqs_file:
            seqs = seqs_file.read_block()

        aligner = TraceAligner()
        traces = aligner.compute_traces(hmm, seqs)

        for trace, seq in zip(traces, seqs):
            self.assertEqual(trace.M, hmm.M)
            self.assertEqual(trace.L, len(seq))

    def test_bool(self):
        self.assertTrue(not self._new_block())
        self.assertFalse(self._new_block())
        t1 = self._new_trace()
        block = self._new_block([t1])
        self.assertTrue(block)

    def test_identity(self):
        self.assertIsNot(self._new_block(), self._new_block())

    def test_len(self):
        t1 = self._new_trace(1)
        t2 = self._new_trace(2)
        t3 = self._new_trace(3)

        self.assertEqual(len(self._new_block([])), 0)
        self.assertEqual(len(self._new_block([t1])), 1)
        self.assertEqual(len(self._new_block([t1, t2, t3])), 3)

    def test_append(self):
        block = self._new_block()
        self.assertEqual(len(block), 0)
        t1 = self._new_trace(1)
        block.append(t1)
        self.assertEqual(len(block), 1)
        self.assertEqual(block[0], t1)
        t2 = self._new_trace(2)
        block.append(t2)
        self.assertEqual(len(block), 2)
        self.assertEqual(block[1], t2)

    def test_iter(self):
        t1 = self._new_trace(1)
        t2 = self._new_trace(2)
        t3 = self._new_trace(3)

        block = self._new_block([t1, t2, t3])
        self.assertEqual(len(block), 3)

        iterator = iter(block)
        self.assertIs(next(iterator), t1)
        self.assertIs(next(iterator), t2)
        self.assertIs(next(iterator), t3)
        self.assertRaises(StopIteration, next, iterator)

    def test_clear(self):
        t1 = self._new_trace(1)
        t2 = self._new_trace(2)

        block = self._new_block([t1, t2])
        self.assertEqual(len(block), 2)
        block.clear()
        self.assertEqual(len(block), 0)
        block.clear()
        self.assertEqual(len(block), 0)

        block = self._new_block()
        self.assertEqual(len(block), 0)
        block.clear()
        self.assertEqual(len(block), 0)

    def test_getitem(self):
        t1 = self._new_trace(1)
        t2 = self._new_trace(2)

        block = self._new_block([t1, t2])
        self.assertEqual(len(block), 2)

        self.assertIs(block[0], t1)
        self.assertIs(block[1], t2)
        self.assertIs(block[-1], t2)
        self.assertIs(block[-2], t1)

    def test_getitem_indexerror(self):
        t1 = self._new_trace()
        t2 = self._new_trace()
        block = self._new_block([t1, t2])

        with self.assertRaises(IndexError):
            block[3]
        with self.assertRaises(IndexError):
            block[-5]

    def test_setitem(self):
        t1 = self._new_trace(1)
        t2 = self._new_trace(2)
        t3 = self._new_trace(3)

        block = self._new_block([t1, t2])
        self.assertEqual(len(block), 2)

        block[0] = t3
        self.assertIs(block[0], t3)
        self.assertIs(block[1], t2)

    def test_remove(self):
        t1 = self._new_trace(1)
        t2 = self._new_trace(2)

        block = self._new_block([t1, t2])
        block.remove(t1)
        self.assertEqual(len(block), 1)
        self.assertIs(block[0], t2)

    def test_index(self):
        t1 = self._new_trace(1)
        t2 = self._new_trace(2)
        t3 = self._new_trace(3)

        block = self._new_block([t1, t2])
        self.assertEqual(len(block), 2)

        self.assertEqual(block.index(t1), 0)
        self.assertEqual(block.index(t2), 1)

        self.assertRaises(ValueError, block.index, t3)
        self.assertRaises(ValueError, block.index, t1, start=1)

    def test_pop(self):
        t1 = self._new_trace(1)
        t2 = self._new_trace(2)
        t3 = self._new_trace(3)

        block = self._new_block([t1, t2, t3])
        self.assertEqual(len(block), 3)

        self.assertIs(block.pop(), t3)
        self.assertIs(block.pop(0), t1)

        self.assertEqual(len(block), 1)
        self.assertIs(block[0], t2)

        self.assertIs(block.pop(-1), t2)
        self.assertEqual(len(block), 0)
        self.assertRaises(IndexError, block.pop)

        empty = self._new_block()
        self.assertRaises(IndexError, block.pop)

    def test_setitem_slice(self):
        t1 = self._new_trace(1)
        t2 = self._new_trace(2)
        t3 = self._new_trace(3)

        block = self._new_block([t1, t2])
        self.assertEqual(len(block), 2)
        self.assertIs(block[0], t1)
        self.assertIs(block[1], t2)

        block[:] = [t1, t2, t3]
        self.assertEqual(len(block), 3)
        self.assertIs(block[0], t1)
        self.assertIs(block[1], t2)
        self.assertIs(block[2], t3)

        block[2:5] = [t3, t3, t3]
        self.assertEqual(len(block), 5)
        self.assertIs(block[3], t3)
        self.assertIs(block[4], t3)

    # def test_delitem_slice(self):
    #     t1 = self._new_trace(1)
    #     t2 = self._new_trace(2)
    #     t3 = self._new_trace(3)

    #     block = self._new_block([t1, t2, t3])
    #     self.assertEqual(len(block), 3)
    #     del block[:]
    #     self.assertEqual(len(block), 0)

    #     block = self._new_block([t1, t2, t3])
    #     self.assertEqual(len(block), 3)
    #     del block[::2]
    #     self.assertEqual(len(block), 1)
    #     self.assertIs(block[0], t2)

    def test_contains(self):
        t1 = self._new_trace(1)
        t2 = self._new_trace(2)
        t3 = self._new_trace(3)

        block = self._new_block([t1, t2])
        self.assertEqual(len(block), 2)
        self.assertIn(t1, block)
        self.assertIn(t2, block)
        self.assertNotIn(t3, block)

        self.assertNotIn(42, block)
        self.assertNotIn(object(), block)

    # def test_copy(self):
    #     t1 = self._new_trace()
    #     t2 = self._new_trace()
    #     t3 = self._new_trace()

    #     block = self._new_block([t1, t2, t3])
    #     block_copy = block.copy()

    #     self.assertEqual(len(block), len(block_copy))
    #     self.assertEqual(block, block_copy)

    #     self.assertIs(block[0], block_copy[0])
    #     self.assertIs(block[1], block_copy[1])
    #     self.assertIs(block[2], block_copy[2])

    #     block.remove(t3)
    #     self.assertEqual(len(block), 2)
    #     self.assertEqual(len(block_copy), 3)

    # def test_pickle(self):
    #     t1 = self._new_trace(1)
    #     t2 = self._new_trace(2)
    #     t3 = self._new_trace(3)

    #     block = self._new_block([t1, t2, t3])
    #     block_pickled = pickle.loads(pickle.dumps(block))
        
    #     self.assertEqual(len(block), len(block_pickled))
    #     self.assertEqual(block[0], block_pickled[0])
    #     self.assertEqual(block[1], block_pickled[1])
    #     self.assertEqual(block[2], block_pickled[2])
