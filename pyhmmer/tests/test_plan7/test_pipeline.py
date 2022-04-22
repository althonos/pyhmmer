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
from pyhmmer.easel import Alphabet, SequenceFile, DigitalSequence, TextSequence, MSAFile, DigitalMSA
from pyhmmer.errors import AlphabetMismatch


class TestSearchPipeline(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.alphabet = Alphabet.amino()

        seq_path = pkg_resources.resource_filename("pyhmmer.tests", "data/seqs/938293.PRJEB85.HG003687.faa")
        with SequenceFile(seq_path, digital=True, alphabet=cls.alphabet) as f:
            cls.references = list(f)

        msa_path = pkg_resources.resource_filename("pyhmmer.tests", "data/msa/LuxC.sto")
        with MSAFile(msa_path, digital=True, alphabet=cls.alphabet) as msa_f:
            cls.msa = next(msa_f)

    def test_search_seq_alphabet_mismatch(self):
        pipeline = Pipeline(alphabet=Alphabet.dna())

        # mismatch between pipeline alphabet and database alphabet
        dsq1 = TextSequence(sequence="ATGC").digitize(pipeline.alphabet)
        self.assertRaises(AlphabetMismatch, pipeline.search_seq, dsq1, self.references)

        # mismatch between pipeline alphabet and query alphabet
        dsq2 = TextSequence(sequence="IRGIY").digitize(self.alphabet)
        self.assertRaises(AlphabetMismatch, pipeline.search_seq, dsq2, self.references)

        # check that all ref sequences are checked, not just the first one
        references = self.references.copy()
        references.append(TextSequence(sequence="ATGC").digitize(Alphabet.dna()))
        pipeline = Pipeline(alphabet=self.alphabet)
        self.assertRaises(AlphabetMismatch, pipeline.search_seq, dsq2, references)

    def test_search_msa_alphabet_mismatch(self):
        pipeline = Pipeline(alphabet=Alphabet.dna())

        # mismatch between pipeline alphabet and query alphabet
        self.assertRaises(AlphabetMismatch, pipeline.search_msa, self.msa, self.references)

        # mismatch between pipeline alphabet and database alphabet
        dsq = TextSequence(sequence="ATGC").digitize(pipeline.alphabet)
        msa = DigitalMSA(pipeline.alphabet, sequences=[dsq], name=b"test")
        self.assertRaises(AlphabetMismatch, pipeline.search_msa, msa, self.references)

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

    def test_Z(self):
        seq = TextSequence(sequence="IRGIYNIIKSVAEDIEIGIIPPSKDHVTISSFKSPRIADT")
        bg = Background(self.alphabet)
        hmm, _, _ = Builder(self.alphabet).build(seq.digitize(self.alphabet), bg)

        # when Z is None, use the number of target sequences
        pipeline = Pipeline(alphabet=self.alphabet)
        self.assertIs(pipeline.Z, None)
        hits = pipeline.search_hmm(hmm, self.references[:100])
        self.assertEqual(hits.Z, 100)
        self.assertIs(pipeline.Z, None)
        # clearing the pipeline will keep the Z number as None
        pipeline.clear()
        self.assertIs(pipeline.Z, None)

        # when Z is not None, use the given number
        pipeline = Pipeline(alphabet=self.alphabet, Z=25)
        self.assertEqual(pipeline.Z, 25)
        hits = pipeline.search_hmm(hmm, self.references[:100])
        self.assertEqual(pipeline.Z, 25)
        self.assertEqual(hits.Z, 25)
        # clearing the pipeline will keep the Z number
        pipeline.clear()
        self.assertEqual(pipeline.Z, 25)

    def test_bit_cutoffs(self):
        seq = TextSequence(sequence="IRGIYNIIKSVAEDIEIGIIPPSKDHVTISSFKSPRIADT")
        bg = Background(self.alphabet)
        hmm, _, _ = Builder(self.alphabet).build(seq.digitize(self.alphabet), bg)
        pipeline = Pipeline(alphabet=self.alphabet, bit_cutoffs="trusted")

        # without cutoffs, the pipeline is not going to be able to process the hmm
        with self.assertRaises(ValueError):
            hits = pipeline.search_hmm(hmm, self.references)

        # with low cutoffs, the domain is found
        hmm.cutoffs.trusted = (10.0, 10.0)
        hits = pipeline.search_hmm(hmm, self.references)
        self.assertEqual(len(hits), 1)
        self.assertLessEqual(hits[0].score, 75.0)

        # with higher cutoff, no domain is found
        hmm.cutoffs.trusted = (75.0, 75.0)
        hits = pipeline.search_hmm(hmm, self.references)
        self.assertEqual(len(hits), 0)


class TestScanPipeline(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.alphabet = Alphabet.amino()

        seq_path = pkg_resources.resource_filename("pyhmmer.tests", "data/seqs/938293.PRJEB85.HG003687.faa")
        with SequenceFile(seq_path, digital=True, alphabet=cls.alphabet) as f:
            cls.references = list(f)

        hmm_file = pkg_resources.resource_filename("pyhmmer.tests", "data/hmms/txt/t2pks.hmm")
        with HMMFile(hmm_file) as f:
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


class TestIteratePipeline(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.alphabet = Alphabet.amino()

        seq_path = pkg_resources.resource_filename("pyhmmer.tests", "data/seqs/938293.PRJEB85.HG003687.faa")
        with SequenceFile(seq_path, digital=True, alphabet=cls.alphabet) as f:
            cls.references = list(f)

        query_path = pkg_resources.resource_filename("pyhmmer.tests", "data/seqs/LuxC.faa")
        with SequenceFile(query_path, digital=True, alphabet=cls.alphabet) as f:
            cls.query = next(seq for seq in f if b"P12748" in seq.name)

    def test_iterate_seq(self):
        # check that `Pipeline.iterate_seq` produces consistent results
        # compared to a 'gold standard' run of `jackhmmer`.
        pipeline = Pipeline(alphabet=self.alphabet, incE=0.001, incdomE=0.001)
        search_iterator = pipeline.iterate_seq(self.query, self.references)

        for n in range(3):
            self.assertEqual(search_iterator.iteration, n)
            iteration = next(search_iterator)

            self.assertEqual(search_iterator.iteration, n+1)
            self.assertEqual(iteration.iteration, n+1)

            msa_path = pkg_resources.resource_filename("pyhmmer.tests", f"data/jackhmmer/p12748-{n+1}.sto")
            with MSAFile(msa_path, digital=True, alphabet=self.alphabet) as msa_file:
                ref_msa = msa_file.read()

            self.assertEqual(len(iteration.msa), len(ref_msa))
            self.assertEqual(len(iteration.msa.sequences), len(ref_msa.sequences))
            self.assertEqual(iteration.msa.name, ref_msa.name)

            self.assertEqual(
                {seq.name for seq in iteration.msa.sequences},
                {seq.name for seq in ref_msa.sequences},
            )

        self.assertEqual(iteration.iteration, 3)
        self.assertTrue(iteration.converged)
