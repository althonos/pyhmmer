import abc
import io
import itertools
import os
import unittest
import random
import tempfile
import threading
import pkg_resources

import pyhmmer
from pyhmmer.easel import Alphabet, SequenceFile, DigitalSequence, DigitalSequenceBlock, TextSequence, TextSequenceBlock, MSAFile, DigitalMSA
from pyhmmer.plan7 import Background, Builder, Pipeline, HMMFile, TopHits, OptimizedProfileBlock, Profile, LongTargetsPipeline
from pyhmmer.errors import AlphabetMismatch


class TestSearchPipeline(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.alphabet = Alphabet.amino()

        cls.reference_path = pkg_resources.resource_filename("pyhmmer.tests", "data/seqs/938293.PRJEB85.HG003687.faa")
        with SequenceFile(cls.reference_path, digital=True, alphabet=cls.alphabet) as f:
            cls.references = f.read_block()

        cls.msa_path = pkg_resources.resource_filename("pyhmmer.tests", "data/msa/LuxC.sto")
        with MSAFile(cls.msa_path, digital=True, alphabet=cls.alphabet) as msa_f:
            cls.msa = msa_f.read()

    def test_search_hmm_alphabet_mismatch(self):
        rng = pyhmmer.easel.Randomness(0)
        pipeline = Pipeline(alphabet=Alphabet.dna())

        # mismatch between pipeline alphabet and database alphabet
        hmm1 = pyhmmer.plan7.HMM.sample(100, pipeline.alphabet, rng)
        self.assertRaises(AlphabetMismatch, pipeline.search_hmm, hmm1, self.references)

        # mismatch between pipeline alphabet and query alphabet
        hmm2 = pyhmmer.plan7.HMM.sample(100, self.alphabet, rng)
        self.assertRaises(AlphabetMismatch, pipeline.search_hmm, hmm2, self.references)

    def test_search_seq_alphabet_mismatch(self):
        pipeline = Pipeline(alphabet=Alphabet.dna())

        # mismatch between pipeline alphabet and database alphabet
        dsq1 = TextSequence(sequence="ATGC").digitize(pipeline.alphabet)
        self.assertRaises(AlphabetMismatch, pipeline.search_seq, dsq1, self.references)

        # mismatch between pipeline alphabet and query alphabet
        dsq2 = TextSequence(sequence="IRGIY").digitize(self.alphabet)
        self.assertRaises(AlphabetMismatch, pipeline.search_seq, dsq2, self.references)

    def test_search_msa_alphabet_mismatch(self):
        pipeline = Pipeline(alphabet=Alphabet.dna())

        # mismatch between pipeline alphabet and query alphabet
        self.assertRaises(AlphabetMismatch, pipeline.search_msa, self.msa, self.references)

        # mismatch between pipeline alphabet and database alphabet
        dsq = TextSequence(sequence="ATGC").digitize(pipeline.alphabet)
        msa = DigitalMSA(pipeline.alphabet, sequences=[dsq], name=b"test")
        self.assertRaises(AlphabetMismatch, pipeline.search_msa, msa, self.references)

    def test_search_hmm_block(self):
        seq = TextSequence(sequence="IRGIYNIIKSVAEDIEIGIIPPSKDHVTISSFKSPRIADT", name=b"seq1")
        bg = Background(self.alphabet)
        hmm, _, _ = Builder(self.alphabet).build(seq.digitize(self.alphabet), bg)
        pipeline = Pipeline(alphabet=self.alphabet)
        hits = pipeline.search_hmm(hmm, self.references)
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits.query_name, hmm.name)
        self.assertEqual(hits.query_accession, hmm.accession)

    def test_search_hmm_file(self):
        seq = TextSequence(sequence="IRGIYNIIKSVAEDIEIGIIPPSKDHVTISSFKSPRIADT", name=b"seq1")
        bg = Background(self.alphabet)
        hmm, _, _ = Builder(self.alphabet).build(seq.digitize(self.alphabet), bg)
        pipeline = Pipeline(alphabet=self.alphabet)
        with SequenceFile(self.reference_path, digital=True, alphabet=self.alphabet) as seqs_file:
            hits = pipeline.search_hmm(hmm, seqs_file)
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits.query_name, hmm.name)
        self.assertEqual(hits.query_accession, hmm.accession)

    def test_search_hmm_unnamed(self):
        # make sure `Pipeline.search_hmm` doesn't crash when given an HMM with no name
        rng = pyhmmer.easel.Randomness()
        hmm = pyhmmer.plan7.HMM.sample(100, self.alphabet, rng)
        hmm.name = hmm.accession = None
        pipeline = Pipeline(alphabet=self.alphabet)
        hits = pipeline.search_hmm(hmm, self.references)
        self.assertIs(hits.query_name, None)
        self.assertIs(hits.query_accession, None)

    def test_search_seq_block(self):
        seq = TextSequence(sequence="IRGIYNIIKSVAEDIEIGIIPPSKDHVTISSFKSPRIADT", name=b"seq1", accession=b"SQ001")
        pipeline = Pipeline(alphabet=self.alphabet)
        hits = pipeline.search_seq(seq.digitize(self.alphabet), self.references)
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits.query_name, seq.name)
        # self.assertEqual(hits.query_accession, seq.accession) # NOTE: p7_SingleBuilder doesn't copy the accession...

    def test_search_seq_file(self):
        seq = TextSequence(sequence="IRGIYNIIKSVAEDIEIGIIPPSKDHVTISSFKSPRIADT", name=b"seq1", accession=b"SQ001")
        pipeline = Pipeline(alphabet=self.alphabet)
        with SequenceFile(self.reference_path, digital=True, alphabet=self.alphabet) as seqs_file:
            hits = pipeline.search_seq(seq.digitize(self.alphabet), seqs_file)
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits.query_name, seq.name)
        # self.assertEqual(hits.query_accession, seq.accession) # NOTE: p7_SingleBuilder doesn't copy the accession...

    def test_search_msa_block(self):
        pipeline = Pipeline(alphabet=self.alphabet)
        hits = pipeline.search_msa(self.msa, self.references)
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits.query_name, self.msa.name)
        self.assertEqual(hits.query_accession, self.msa.accession)

    def test_search_msa_file(self):
        pipeline = Pipeline(alphabet=self.alphabet)
        with SequenceFile(self.reference_path, digital=True, alphabet=self.alphabet) as seqs_file:
            hits = pipeline.search_msa(self.msa, seqs_file)
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits.query_name, self.msa.name)
        self.assertEqual(hits.query_accession, self.msa.accession)

    def test_Z(self):
        seq = TextSequence(sequence="IRGIYNIIKSVAEDIEIGIIPPSKDHVTISSFKSPRIADT", name=b"seq1")
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
            cls.references = f.read_block()

        hmm_file = pkg_resources.resource_filename("pyhmmer.tests", "data/hmms/txt/t2pks.hmm")
        with HMMFile(hmm_file) as f:
            cls.hmms = list(f)

    def test_scan_seq_alphabet_mismatch(self):
        pipeline = Pipeline(alphabet=Alphabet.dna())

        # mismatch between pipeline alphabet and query alphabet
        dsq = TextSequence(sequence="IRGIY").digitize(self.alphabet)
        oprofiles = OptimizedProfileBlock(pipeline.alphabet)
        self.assertRaises(AlphabetMismatch, pipeline.scan_seq, dsq, oprofiles)

        # mismatch between pipeline alphabet and database alphabet
        dsq = TextSequence(sequence="ATGC").digitize(pipeline.alphabet)
        oprofiles = OptimizedProfileBlock(self.alphabet)
        self.assertRaises(AlphabetMismatch, pipeline.scan_seq, dsq, oprofiles)

    def test_scan_seq_block(self):
        seq = next(x for x in self.references if x.name == b"938293.PRJEB85.HG003687_188")

        oprofiles = OptimizedProfileBlock(seq.alphabet)
        background = Background(seq.alphabet)
        for hmm in self.hmms:
            profile = Profile(hmm.M, hmm.alphabet)
            profile.configure(hmm, background)
            oprofiles.append(profile.optimized())

        pipeline = Pipeline(alphabet=self.alphabet)
        hits = pipeline.scan_seq(seq, oprofiles)
        self.assertEqual(len(hits), 6)  # number found with `hmmscan`

    def test_scan_seq_file(self):
        seq = next(x for x in self.references if x.name == b"938293.PRJEB85.HG003687_188")

        hmm_file = pkg_resources.resource_filename("pyhmmer.tests", "data/hmms/db/t2pks.hmm")
        with HMMFile(hmm_file) as f:
            pipeline = Pipeline(alphabet=self.alphabet)
            hits = pipeline.scan_seq(seq, f.optimized_profiles())

        self.assertEqual(len(hits), 6)  # number found with `hmmscan`


class TestIteratePipeline(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.alphabet = Alphabet.amino()

        seq_path = pkg_resources.resource_filename("pyhmmer.tests", "data/seqs/938293.PRJEB85.HG003687.faa")
        with SequenceFile(seq_path, digital=True, alphabet=cls.alphabet) as f:
            cls.references = f.read_block()

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


class TestLongTargetsPipeline(unittest.TestCase):

    @staticmethod
    def _random_sequence(alphabet, name, length=5000):
        symbols = alphabet.symbols[:alphabet.K]
        return TextSequence(
            name=name,
            sequence="".join(random.sample(symbols, 1)[0] for i in range(length))
        ).digitize(alphabet)

    def test_search_hmm_block(self):
        dna = Alphabet.dna()
        rng = pyhmmer.easel.Randomness(0)

        targets = DigitalSequenceBlock(dna, [
            self._random_sequence(dna, b"seq1"),
            self._random_sequence(dna, b"seq2"),
        ])

        hmm = pyhmmer.plan7.HMM.sample(100, dna, rng)
        hmm.name = b"test_one"
        hmm.accession = b"TST001"

        pipeline = LongTargetsPipeline(alphabet=dna)
        hits = pipeline.search_hmm(hmm, targets)
        self.assertEqual(hits.query_name, hmm.name)
        self.assertEqual(hits.query_accession, hmm.accession)

    def test_search_hmm_file(self):
        dna = Alphabet.dna()
        rng = pyhmmer.easel.Randomness(0)

        hmm = pyhmmer.plan7.HMM.sample(100, dna, rng)
        hmm.name = b"test_one"
        hmm.accession = b"TST001"

        with tempfile.NamedTemporaryFile(mode="wb", suffix=".fna") as f:
            self._random_sequence(dna, b"seq1").write(f)
            self._random_sequence(dna, b"seq2").write(f)
            f.flush()

            with pyhmmer.easel.SequenceFile(f.name, digital=True, alphabet=dna) as targets:
                pipeline = LongTargetsPipeline(alphabet=dna)
                hits = pipeline.search_hmm(hmm, targets)

        self.assertEqual(hits.query_name, hmm.name)
        self.assertEqual(hits.query_accession, hmm.accession)

    def test_search_hmm_alphabet_mismatch(self):
        dna = Alphabet.dna()
        amino = Alphabet.amino()
        rng = pyhmmer.easel.Randomness(0)

        # mismatch between pipeline alphabet and database alphabet
        pipeline = LongTargetsPipeline(alphabet=dna)
        targets = DigitalSequenceBlock(amino)
        hmm = pyhmmer.plan7.HMM.sample(100, dna, rng)
        self.assertRaises(AlphabetMismatch, pipeline.search_hmm, hmm, targets)

        # mismatch between pipeline alphabet and query alphabet
        pipeline = LongTargetsPipeline(alphabet=dna)
        targets = DigitalSequenceBlock(dna)
        hmm = pyhmmer.plan7.HMM.sample(100, amino, rng)
        self.assertRaises(AlphabetMismatch, pipeline.search_hmm, hmm, targets)