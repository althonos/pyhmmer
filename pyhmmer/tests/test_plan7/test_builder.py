import copy
import datetime
import io
import itertools
import os
import shutil
import sys
import unittest
import tempfile
import pkg_resources
from unittest import mock

import pyhmmer
from pyhmmer.errors import EaselError, AlphabetMismatch
from pyhmmer.easel import Alphabet, SequenceFile, MSAFile
from pyhmmer.plan7 import HMMFile, Pipeline, Builder, Background

from ..utils import HMMER_FOLDER

class _TestBuilderBase(object):

    @classmethod
    def setUpClass(cls):
        cls.testdata = os.path.join(HMMER_FOLDER, "testsuite")

        cls.faa_path = pkg_resources.resource_filename("pyhmmer.tests", "data/seqs/PKSI.faa")
        with SequenceFile(cls.faa_path, digital=True) as seqf:
            cls.proteins = list(seqf)

        if os.path.exists(cls.testdata):
            cls.ecori_fa = os.path.join(cls.testdata, "ecori.fa")
            with SequenceFile(cls.ecori_fa) as seqf:
                cls.dna = next(seqf)

            cls.ecori_hmm = cls.ecori_fa.replace(".fa", ".hmm")
            with HMMFile(cls.ecori_hmm) as hmmf:
                cls.ecori = next(hmmf)

            cls.msa_path = os.path.join(cls.testdata, "3box.sto")
            with MSAFile(cls.msa_path, digital=True) as msa_file:
                cls.msa = next(msa_file)
                cls.msa.name = b"3box"


class TestBuilder(unittest.TestCase):

    def test_invalid_window_length(self):
        abc = Alphabet.amino()
        builder = Builder(alphabet=abc)
        with self.assertRaises(ValueError):
            builder.window_length = 2

    def test_init_error(self):
        abc = Alphabet.amino()
        self.assertRaises(ValueError, Builder, alphabet=abc, architecture="nonsense")
        self.assertRaises(ValueError, Builder, alphabet=abc, weighting="nonsense")
        self.assertRaises(ValueError, Builder, alphabet=abc, effective_number="nonsense")
        self.assertRaises(ValueError, Builder, alphabet=abc, prior_scheme="nonsense")

    def test_init_effective_number(self):
        abc = Alphabet.amino()
        builder = Builder(alphabet=abc, effective_number=2.0)
        self.assertEqual(builder.effective_number, 2.0)
        self.assertRaises(TypeError, Builder, alphabet=abc, effective_number=[])

    def test_copy(self):
        abc = Alphabet.amino()
        builder = Builder(alphabet=abc)
        bc1 = builder.copy()
        self.assertEqual(builder.alphabet, bc1.alphabet)
        self.assertEqual(builder.seed, bc1.seed)
        bc2 = copy.copy(bc1)
        self.assertEqual(builder.alphabet, bc2.alphabet)
        self.assertEqual(builder.seed, bc2.seed)


class TestBuilderSingle(_TestBuilderBase, unittest.TestCase):

    def test_alphabet_mismatch(self):
        abc = Alphabet.dna()
        amino = Alphabet.amino()

        builder = Builder(alphabet=abc)
        bg = Background(alphabet=abc)

        seq = self.proteins[0]
        self.assertRaises(AlphabetMismatch, builder.build, seq, bg)

    def test_command_line(self):
        abc = Alphabet.amino()
        builder = Builder(alphabet=abc)
        seq = self.proteins[1]

        argv = ["mycommand", "--param", "1"]
        with mock.patch.object(sys, "argv", argv):
            hmm, profile, opt = builder.build(seq, Background(abc))
        self.assertEqual(hmm.command_line, " ".join(argv))

    def test_creation_time(self):
        abc = Alphabet.amino()
        builder = Builder(alphabet=abc)
        seq = self.proteins[1]

        hmm, profile, opt = builder.build(seq, Background(abc))
        self.assertIsInstance(hmm.creation_time, datetime.datetime)

    def test_protein(self):
        abc = Alphabet.amino()
        builder = Builder(alphabet=abc)
        seq = self.proteins[1]
        hmm, profile, opt = builder.build(seq, Background(abc))
        self.assertEqual(hmm.name, seq.name)
        self.assertEqual(hmm.alphabet, abc)
        self.assertEqual(profile.alphabet, abc)
        self.assertEqual(opt.alphabet, abc)

    @unittest.skipUnless(os.path.exists(HMMER_FOLDER), "test data not available")
    def test_dna(self):
        abc = Alphabet.dna()
        builder = Builder(alphabet=abc)
        hmm, profile, opt = builder.build(self.dna.digitize(abc), Background(abc))
        self.assertEqual(hmm.name, self.dna.name)
        self.assertEqual(hmm.alphabet, abc)
        self.assertEqual(profile.alphabet, abc)
        self.assertEqual(opt.alphabet, abc)
        self.assertEqual(hmm.M, self.ecori.M)

    @unittest.skipUnless(os.path.exists(HMMER_FOLDER), "test data not available")
    def test_score_system_change(self):
        abc = Alphabet.amino()
        bg = Background(abc)
        builder = Builder(alphabet=abc)
        seq = self.proteins[1]
        # build first HMM with mx=BLOSUM62, popen=0.02, pextend=0.4
        self.assertEqual(builder.score_matrix, "BLOSUM62")
        hmm1, _, _ = builder.build(seq, bg)
        hmm1.command_line = hmm1.creation_time = None
        # build second HMM with mx=BLOSUM90, popen=0.02, pextend=0.4
        # --> the transition probabilities won't change, match emissions will
        builder.score_matrix = "BLOSUM90"
        self.assertEqual(builder.score_matrix, "BLOSUM90")
        hmm2, _, _ = builder.build(seq, bg)
        hmm2.command_line = hmm1.creation_time = None
        self.assertEqual(hmm1.transition_probabilities, hmm2.transition_probabilities)
        self.assertNotEqual(hmm1.match_emissions, hmm2.match_emissions)
        # build third HMM with mx=BLOSUM90, popen=0.04, pextend=0.2
        # --> the transition probabilities will change, match emissions won't
        builder.popen = 0.04
        builder.pextend = 0.2
        hmm3, _, _ = builder.build(seq, bg)
        hmm3.command_line = hmm1.creation_time = None
        self.assertNotEqual(hmm2.transition_probabilities, hmm3.transition_probabilities)
        self.assertEqual(hmm2.match_emissions, hmm3.match_emissions)


class TestBuilderMSA(_TestBuilderBase, unittest.TestCase):

    @unittest.skipUnless(os.path.exists(HMMER_FOLDER), "test data not available")
    def test_command_line(self):
        abc = Alphabet.dna()
        builder = Builder(alphabet=abc)

        with MSAFile(os.path.join(self.testdata, "3box.sto"), digital=True) as msa_file:
            msa = next(msa_file)
            msa.name = b"3box"

        with HMMFile(os.path.join(self.testdata, "3box.hmm")) as hmm_file:
            hmm_exp = next(hmm_file)

        argv = ["mycommand", "--param", "1"]
        with mock.patch.object(sys, "argv", argv):
            hmm, profile, opt = builder.build_msa(msa, Background(abc))
        self.assertEqual(hmm.command_line, " ".join(argv))

    @unittest.skipUnless(os.path.exists(HMMER_FOLDER), "test data not available")
    def test_creation_time(self):
        abc = Alphabet.dna()
        builder = Builder(alphabet=abc)

        hmm, profile, opt = builder.build_msa(self.msa, Background(abc))

        with HMMFile(os.path.join(self.testdata, "3box.hmm")) as hmm_file:
            hmm_exp = next(hmm_file)

        self.assertIsInstance(hmm.creation_time, datetime.datetime)

    @unittest.skipUnless(os.path.exists(HMMER_FOLDER), "test data not available")
    def test_dna(self):
        abc = Alphabet.dna()
        builder = Builder(alphabet=abc)

        hmm, profile, opt = builder.build_msa(self.msa, Background(abc))

        with HMMFile(os.path.join(self.testdata, "3box.hmm")) as hmm_file:
            hmm_exp = next(hmm_file)

        self.assertEqual(hmm.name, b"3box")
        self.assertEqual(hmm.M, hmm_exp.M)

    @unittest.skipUnless(os.path.exists(HMMER_FOLDER), "test data not available")
    def test_alphabet_mismatch(self):
        dna = Alphabet.dna()
        amino = Alphabet.amino()

        builder = Builder(alphabet=dna)
        bg = Background(alphabet=amino)
        self.assertRaises(AlphabetMismatch, builder.build_msa, self.msa, bg)

        builder = Builder(alphabet=amino)
        bg = Background(alphabet=amino)
        self.assertRaises(AlphabetMismatch, builder.build_msa, self.msa, bg)

    def test_laccase(self):
        # create a new protein builder
        abc = Alphabet.amino()
        builder = Builder(alphabet=abc)
        # open the MSA
        msa_path = pkg_resources.resource_filename("pyhmmer.tests", "data/msa/laccase.clw")
        with MSAFile(msa_path, digital=True, alphabet=abc) as msa_file:
            msa = msa_file.read()
            msa.name = b"laccase"
        # read the expected HMM
        hmm_path = pkg_resources.resource_filename("pyhmmer.tests", "data/hmms/txt/laccase.hmm")
        with HMMFile(hmm_path) as hmm_file:
            hmm_hmmbuild = next(hmm_file)
            hmm_hmmbuild.creation_time = None
        # build the HMM from the MSA
        hmm_pyhmmer, _, _ = builder.build_msa(msa, Background(abc))
        hmm_pyhmmer.command_line = hmm_pyhmmer.creation_time = None
        # check M and the transition probabilities (which were the issue in v0.4.9)
        # NOTE: exact equality cannot be checked because the `hmmbuild` HMM was
        #       serialized, which drops a bit of numerical accuracy.
        self.assertEqual(hmm_pyhmmer.M, hmm_hmmbuild.M)
        for i in range(hmm_pyhmmer.M):
            for t1, t2 in zip(hmm_pyhmmer.transition_probabilities[i], hmm_hmmbuild.transition_probabilities[i]):
                self.assertAlmostEqual(t1, t2, places=5)
