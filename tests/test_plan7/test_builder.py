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


class test_build_msa_laccase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.testdata = os.path.realpath(os.path.join(
            __file__,
            os.pardir,
            os.pardir,
            os.pardir,
            "vendor",
            "hmmer",
            "testsuite",
        ))

        cls.ecori_fa = os.path.join(cls.testdata, "ecori.fa")
        with SequenceFile(cls.ecori_fa) as seqf:
            cls.dna = next(seqf)

        cls.ecori_hmm = cls.ecori_fa.replace(".fa", ".hmm")
        with HMMFile(cls.ecori_hmm) as hmmf:
            cls.ecori = next(hmmf)

        cls.faa_path = pkg_resources.resource_filename("tests", "data/seqs/PKSI.faa")
        with SequenceFile(cls.faa_path) as seqf:
            cls.proteins = list(seqf)

        cls.msa_path = os.path.join(cls.testdata, "3box.sto")
        with MSAFile(cls.msa_path) as msa_file:
            msa_file.set_digital(msa_file.guess_alphabet())
            cls.msa = next(msa_file)
            cls.msa.name = b"3box"

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

    def test_build_alphabet_mismatch(self):
        abc = Alphabet.dna()
        amino = Alphabet.amino()

        builder = Builder(alphabet=abc)
        bg = Background(alphabet=abc)

        seq = self.proteins[0].digitize(amino)
        self.assertRaises(AlphabetMismatch, builder.build, seq, bg)

    def test_build_command_line(self):
        abc = Alphabet.amino()
        builder = Builder(alphabet=abc)
        seq = self.proteins[1].digitize(abc)

        argv = ["mycommand", "--param", "1"]
        with mock.patch.object(sys, "argv", argv):
            hmm, profile, opt = builder.build(seq, Background(abc))
        self.assertEqual(hmm.command_line, " ".join(argv))

    def test_build_creation_time(self):
        abc = Alphabet.amino()
        builder = Builder(alphabet=abc)
        seq = self.proteins[1].digitize(abc)

        hmm, profile, opt = builder.build(seq, Background(abc))
        self.assertIsInstance(hmm.creation_time, datetime.datetime)

    def test_build_protein(self):
        abc = Alphabet.amino()
        builder = Builder(alphabet=abc)
        seq = self.proteins[1].digitize(abc)
        hmm, profile, opt = builder.build(seq, Background(abc))
        self.assertEqual(hmm.name, seq.name)
        self.assertEqual(hmm.alphabet, abc)
        self.assertEqual(profile.alphabet, abc)
        self.assertEqual(opt.alphabet, abc)

    def test_build_dna(self):
        abc = Alphabet.dna()
        builder = Builder(alphabet=abc)
        hmm, profile, opt = builder.build(self.dna.digitize(abc), Background(abc))
        self.assertEqual(hmm.name, self.dna.name)
        self.assertEqual(hmm.alphabet, abc)
        self.assertEqual(profile.alphabet, abc)
        self.assertEqual(opt.alphabet, abc)
        self.assertEqual(hmm.M, self.ecori.M)

    def test_build_msa_command_line(self):
        abc = Alphabet.dna()
        builder = Builder(alphabet=abc)

        with MSAFile(os.path.join(self.testdata, "3box.sto")) as msa_file:
            msa_file.set_digital(abc)
            msa = next(msa_file)
            msa.name = b"3box"

        with HMMFile(os.path.join(self.testdata, "3box.hmm")) as hmm_file:
            hmm_exp = next(hmm_file)

        argv = ["mycommand", "--param", "1"]
        with mock.patch.object(sys, "argv", argv):
            hmm, profile, opt = builder.build_msa(msa, Background(abc))
        self.assertEqual(hmm.command_line, " ".join(argv))

    def test_build_msa_creation_time(self):
        abc = Alphabet.dna()
        builder = Builder(alphabet=abc)

        hmm, profile, opt = builder.build_msa(self.msa, Background(abc))

        with HMMFile(os.path.join(self.testdata, "3box.hmm")) as hmm_file:
            hmm_exp = next(hmm_file)

        self.assertIsInstance(hmm.creation_time, datetime.datetime)

    def test_build_msa_dna(self):
        abc = Alphabet.dna()
        builder = Builder(alphabet=abc)

        hmm, profile, opt = builder.build_msa(self.msa, Background(abc))

        with HMMFile(os.path.join(self.testdata, "3box.hmm")) as hmm_file:
            hmm_exp = next(hmm_file)

        self.assertEqual(hmm.name, b"3box")
        self.assertEqual(hmm.M, hmm_exp.M)

    def test_build_msa_alphabet_mismatch(self):
        dna = Alphabet.dna()
        amino = Alphabet.amino()

        builder = Builder(alphabet=dna)
        bg = Background(alphabet=amino)
        self.assertRaises(AlphabetMismatch, builder.build_msa, self.msa, bg)

        builder = Builder(alphabet=amino)
        bg = Background(alphabet=amino)
        self.assertRaises(AlphabetMismatch, builder.build_msa, self.msa, bg)

    def test_build_msa_laccase(self):
        # create a new protein builder
        abc = Alphabet.amino()
        builder = Builder(alphabet=abc)
        # open the MSA
        msa_path = pkg_resources.resource_filename("tests", "data/msa/laccase.clw")
        with MSAFile(msa_path) as msa_file:
            msa_file.set_digital(abc)
            msa = msa_file.read()
            msa.name = b"laccase"
        # read the expected HMM
        hmm_path = pkg_resources.resource_filename("tests", "data/hmms/txt/laccase.hmm")
        with HMMFile(hmm_path) as hmm_file:
            hmm_hmmbuild = next(hmm_file)
            hmm_hmmbuild.creation_time = None
        # build the HMM from the MSA
        hmm_pyhmmer, _, _ = builder.build_msa(msa, Background(abc))
        hmm_pyhmmer.command_line = hmm_pyhmmer.creation_time = None
        # NOTE: since saving a HMM loses some numerical precision, it would
        #       be hard to compare for exact equality between the serialized
        #       HMM obtained from `hmmbuild`, and the one just obtained with
        #       the builder. To avoid this, we perform the comparison after
        #       having serialized the two HMMs
        # write pyHMMER HMM
        buffer = io.BytesIO()
        hmm_pyhmmer.write(buffer)
        # write HMMER HMM
        expected = io.BytesIO()
        hmm_hmmbuild.write(expected)
        # compare serialized HMMs
        self.assertMultiLineEqual(
            buffer.getvalue().decode('utf-8'),
            expected.getvalue().decode('utf-8'),
        )
