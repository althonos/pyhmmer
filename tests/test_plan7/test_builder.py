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
from pyhmmer.errors import EaselError
from pyhmmer.easel import Alphabet, SequenceFile, MSAFile
from pyhmmer.plan7 import HMMFile, Pipeline, Builder, Background


class TestBuilder(unittest.TestCase):

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

    def test_copy(self):
        abc = Alphabet.amino()
        builder = Builder(alphabet=abc)
        copy = builder.copy()
        self.assertEqual(builder.alphabet, copy.alphabet)
        self.assertEqual(builder.seed, copy.seed)

    def test_build_command_line(self):
        abc = Alphabet.amino()
        builder = Builder(alphabet=abc)
        seq = self.proteins[1].digitize(abc)

        argv = ["mycommand", "--param", "1"]
        with mock.patch.object(sys, "argv", argv):
            hmm, profile, opt = builder.build(seq, Background(abc))
        self.assertEqual(hmm.command_line, " ".join(argv))

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

    def test_build_msa_dna(self):
        abc = Alphabet.dna()
        builder = Builder(alphabet=abc)

        with MSAFile(os.path.join(self.testdata, "3box.sto")) as msa_file:
            msa_file.set_digital(abc)
            msa = next(msa_file)
            msa.name = b"3box"

        with HMMFile(os.path.join(self.testdata, "3box.hmm")) as hmm_file:
            hmm_exp = next(hmm_file)

        hmm, profile, opt = builder.build_msa(msa, Background(abc))
        self.assertEqual(hmm.name, b"3box")
        self.assertEqual(hmm.M, hmm_exp.M)
