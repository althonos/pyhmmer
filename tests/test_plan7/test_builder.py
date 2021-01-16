import io
import itertools
import os
import shutil
import unittest
import tempfile
import pkg_resources

import pyhmmer
from pyhmmer.errors import EaselError
from pyhmmer.easel import Alphabet, SequenceFile
from pyhmmer.plan7 import HMMFile, Pipeline, Builder, Background


class TestBuilder(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.ecori_fa = os.path.realpath(os.path.join(
            __file__,
            os.pardir,
            os.pardir,
            os.pardir,
            "vendor",
            "hmmer",
            "testsuite",
            "ecori.fa"
        ))
        with SequenceFile(cls.ecori_fa) as seqf:
            cls.dna = next(seqf)

        cls.ecori_hmm = cls.ecori_fa.replace(".fa", ".hmm")
        with HMMFile(cls.ecori_hmm) as hmmf:
            cls.ecori = next(hmmf)

        cls.faa_path = pkg_resources.resource_filename("tests", "data/seqs/PKSI.faa")
        with SequenceFile(cls.faa_path) as seqf:
            cls.proteins = list(seqf)

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
