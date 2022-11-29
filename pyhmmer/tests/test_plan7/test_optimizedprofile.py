import copy
import io
import itertools
import os
import shutil
import unittest
import tempfile
import pkg_resources

import pyhmmer
from pyhmmer.errors import EaselError, AlphabetMismatch
from pyhmmer.easel import Alphabet, SequenceFile, Randomness
from pyhmmer.plan7 import HMM, HMMFile, Profile, Background


class TestOptimizedProfile(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.rng = Randomness(seed=0)
        cls.alphabet = Alphabet.amino()
        cls.background = Background(cls.alphabet)

    @classmethod
    def _random_hmm(cls, name, M=100):
        hmm = HMM.sample(M, cls.alphabet, cls.rng)
        hmm.name = name
        return hmm

    @classmethod
    def _random_profile(cls, name, M=100):
        hmm = cls._random_hmm(name, M=M)
        profile = Profile(hmm.M, hmm.alphabet)
        profile.configure(hmm, cls.background, 200)
        return profile

    @classmethod
    def _random_optimized_profile(cls, name, M=100):
        return cls._random_profile(name, M=M).optimized()

    def test_eq(self):
        profile1 = self._random_profile(b"profile1")
        om1 = profile1.optimized()
        om2 = profile1.optimized()
        self.assertEqual(om1, om2)
        profile2 = self._random_profile(b"profile2")
        om3 = profile2.optimized()
        self.assertNotEqual(om1, om3)

    def test_copy(self):
        om1 = self._random_optimized_profile(b"profile1")
        om2 = copy.copy(om1)
        self.assertEqual(om1, om2)
        self.assertIsNot(om1, om2)

    def test_M(self):
        om = self._random_optimized_profile(b"profile1", 200)
        self.assertEqual(om.M, 200)

    def test_L(self):
        hmm = self._random_hmm(b"profile1")
        profile = Profile(hmm.M, hmm.alphabet)
        profile.configure(hmm, self.background, 200)
        om = profile.optimized()
        self.assertEqual(om.L, 200)
        om.L = 300
        self.assertEqual(om.L, 300)

    def test_accession(self):
        hmm = self._random_hmm(b"profile1")
        hmm.accession = b"TST001"
        profile = Profile(hmm.M, hmm.alphabet)
        profile.configure(hmm, self.background, 200)
        om = profile.optimized()
        self.assertEqual(om.accession, hmm.accession)

    def test_description(self):
        hmm = self._random_hmm(b"profile1")
        hmm.description = b"test profile one"
        profile = Profile(hmm.M, hmm.alphabet)
        profile.configure(hmm, self.background, 200)
        om = profile.optimized()
        self.assertEqual(om.description, hmm.description)

    def test_name(self):
        hmm = self._random_hmm(b"profile1")
        profile = Profile(hmm.M, hmm.alphabet)
        profile.configure(hmm, self.background, 200)
        om = profile.optimized()
        self.assertEqual(om.name, hmm.name)

    def test_consensus(self):
        hmm = self._random_hmm(b"profile1")
        profile = Profile(hmm.M, hmm.alphabet)
        profile.configure(hmm, self.background, 200)
        om1 = profile.optimized()
        self.assertEqual(om1.consensus, hmm.consensus)
        profile = Profile(hmm.M, hmm.alphabet)
        om2 = profile.optimized()
        self.assertIs(om2.consensus, None)

    def test_consensus_structure(self):
        profile = Profile(100, self.alphabet)
        om = profile.optimized()
        self.assertIs(om.consensus_structure, None)

    def test_offsets(self):
        om = self._random_optimized_profile(b"profile1")
        self.assertIs(om.offsets.model, None)
        self.assertIs(om.offsets.profile, None)
        self.assertIs(om.offsets.filter, None)
