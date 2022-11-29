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


class TestProfile(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        rng = Randomness(seed=0)
        cls.alphabet = Alphabet.amino()
        cls.hmm = HMM.sample(100, cls.alphabet, rng)
        cls.hmm.name = b"hmm_one"
        cls.hmm.accession = b"HMM1"
        cls.hmm.description = b"first HMM"
        cls.background = Background(cls.alphabet)
        cls.profile = Profile(cls.hmm.M, cls.alphabet)
        cls.profile.configure(cls.hmm, cls.background, 200)

    def test_configure(self):
        profile = Profile(self.hmm.M, Alphabet.dna())
        bg = Background(profile.alphabet)
        self.assertNotEqual(profile.alphabet, self.hmm.alphabet)
        self.assertRaises(AlphabetMismatch, profile.configure, self.hmm, bg, 200)
        hmm = HMM(100, profile.alphabet)
        self.assertNotEqual(profile.alphabet, self.background.alphabet)
        self.assertRaises(AlphabetMismatch, profile.configure, hmm, self.background, 200)

    def test_eq(self):
        profile2 = Profile(self.hmm.M, self.alphabet)
        profile2.configure(self.hmm, self.background, 200)
        self.assertEqual(self.profile, profile2)
        self.assertNotEqual(self.profile, 1)

    def test_copy(self):
        profile2 = copy.copy(self.profile)
        self.assertEqual(self.profile, profile2)

    def test_profile_modes(self):
        profile = self.profile.copy()
        for multihit in (False, True):
            for local in (False, True):
                profile.configure(
                    self.hmm,
                    self.background,
                    200,
                    multihit=multihit,
                    local=local
                )
                self.assertEqual(profile.is_multihit(), multihit)
                self.assertEqual(profile.is_local(), local)

    def test_M(self):
        self.assertEqual(self.profile.M, self.hmm.M)

    def test_L(self):
        profile = Profile(self.hmm.M, self.alphabet)
        profile.configure(self.hmm, self.background, 200)
        self.assertEqual(profile.L, 200)
        profile.L = 300
        self.assertEqual(profile.L, 300)

    def test_accession(self):
        self.assertEqual(self.profile.accession, self.hmm.accession)

    def test_description(self):
        self.assertEqual(self.profile.description, self.hmm.description)

    def test_name(self):
        self.assertEqual(self.profile.name, self.hmm.name)

    def test_clear(self):
        profile = self.profile.copy()
        self.assertNotEqual(profile.M, 0)
        profile.clear()
        self.assertEqual(profile.M, 0)

    def test_consensus(self):
        self.assertEqual(self.profile.consensus, self.hmm.consensus)
        profile = Profile(self.hmm.M, self.alphabet)
        self.assertIs(profile.consensus, None)

    def test_consensus_structure(self):
        self.assertEqual(self.profile.consensus, self.hmm.consensus)
        profile = Profile(self.hmm.M, self.alphabet)
        self.assertIs(profile.consensus, None)

    def test_offsets(self):
        self.assertIs(self.profile.offsets.model, None)
        self.assertIs(self.profile.offsets.profile, None)
        self.assertIs(self.profile.offsets.filter, None)
