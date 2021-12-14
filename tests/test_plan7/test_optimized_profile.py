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
from pyhmmer.easel import Alphabet, SequenceFile
from pyhmmer.plan7 import HMM, HMMFile, Profile, Background


class TestOptimizedProfile(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.hmm_path = pkg_resources.resource_filename("tests", "data/hmms/db/Thioesterase.hmm")
        with HMMFile(cls.hmm_path) as hmm_file:
            cls.hmm = next(hmm_file)
        cls.alphabet = cls.hmm.alphabet
        cls.background = Background(cls.alphabet)
        cls.profile = Profile(cls.hmm.M, cls.alphabet)
        cls.profile.configure(cls.hmm, cls.background, 200)
        cls.optimized_profile = cls.profile.optimized()

    # def test_configure(self):
    #     profile = Profile(self.hmm.M, Alphabet.dna())
    #     bg = Background(profile.alphabet)
    #     self.assertNotEqual(profile.alphabet, self.hmm.alphabet)
    #     self.assertRaises(AlphabetMismatch, profile.configure, self.hmm, bg, 200)
    #     hmm = HMM(100, profile.alphabet)
    #     self.assertNotEqual(profile.alphabet, self.background.alphabet)
    #     self.assertRaises(AlphabetMismatch, profile.configure, hmm, self.background, 200)
    #
    # def test_eq(self):
    #     profile2 = Profile(self.hmm.M, self.alphabet)
    #     profile2.configure(self.hmm, self.background, 200)
    #     self.assertEqual(self.profile, profile2)
    #     self.assertNotEqual(self.profile, 1)

    # def test_copy(self):
    #     profile2 = copy.copy(self.profile)
    #     self.assertEqual(self.profile, profile2)
    #
    # def test_profile_modes(self):
    #     profile = self.profile.copy()
    #     for multihit in (False, True):
    #         for local in (False, True):
    #             profile.configure(
    #                 self.hmm,
    #                 self.background,
    #                 200,
    #                 multihit=multihit,
    #                 local=local
    #             )
    #             self.assertEqual(profile.is_multihit(), multihit)
    #             self.assertEqual(profile.is_local(), local)

    def test_M(self):
        self.assertEqual(self.optimized_profile.M, self.hmm.M)

    def test_L(self):
        profile = Profile(self.hmm.M, self.alphabet)
        profile.configure(self.hmm, self.background, 200)
        optimized_profile = profile.optimized()
        self.assertEqual(optimized_profile.L, 200)
        optimized_profile.L = 300
        self.assertEqual(optimized_profile.L, 300)

    def test_accession(self):
        self.assertEqual(self.optimized_profile.accession, self.hmm.accession)

    def test_description(self):
        self.assertEqual(self.optimized_profile.description, self.hmm.description)

    def test_name(self):
        self.assertEqual(self.optimized_profile.name, self.hmm.name)

    def test_consensus(self):
        self.assertEqual(self.optimized_profile.consensus, self.hmm.consensus)
        profile = Profile(self.hmm.M, self.alphabet)
        optimized_profile = profile.optimized()
        self.assertIs(optimized_profile.consensus, None)

    def test_consensus_structure(self):
        self.assertEqual(self.optimized_profile.consensus, self.hmm.consensus)
        profile = Profile(self.hmm.M, self.alphabet)
        optimized_profile = profile.optimized()
        self.assertIs(optimized_profile.consensus, None)

    def test_offsets(self):
        self.assertIs(self.optimized_profile.offsets.model, None)
        self.assertIs(self.optimized_profile.offsets.profile, None)
        self.assertIs(self.optimized_profile.offsets.filter, None)
