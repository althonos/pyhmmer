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


class TestProfile(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.hmm_path = pkg_resources.resource_filename("tests", "data/hmms/txt/Thioesterase.hmm")
        with HMMFile(cls.hmm_path) as hmm_file:
            cls.hmm = next(hmm_file)

    def setUp(self):
        self.alphabet = self.hmm.alphabet
        self.background = Background(self.alphabet)
        self.profile = Profile(self.hmm.M, self.alphabet)
        self.profile.configure(self.hmm, self.background, 200)

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
        for multihit in (False, True):
            for local in (False, True):
                self.profile.configure(
                    self.hmm,
                    self.background,
                    200,
                    multihit=multihit,
                    local=local
                )
                self.assertEqual(self.profile.is_multihit(), multihit)
                self.assertEqual(self.profile.is_local(), local)

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
        self.assertNotEqual(self.profile.M, 0)
        self.profile.clear()
        self.assertEqual(self.profile.M, 0)
