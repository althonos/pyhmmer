import copy
import io
import itertools
import os
import pickle
import shutil
import unittest
import tempfile
import pkg_resources

import pyhmmer
from pyhmmer.errors import EaselError
from pyhmmer.easel import Alphabet, SequenceFile, VectorF
from pyhmmer.plan7 import HMM, HMMFile, Pipeline


class TestHMM(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.hmm_path = pkg_resources.resource_filename("pyhmmer.tests", "data/hmms/txt/Thioesterase.hmm")
        with HMMFile(cls.hmm_path) as hmm_file:
            cls.hmm = next(hmm_file)

    def test_checksum(self):
        abc = Alphabet.dna()
        hmm = HMM(100, abc)
        self.assertIs(hmm.checksum, None)

        hmm_path = pkg_resources.resource_filename("pyhmmer.tests", "data/hmms/txt/LuxC.hmm")
        with HMMFile(hmm_path) as hmm_file:
            hmm = next(hmm_file)
        self.assertEqual(hmm.checksum, 395769323)

    def test_composition(self):
        abc = Alphabet.dna()
        hmm = HMM(100, abc)
        self.assertIs(hmm.composition, None)
        hmm.set_composition()
        self.assertEqual(len(hmm.composition), hmm.alphabet.K)
        self.assertEqual(hmm.composition.shape, (hmm.alphabet.K,))

    def test_copy(self):
        hmm = copy.copy(self.hmm)
        self.assertEqual(self.hmm.name, hmm.name)
        self.assertEqual(self.hmm.accession, hmm.accession)
        self.assertEqual(self.hmm.description, hmm.description)
        self.assertEqual(self.hmm.consensus, hmm.consensus)

    def test_eq(self):
        abc = Alphabet.amino()
        self.assertEqual(self.hmm, self.hmm)
        self.assertEqual(self.hmm, self.hmm.copy())
        self.assertNotEqual(self.hmm, HMM(100, abc))
        self.assertNotEqual(self.hmm, 1)

    def test_insert_emissions(self):
        emissions = self.hmm.insert_emissions
        self.assertEqual(emissions.shape, (self.hmm.M+1, self.hmm.alphabet.K))

        row = emissions[0]
        self.assertAlmostEqual(row.sum(), 1.0, places=5)

    def test_command_line(self):
        abc = Alphabet.amino()
        hmm = HMM(100, abc)
        self.assertIs(hmm.command_line, None)
        hmm.command_line = "hmmbuild x"
        self.assertEqual(hmm.command_line, "hmmbuild x")
        hmm.command_line = ""
        self.assertEqual(hmm.command_line, "")
        hmm.command_line = None
        self.assertIs(hmm.command_line, None)

    def test_name(self):
        abc = Alphabet.amino()
        hmm = HMM(100, abc)

        self.assertIs(hmm.name, None)
        hmm.name = b"Test"
        self.assertEqual(hmm.name, b"Test")
        hmm.name = b""
        self.assertEqual(hmm.name, b"")
        hmm.name = None
        self.assertIs(hmm.name, None)

    def test_accession(self):
        abc = Alphabet.amino()
        hmm = HMM(100, abc)

        self.assertIs(hmm.accession, None)
        hmm.accession = b"TST001"
        self.assertEqual(hmm.accession, b"TST001")
        hmm.accession = b""
        self.assertEqual(hmm.accession, b"")
        hmm.accession = None
        self.assertIs(hmm.accession, None)

    def test_description(self):
        abc = Alphabet.amino()
        hmm = HMM(100, abc)

        self.assertIs(hmm.description, None)
        hmm.description = b"A very important HMM"
        self.assertEqual(hmm.description, b"A very important HMM")
        hmm.description = b""
        self.assertEqual(hmm.description, b"")
        hmm.description = None
        self.assertIs(hmm.description, None)

    def test_consensus(self):
        # if not set, HMM is None
        abc = Alphabet.amino()
        hmm = HMM(100, abc)
        self.assertIs(hmm.consensus, None)
        # if the HMM is fully configured, the consensus should be as
        # long as the number of nodes
        self.assertEqual(len(self.hmm.consensus), self.hmm.M)

    def test_write(self):
        buffer = io.BytesIO()
        self.hmm.write(buffer)

        self.assertNotEqual(buffer.tell(), 0)
        buffer.seek(0)

        with open(self.hmm_path, "rb") as f:
            # 1st line is skipped, cause it contains the format date, which will
            # obviously not be the same as the reference
            for line_written, line_ref in itertools.islice(zip(buffer, f), 1, None):
                self.assertEqual(line_written, line_ref)

    def test_write_empty(self):
        buffer = io.BytesIO()
        hmm = HMM(10, Alphabet.amino())
        hmm.write(buffer)
        self.assertGreater(len(buffer.getvalue()), 0)

    def test_write_error(self):
        buffer = io.StringIO()

        with self.assertRaises(EaselError) as ctx:
            self.hmm.write(buffer)
        self.assertEqual(ctx.exception.code, 27)
        self.assertIsInstance(ctx.exception.__cause__, TypeError)

    def test_renormalize(self):
        dna = Alphabet.dna()
        hmm = HMM(10, dna)
        for i in range(1, hmm.M+1):
            for j in range(dna.K):
                hmm.match_emissions[i, j] = 25.0

        hmm.renormalize()
        for i in range(hmm.M+1):
            self.assertAlmostEqual(hmm.match_emissions[i].sum(), 1.0, places=5)

    def test_scale(self):
        dna = Alphabet.dna()
        hmm = HMM(10, dna)
        for i in range(1, hmm.M+1):
            for j in range(dna.K):
                hmm.match_emissions[i, j] = 25.0

        hmm.scale(2.0)
        for i in range(1, hmm.M+1):
            for j in range(dna.K):
                self.assertAlmostEqual(hmm.match_emissions[i, j], 50.0, places=5)

    def test_zero(self):
        dna = Alphabet.dna()
        hmm = HMM(10, dna)
        hmm.name = b"custom"
        for i in range(1, hmm.M+1):
            for j in range(dna.K):
                hmm.match_emissions[i, j] = 25.0

        hmm.zero()
        for i in range(hmm.M+1):
            for j in range(dna.K):
                self.assertEqual(hmm.match_emissions[i, j], 0.0)
        self.assertEqual(hmm.name, b"custom")

    def test_pickle(self):
        h1 = self.hmm
        h2 = pickle.loads(pickle.dumps(h1))

        self.assertEqual(h1.alphabet, h2.alphabet)
        self.assertEqual(h1.M, h2.M)
        self.assertEqual(h1.transition_probabilities, h2.transition_probabilities)
        self.assertEqual(h1.insert_emissions, h2.insert_emissions)
        self.assertEqual(h1.match_emissions, h2.match_emissions)
        self.assertEqual(h1.name, h2.name)
        self.assertEqual(h1.accession, h2.accession)
        self.assertEqual(h1.description, h2.description)
        self.assertEqual(h1.checksum, h2.checksum)
        self.assertEqual(h1.composition, h2.composition)
        self.assertEqual(h1.consensus, h2.consensus)
        self.assertEqual(h1.consensus_structure, h2.consensus_structure)
        self.assertEqual(h1.consensus_accessibility, h2.consensus_accessibility)
        self.assertEqual(h1.reference, h2.reference)
        self.assertEqual(h1.model_mask, h2.model_mask)
        self.assertEqual(h1.command_line, h2.command_line)
        self.assertEqual(h1.nseq, h2.nseq)
        self.assertEqual(h1.nseq_effective, h2.nseq_effective)

        self.assertEqual(h1.creation_time, h2.creation_time)
        self.assertEqual(h1.cutoffs, h2.cutoffs)
        self.assertEqual(h1.evalue_parameters, h2.evalue_parameters)
        self.assertEqual(h1.composition, h2.composition)

        self.assertEqual(h1, h2)

    def test_no_cutoffs(self):
        hmm_path = pkg_resources.resource_filename("pyhmmer.tests", "data/hmms/txt/Thioesterase.hmm")
        with HMMFile(hmm_path) as hmm_file:
            hmm = next(hmm_file)

        self.assertFalse(hmm.cutoffs.gathering_available())
        self.assertIs(hmm.cutoffs.gathering, None)
        self.assertIs(hmm.cutoffs.gathering1, None)
        self.assertIs(hmm.cutoffs.gathering2, None)
        self.assertFalse(hmm.cutoffs.noise_available())
        self.assertIs(hmm.cutoffs.noise, None)
        self.assertIs(hmm.cutoffs.noise1, None)
        self.assertIs(hmm.cutoffs.noise2, None)
        self.assertFalse(hmm.cutoffs.trusted_available())
        self.assertIs(hmm.cutoffs.trusted, None)
        self.assertIs(hmm.cutoffs.trusted1, None)
        self.assertIs(hmm.cutoffs.trusted2, None)

        hmm.cutoffs.gathering = (10.0, 12.0)
        self.assertTrue(hmm.cutoffs.gathering_available())
        self.assertFalse(hmm.cutoffs.noise_available())
        self.assertFalse(hmm.cutoffs.trusted_available())
        self.assertEqual(hmm.cutoffs.gathering, (10.0, 12.0))
        self.assertEqual(hmm.cutoffs.gathering1, 10.0)
        self.assertEqual(hmm.cutoffs.gathering2, 12.0)

        hmm.cutoffs.noise = (8.0, 5.0)
        self.assertTrue(hmm.cutoffs.gathering_available())
        self.assertTrue(hmm.cutoffs.noise_available())
        self.assertFalse(hmm.cutoffs.trusted_available())
        self.assertEqual(hmm.cutoffs.noise, (8.0, 5.0))
        self.assertEqual(hmm.cutoffs.noise1, 8.0)
        self.assertEqual(hmm.cutoffs.noise2, 5.0)

        hmm.cutoffs.trusted = (15.0, 14.0)
        self.assertTrue(hmm.cutoffs.gathering_available())
        self.assertTrue(hmm.cutoffs.noise_available())
        self.assertTrue(hmm.cutoffs.trusted_available())
        self.assertEqual(hmm.cutoffs.trusted, (15.0, 14.0))
        self.assertEqual(hmm.cutoffs.trusted1, 15.0)
        self.assertEqual(hmm.cutoffs.trusted2, 14.0)

    def test_cutoffs(self):
        hmm_path = pkg_resources.resource_filename("pyhmmer.tests", "data/hmms/txt/PF02826.hmm")
        with HMMFile(hmm_path) as hmm_file:
            hmm = next(hmm_file)

        self.assertTrue(hmm.cutoffs.gathering_available())
        self.assertTrue(hmm.cutoffs.noise_available())
        self.assertTrue(hmm.cutoffs.trusted_available())

        hmm.cutoffs.gathering = None
        self.assertFalse(hmm.cutoffs.gathering_available())
        self.assertTrue(hmm.cutoffs.noise_available())
        self.assertTrue(hmm.cutoffs.trusted_available())
        self.assertIs(hmm.cutoffs.gathering, None)
        self.assertIs(hmm.cutoffs.gathering1, None)
        self.assertIs(hmm.cutoffs.gathering2, None)

        del hmm.cutoffs.noise
        self.assertFalse(hmm.cutoffs.gathering_available())
        self.assertFalse(hmm.cutoffs.noise_available())
        self.assertTrue(hmm.cutoffs.trusted_available())
        self.assertIs(hmm.cutoffs.noise, None)
        self.assertIs(hmm.cutoffs.noise1, None)
        self.assertIs(hmm.cutoffs.noise2, None)

        hmm.cutoffs.trusted = None
        self.assertFalse(hmm.cutoffs.gathering_available())
        self.assertFalse(hmm.cutoffs.noise_available())
        self.assertFalse(hmm.cutoffs.trusted_available())
        self.assertIs(hmm.cutoffs.trusted, None)
        self.assertIs(hmm.cutoffs.trusted1, None)
        self.assertIs(hmm.cutoffs.trusted2, None)
