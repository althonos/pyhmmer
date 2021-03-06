import copy
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
from pyhmmer.plan7 import HMM, HMMFile, Pipeline


class TestHMM(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.hmm_path = pkg_resources.resource_filename("tests", "data/hmms/txt/Thioesterase.hmm")
        with HMMFile(cls.hmm_path) as hmm_file:
            cls.hmm = next(hmm_file)

    def test_copy(self):
        hmm2 = copy.copy(self.hmm)
        self.assertEqual(self.hmm.name, hmm2.name)
        self.assertEqual(self.hmm.accession, hmm2.accession)
        self.assertEqual(self.hmm.description, hmm2.description)
        self.assertEqual(self.hmm.consensus, hmm2.consensus)

    def test_eq(self):
        abc = Alphabet.amino()
        self.assertEqual(self.hmm, self.hmm)
        self.assertEqual(self.hmm, self.hmm.copy())
        self.assertNotEqual(self.hmm, HMM(100, abc))
        self.assertNotEqual(self.hmm, 1)

    def test_name(self):
        abc = Alphabet.amino()
        hmm = HMM(100, abc)

        self.assertIs(hmm.name, None)
        hmm.name = b"Test"
        self.assertEqual(hmm.name, b"Test")
        hmm.name = b""
        self.assertEqual(hmm.name, b"")
        hmm.name = None
        self.assertEqual(hmm.name, None)

    def test_accession(self):
        abc = Alphabet.amino()
        hmm = HMM(100, abc)

        self.assertIs(hmm.accession, None)
        hmm.accession = b"TST001"
        self.assertEqual(hmm.accession, b"TST001")
        hmm.accession = b""
        self.assertEqual(hmm.accession, b"")
        hmm.accession = None
        self.assertEqual(hmm.accession, None)

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

    def test_write_error(self):
        buffer = io.StringIO()

        with self.assertRaises(EaselError) as ctx:
            self.hmm.write(buffer)
        self.assertEqual(ctx.exception.code, 27)
        self.assertIsInstance(ctx.exception.__cause__, TypeError)
