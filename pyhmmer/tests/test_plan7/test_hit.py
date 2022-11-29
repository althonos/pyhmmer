import io
import itertools
import os
import shutil
import unittest
import tempfile
import pkg_resources

import pyhmmer
from pyhmmer.errors import EaselError
from pyhmmer.easel import SequenceFile, TextMSA, DigitalMSA
from pyhmmer.plan7 import HMMFile, Pipeline, TopHits


class TestTopHits(unittest.TestCase):

    def setUp(self):
        hmm_path = pkg_resources.resource_filename("pyhmmer.tests", "data/hmms/txt/PF02826.hmm")
        seqs_path = pkg_resources.resource_filename("pyhmmer.tests", "data/seqs/938293.PRJEB85.HG003687.faa")

        with HMMFile(hmm_path) as hmm_file:
            self.hmm = hmm_file.read()
        with SequenceFile(seqs_path, digital=True, alphabet=self.hmm.alphabet) as seqs_file:
            self.seqs = seqs_file.read_block()

        self.pipeline = Pipeline(alphabet=self.hmm.alphabet)
        self.hits = self.pipeline.search_hmm(self.hmm, self.seqs)
        self.pipeline.clear()

    def test_name_property(self):
        hits = self.hits.copy()
        hit = hits[-1]

        self.assertEqual(hit.name, b"938293.PRJEB85.HG003687_187")
        self.assertEqual(hits[-1].name, b"938293.PRJEB85.HG003687_187")
        self.assertEqual(self.hits[-1].name, b"938293.PRJEB85.HG003687_187")

        hits[-1].name = b"new name"

        self.assertEqual(hit.name, b"new name")
        self.assertEqual(hits[-1].name, b"new name")
        self.assertEqual(self.hits[-1].name, b"938293.PRJEB85.HG003687_187")

        with self.assertRaises(TypeError):
            hit.name = None

    def test_accession_property(self):
        hits = self.hits.copy()
        hit = hits[-1]

        self.assertIs(hit.accession, None)
        self.assertIs(hits[-1].accession, None)
        self.assertIs(self.hits[-1].accession, None)

        hits[-1].accession = b"NEW"

        self.assertEqual(hit.accession, b"NEW")
        self.assertEqual(hits[-1].accession, b"NEW")
        self.assertIs(self.hits[-1].accession, None)

        hits[-1].accession = None

        self.assertIs(hit.accession, None)
        self.assertIs(hits[-1].accession, None)
        self.assertIs(self.hits[-1].accession, None)

    def test_description_property(self):
        hits = self.hits.copy()
        hit = hits[-1]

        self.assertTrue(hit.description.startswith(b'# 202177 # 203019 #'))
        self.assertTrue(hits[-1].description.startswith(b'# 202177 # 203019 #'))
        self.assertTrue(self.hits[-1].description.startswith(b'# 202177 # 203019 #'))

        hits[-1].description = None

        self.assertIs(hit.description, None)
        self.assertIs(hits[-1].description, None)
        self.assertTrue(self.hits[-1].description.startswith(b'# 202177 # 203019 #'))

        hits[-1].description = b"NEW"

        self.assertEqual(hit.description, b"NEW")
        self.assertEqual(hits[-1].description, b"NEW")
        self.assertTrue(self.hits[-1].description.startswith(b'# 202177 # 203019 #'))
