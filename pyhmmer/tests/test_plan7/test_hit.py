import io
import itertools
import os
import shutil
import unittest
import tempfile

import pyhmmer
from pyhmmer.errors import EaselError
from pyhmmer.easel import SequenceFile, TextMSA, DigitalMSA
from pyhmmer.plan7 import HMMFile, Pipeline, TopHits

from .. import __name__ as __package__
from .utils import HMMER_FOLDER, resource_files


@unittest.skipUnless(resource_files, "importlib.resources.files not available")
class TestHit(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        hmm_path = resource_files(__package__).joinpath("data", "hmms", "txt", "PF02826.hmm")
        with HMMFile(hmm_path) as hmm_file:
            hmm = hmm_file.read()
        seqs_path = resource_files(__package__).joinpath("data", "seqs", "938293.PRJEB85.HG003687.faa")
        with SequenceFile(seqs_path, digital=True, alphabet=hmm.alphabet) as seqs_file:
            seqs = seqs_file.read_block()
        pipeline = Pipeline(alphabet=hmm.alphabet)
        cls._hits = pipeline.search_hmm(hmm, seqs)
    
    def setUp(self):
        self.hits = self._hits.copy()

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

    def test_length_property(self):
        hits = self.hits.copy()
        hit = hits[-1]

        self.assertEqual(hit.length, 281)
        self.assertEqual(hits[-1].length, 281)
        self.assertEqual(self.hits[-1].length, 281)

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

    def test_manual_inclusion(self):
        hit = self.hits[-1]

        hit.included = hit.reported = hit.dropped = hit.duplicate = False
        self.assertFalse(hit.dropped)
        self.assertFalse(hit.duplicate)
        self.assertFalse(hit.reported)
        self.assertFalse(hit.included)

        hit.included = True
        self.assertFalse(hit.duplicate)
        self.assertFalse(hit.dropped)
        self.assertTrue(hit.included)
        self.assertTrue(hit.reported) # included hits are always reported

    def test_manual_drop(self):
        hit = self.hits[-1]

        hit.dropped = False
        hit.included = hit.reported = True
        self.assertFalse(hit.dropped)
        self.assertTrue(hit.reported)
        self.assertTrue(hit.included)

        hit.dropped = True
        self.assertTrue(hit.dropped)
        self.assertTrue(hit.reported) # dropped hits may be reported
        self.assertFalse(hit.included) # dropped hits are never included
        
        hit.dropped = False
        self.assertFalse(hit.dropped)
        self.assertTrue(hit.reported)
        self.assertFalse(hit.included)

    def test_manual_duplicate(self):
        hit = self.hits[-1]

        hit.duplicate = hit.dropped = False
        hit.included = hit.reported = True
        self.assertFalse(hit.duplicate)
        self.assertTrue(hit.reported)
        self.assertTrue(hit.included)

        hit.duplicate = True
        self.assertTrue(hit.duplicate)
        self.assertFalse(hit.reported) # duplicate hits are never reported
        self.assertFalse(hit.included) # duplicate hits are never included
        
        hit.duplicate = False
        self.assertFalse(hit.duplicate)
        self.assertFalse(hit.reported) # duplicate hits are never reported
        self.assertFalse(hit.included) # duplicate hits are never included
        
