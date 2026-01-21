import platform
import sys
import unittest

import pyhmmer
from pyhmmer.errors import EaselError
from pyhmmer.easel import SequenceFile, TextMSA, DigitalMSA
from pyhmmer.plan7 import HMMFile, Pipeline, TopHits

from .. import __name__ as __package__
from .utils import HMMER_FOLDER, resource_files


@unittest.skipUnless(resource_files, "importlib.resources.files not available")
class TestAlignment(unittest.TestCase):

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
        self.ali = self.hits[0].best_domain.alignment

    @unittest.skipIf(platform.system() == "Windows", "writing to fileobj unsupported on Windows")
    def test_str(self):
        rendered = str(self.ali)
        lines = rendered.splitlines()
        self.assertEqual(len(lines), 5)
        self.assertTrue(lines[1].strip().startswith(self.hits.query.name))
        self.assertTrue(lines[3].strip().startswith(self.hits[0].name))

    @unittest.skipIf(sys.implementation.name == "pypy", "`getsizeof` not supported on PyPY")
    def test_sizeof(self):
        self.assertGreater(sys.getsizeof(self.ali), 0)