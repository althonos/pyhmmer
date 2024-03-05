import abc
import collections
import math
import io
import itertools
import os
import unittest
import tempfile
import threading

import pyhmmer
from pyhmmer.plan7 import Pipeline, HMMFile, HMMPressedFile, TopHits, Hit
from pyhmmer.easel import Alphabet, DigitalMSA, MSAFile, SequenceFile, TextSequence

from .utils import resource_files

class _TestSearch(metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def get_hits(self, hmm, sequences):
        return NotImplemented

    @staticmethod
    def table(name):
        path = resource_files(__package__).joinpath("data", "tables", name)
        return path.open()

    @staticmethod
    def hmm_file(name):
        path = resource_files(__package__).joinpath("data", "hmms", "txt", "{}.hmm".format(name))
        return HMMFile(path)

    @staticmethod
    def seqs_file(name, digital=False):
        seqs_path = resource_files(__package__).joinpath("data", "seqs", "{}.faa".format(name))
        return SequenceFile(seqs_path, digital=digital)

    @unittest.skipUnless(resource_files, "importlib.resources not available")
    def test_thioestherase(self):
        # $ hmmsearch data/hmm/Thioesterase.hmm data/seqs/938293.PRJEB85.HG003687.faa
        #
        # Query:       Thioesterase  [M=243]
        # Scores for complete sequences (score includes all domains):
        #    --- full sequence ---   --- best 1 domain ---    -#dom-
        #     E-value  score  bias    E-value  score  bias    exp  N  Sequence                    Description
        #     ------- ------ -----    ------- ------ -----   ---- --  --------                    -----------
        #   ------ inclusion threshold ------
        #       0.096    8.6   1.5       0.14    8.1   1.5    1.2  1  938293.PRJEB85.HG003687_113  # 108502 # 109284 # 1 # ID=8_113;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.341
        #
        # Domain annotation for each sequence (and alignments):
        # >> 938293.PRJEB85.HG003687_113  # 108502 # 109284 # 1 # ID=8_113;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.341
        #    #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
        #  ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
        #    1 ?    8.1   1.5   6.5e-05      0.14      79      93 ..     115     129 ..     115     129 .. 0.96
        #
        #   Alignments for each domain:
        #   == domain 1  score: 8.1 bits;  conditional E-value: 6.5e-05
        #                  Thioesterase  79 GWSfGGvlAyEmArq 93
        #                                   G+S+GG +A ++A++
        #   938293.PRJEB85.HG003687_113 115 GHSMGGSVAVAIAHE 129
        #                                   9************96 PP
        with self.hmm_file("Thioesterase") as hmm_file:
            hmm = hmm_file.read()
        with self.seqs_file("938293.PRJEB85.HG003687", digital=True) as seqs_file:
            seqs = seqs_file.read_block()

        hits = self.get_hits(hmm, seqs)
        self.assertEqual(len(hits), 1)

        hits.sort()

        hit = hits[0]
        self.assertEqual(hit.name, b"938293.PRJEB85.HG003687_113")
        self.assertAlmostEqual(hit.score, 8.6, delta=0.1)  # printed with %6.1f
        self.assertAlmostEqual(hit.bias, 1.5, delta=0.1)  # printed with  %5.1f
        self.assertAlmostEqual(hit.evalue, 0.096, delta=0.01)  # printed with %9.2g
        self.assertEqual(len(hit.domains), 1)

        domain = hit.domains[0]
        self.assertAlmostEqual(domain.score, 8.1, delta=0.1)
        self.assertAlmostEqual(domain.bias, 1.5, delta=0.1)
        self.assertAlmostEqual(domain.i_evalue, 0.14, places=2)  # printed with %9.2g
        self.assertAlmostEqual(domain.c_evalue, 6.5e-05, places=2)  # printed with %9.2g
        self.assertEqual(domain.alignment.target_from, 115)
        self.assertEqual(domain.alignment.target_to, 129)
        self.assertEqual(domain.alignment.target_length, 261)
        self.assertEqual(domain.alignment.hmm_from, 79)
        self.assertEqual(domain.alignment.hmm_to, 93)
        self.assertEqual(domain.alignment.hmm_length, 243)
        self.assertEqual(domain.env_from, 115)
        self.assertEqual(domain.env_to, 129)

        last = hit.domains[-1]
        self.assertEqual(domain.alignment.hmm_name, last.alignment.hmm_name)

    @unittest.skipUnless(resource_files, "importlib.resources not available")
    def test_pf02826_file(self):
        with self.hmm_file("PF02826") as hmm_file:
            hmm = next(hmm_file)
        with self.seqs_file("938293.PRJEB85.HG003687", digital=True) as seqs_file:
            hits = self.get_hits(hmm, seqs_file)
            self.assertEqual(len(hits), 22)

        hits.sort()

        with self.table("PF02826.tbl") as table:
            lines = filter(lambda line: not line.startswith("#"), table)
            for line, hit in itertools.zip_longest(lines, hits):
                fields = list(filter(None, line.strip().split(" ")))
                self.assertIsNot(line, None)
                self.assertIsNot(hit, None)
                self.assertEqual(hit.name.decode(), fields[0])
                if fields[1] == "-":
                    self.assertIs(hit.accession, None)
                else:
                    self.assertEqual(hit.accession.decode(), fields[1])
                self.assertAlmostEqual(hit.score, float(fields[5]), delta=0.1)
                self.assertAlmostEqual(hit.bias, float(fields[6]), delta=0.1)
                self.assertAlmostEqual(hit.evalue, float(fields[4]), delta=0.1)

    @unittest.skipUnless(resource_files, "importlib.resources not available")
    def test_pf02826_block(self):
        with self.hmm_file("PF02826") as hmm_file:
            hmm = next(hmm_file)
        with self.seqs_file("938293.PRJEB85.HG003687", digital=True) as seqs_file:
            seqs = seqs_file.read_block()

        hits = self.get_hits(hmm, seqs)
        self.assertEqual(len(hits), 22)

        hits.sort()

        with self.table("PF02826.tbl") as table:
            lines = filter(lambda line: not line.startswith("#"), table)
            for line, hit in itertools.zip_longest(lines, hits):
                fields = list(filter(None, line.strip().split(" ")))
                self.assertIsNot(line, None)
                self.assertIsNot(hit, None)
                self.assertEqual(hit.name.decode(), fields[0])
                if fields[1] == "-":
                    self.assertIs(hit.accession, None)
                else:
                    self.assertEqual(hit.accession.decode(), fields[1])
                self.assertAlmostEqual(hit.score, float(fields[5]), delta=0.1)
                self.assertAlmostEqual(hit.bias, float(fields[6]), delta=0.1)
                self.assertAlmostEqual(hit.evalue, float(fields[4]), delta=0.1)

    @unittest.skipUnless(resource_files, "importlib.resources not available")
    def test_t2pks(self):
        with self.seqs_file("938293.PRJEB85.HG003687", digital=True) as seqs_file:
            seqs = seqs_file.read_block()

        with self.hmm_file("t2pks") as hmm_file:
            all_hits = [self.get_hits(hmm, seqs) for hmm in hmm_file]

        hits = (hit for hits in all_hits for hit in hits)
        with self.table("t2pks.tbl") as table:
            lines = filter(lambda line: not line.startswith("#"), table)
            for line, hit in itertools.zip_longest(lines, hits):
                fields = list(filter(None, line.strip().split(" ")))
                self.assertIsNot(line, None)
                self.assertIsNot(hit, None)
                self.assertEqual(hit.name.decode(), fields[0])
                if fields[1] == "-":
                    self.assertIs(hit.accession, None)
                else:
                    self.assertEqual(hit.accession.decode(), fields[1])
                self.assertAlmostEqual(hit.score, float(fields[5]), delta=0.1)
                self.assertAlmostEqual(hit.bias, float(fields[6]), delta=0.1)
                self.assertAlmostEqual(hit.evalue, float(fields[4]), delta=0.1)

        domains = [
            domain for hits in all_hits for hit in hits for domain in hit.domains
        ]
        with self.table("t2pks.domtbl") as table:
            lines = filter(lambda line: not line.startswith("#"), table)
            for line, domain in itertools.zip_longest(lines, domains):
                fields = list(filter(None, line.strip().split(" ")))
                self.assertIsNot(line, None)
                self.assertIsNot(domain, None)
                self.assertEqual(domain.hit.hits.Z, len(seqs))
                self.assertEqual(domain.hit.name.decode(), fields[0])
                self.assertAlmostEqual(domain.score, float(fields[13]), places=1)
                self.assertAlmostEqual(domain.bias, float(fields[14]), places=1)
                # FIXME: it looks like domZ is not extracted properly
                # self.assertEqual(f"{domain.c_evalue:9.2g}", f"{float(fields[11]):9.2g}")
                self.assertEqual(f"{domain.i_evalue:9.2g}", f"{float(fields[12]):9.2g}")

    @unittest.skipUnless(resource_files, "importlib.resources not available")
    def test_hmm_vs_optimized_profile(self):
        with self.hmm_file("PF02826") as hmm_file:
            hmm = next(hmm_file)
            bg = pyhmmer.plan7.Background(hmm.alphabet)
            profile = pyhmmer.plan7.Profile(hmm.M, hmm.alphabet)
            profile.configure(hmm, bg, 100)
        with self.seqs_file("938293.PRJEB85.HG003687", digital=True) as seqs_file:
            seqs = seqs_file.read_block()

        hits_hmm = self.get_hits(hmm, seqs)
        self.assertEqual(len(hits_hmm), 22)

        hits_oprofile = self.get_hits(profile.to_optimized(), seqs)
        self.assertEqual(len(hits_oprofile), 22)

        for hit_hmm, hit_oprofile in itertools.zip_longest(hits_hmm, hits_oprofile):
            self.assertEqual(hit_hmm.name, hit_oprofile.name)
            self.assertEqual(hit_hmm.accession, hit_oprofile.accession)
            self.assertEqual(hit_hmm.description, hit_oprofile.description)
            self.assertEqual(hit_hmm.score, hit_oprofile.score)
            self.assertEqual(hit_hmm.pre_score, hit_oprofile.pre_score)
            self.assertEqual(hit_hmm.sum_score, hit_oprofile.sum_score)
            self.assertEqual(hit_hmm.bias, hit_oprofile.bias)
            self.assertEqual(hit_hmm.evalue, hit_oprofile.evalue)
            self.assertEqual(hit_hmm.pvalue, hit_oprofile.pvalue)
            for dom_hmm, dom_oprofile in itertools.zip_longest(
                hit_hmm.domains, hit_oprofile.domains
            ):
                self.assertEqual(dom_hmm.env_from, dom_oprofile.env_from)
                self.assertEqual(dom_hmm.env_to, dom_oprofile.env_to)
                self.assertEqual(dom_hmm.score, dom_oprofile.score)
                self.assertEqual(dom_hmm.bias, dom_oprofile.bias)
                self.assertEqual(dom_hmm.correction, dom_oprofile.correction)
                self.assertEqual(dom_hmm.envelope_score, dom_oprofile.envelope_score)
                self.assertEqual(dom_hmm.c_evalue, dom_oprofile.c_evalue)
                self.assertEqual(dom_hmm.i_evalue, dom_oprofile.i_evalue)
                self.assertEqual(dom_hmm.pvalue, dom_oprofile.pvalue)


class TestHmmsearch(_TestSearch, unittest.TestCase):
    def get_hits(self, hmm, seqs):
        return list(pyhmmer.hmmsearch(hmm, seqs))[0]

    @unittest.skipUnless(resource_files, "importlib.resources not available")
    def test_callback_error(self):

        class MyException(Exception):
            pass

        def callback(hmm, total):
            raise MyException("oopsie")

        with self.hmm_file("Thioesterase") as hmm_file:
            hmm = hmm_file.read()
        with self.seqs_file("938293.PRJEB85.HG003687", digital=True) as seqs_file:
            seqs = seqs_file.read_block()

        hits = pyhmmer.hmmsearch(hmm, seqs, cpus=1, callback=callback)
        with self.assertRaises(MyException):
            hit = next(hits)

        hits = pyhmmer.hmmsearch(hmm, seqs, cpus=2, callback=callback)
        with self.assertRaises(MyException):
            hit = next(hits)

    @unittest.skipUnless(resource_files, "importlib.resources not available")
    def test_background_error(self):
        # check that errors occuring in worker threads are recovered and raised
        # in the main threads (a common error is mismatching the HMM and the
        # sequence alphabets).
        seqs = [TextSequence().digitize(Alphabet.dna())]
        with self.hmm_file("PF02826") as hmm_file:
            hmm = next(hmm_file)
        self.assertRaises(ValueError, self.get_hits, hmm, seqs)

    @unittest.skipUnless(resource_files, "importlib.resources not available")
    def test_no_queries(self):
        with self.seqs_file("938293.PRJEB85.HG003687", digital=True) as seqs_file:
            seqs = list(seqs_file)
        hits = pyhmmer.hmmsearch([], seqs)
        self.assertIs(None, next(hits, None))


class TestHmmsearchSingle(TestHmmsearch, unittest.TestCase):
    def get_hits(self, hmm, seqs):
        return list(pyhmmer.hmmsearch(hmm, seqs, cpus=1))[0]

    @unittest.skipUnless(resource_files, "importlib.resources not available")
    def test_no_queries(self):
        with self.seqs_file("938293.PRJEB85.HG003687", digital=True) as seqs_file:
            seqs = list(seqs_file)
        hits = pyhmmer.hmmsearch([], seqs)
        self.assertIs(None, next(hits, None))


class TestPipelinesearch(_TestSearch, unittest.TestCase):
    def get_hits(self, hmm, seqs):
        pipeline = Pipeline(alphabet=hmm.alphabet)
        hits = pipeline.search_hmm(hmm, seqs)
        return hits


class TestHmmpress(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.NamedTemporaryFile(suffix=".hmm", delete=False).name

    def tearDown(self):
        for ext in ["", ".h3m", ".h3p", ".h3i", ".h3f"]:
            if os.path.exists(self.tmp + ext):
                os.remove(self.tmp + ext)

    @unittest.skipUnless(resource_files, "importlib.resources not available")
    def test_roundtrip(self):
        db_folder = resource_files(__package__).joinpath("data", "hmms", "db")
        self.hmm = db_folder.joinpath("Thioesterase.hmm")
        self.h3p = db_folder.joinpath("Thioesterase.hmm.h3p")
        self.h3m = db_folder.joinpath("Thioesterase.hmm.h3m")
        self.h3f = db_folder.joinpath("Thioesterase.hmm.h3f")
        self.h3i = db_folder.joinpath("Thioesterase.hmm.h3f")
        with HMMFile(self.hmm) as hmms:
            n = pyhmmer.hmmer.hmmpress(hmms, self.tmp)
            self.assertEqual(n, 1)
        with HMMFile(self.tmp) as hmm_file:
            hmm = next(hmm_file)
            self.assertEqual(hmm.name, b"Thioesterase")


class TestPhmmer(unittest.TestCase):
    @staticmethod
    def table(name):
        path = resource_files(__package__).joinpath("data", "tables", name)
        return path.open()

    @unittest.skipUnless(resource_files, "importlib.resources not available")
    def test_callback_error(self):

        class MyException(Exception):
            pass

        def callback(hmm, total):
            raise MyException("oopsie")

        alphabet = Alphabet.amino()
        path = resource_files(__package__).joinpath("data", "seqs", "PKSI.faa")
        with SequenceFile(path, digital=True, alphabet=alphabet) as seqs_file:
            seqs = seqs_file.read_block()

        hits = pyhmmer.phmmer(seqs[-1:], seqs, cpus=1, callback=callback)
        with self.assertRaises(MyException):
            hit = next(hits)

        hits = pyhmmer.phmmer(seqs[-1:], seqs, cpus=2, callback=callback)
        with self.assertRaises(MyException):
            hit = next(hits)

    @unittest.skipUnless(resource_files, "importlib.resources not available")
    def test_no_queries(self):
        alphabet = Alphabet.amino()
        path = resource_files(__package__).joinpath("data", "seqs", "PKSI.faa")
        with SequenceFile(path, digital=True, alphabet=alphabet) as seqs_file:
            seqs = seqs_file.read_block()
        hits = pyhmmer.phmmer([], seqs, cpus=1)
        self.assertIs(None, next(hits, None))

    @unittest.skipUnless(resource_files, "importlib.resources not available")
    def test_pksi(self):
        alphabet = Alphabet.amino()
        path = resource_files(__package__).joinpath("data", "seqs", "PKSI.faa")
        with SequenceFile(path, digital=True, alphabet=alphabet) as seqs_file:
            seqs = seqs_file.read_block()
        hits = next(pyhmmer.phmmer(seqs[-1:], seqs, cpus=1))
        hits.sort()

        with self.table("A0A089QRB9.domtbl") as table:
            lines = iter(filter(lambda line: not line.startswith("#"), table))
            it = ((hit, domain) for hit in hits for domain in hit.domains)
            for line, (hit, domain) in itertools.zip_longest(lines, it):
                self.assertIsNot(line, None)
                self.assertIsNot(hit, None)
                fields = list(filter(None, line.strip().split(" ")))

                self.assertEqual(hit.name.decode(), fields[0])
                self.assertAlmostEqual(hit.score, float(fields[7]), delta=0.1)
                self.assertAlmostEqual(hit.bias, float(fields[8]), delta=0.1)
                self.assertAlmostEqual(hit.evalue, float(fields[6]), delta=0.1)

                self.assertAlmostEqual(domain.i_evalue, float(fields[12]), delta=0.1)
                self.assertAlmostEqual(domain.score, float(fields[13]), delta=0.1)

                self.assertEqual(domain.alignment.hmm_from, int(fields[15]))
                self.assertEqual(domain.alignment.hmm_to, int(fields[16]))
                self.assertEqual(domain.alignment.target_from, int(fields[17]))
                self.assertEqual(domain.alignment.target_to, int(fields[18]))
                self.assertEqual(domain.env_from, int(fields[19]))
                self.assertEqual(domain.env_to, int(fields[20]))


@unittest.skipUnless(resource_files, "importlib.resources not available")
class TestJackhmmer(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        with cls.seqs_file("PKSI", digital=True) as seqs_file:
            cls.pksi = seqs_file.read_block()

    @staticmethod
    def table(name):
        path = resource_files(__package__).joinpath("data", "tables", name)
        return path.open()

    @staticmethod
    def seqs_file(name, digital=False):
        seqs_path = resource_files(__package__).joinpath("data", "seqs", "{}.faa".format(name))
        return SequenceFile(seqs_path, digital=digital)

    @staticmethod
    def hmm_file(name):
        path = resource_files(__package__).joinpath("data", "hmms", "txt", "{}.hmm".format(name))
        return HMMFile(path)

    def test_callback_error(self):

        class MyException(Exception):
            pass

        def callback(hmm, total):
            raise MyException("oopsie")

        seqs = self.pksi
        results = pyhmmer.jackhmmer(seqs[-1:], seqs, cpus=1, max_iterations=1, callback=callback)
        with self.assertRaises(MyException):
            result = next(results)

        results = pyhmmer.jackhmmer(seqs[-1:], seqs, cpus=2, max_iterations=1, callback=callback)
        with self.assertRaises(MyException):
            result = next(results)

    def test_no_queries(self):
        seqs = self.pksi
        result = pyhmmer.jackhmmer([], seqs, cpus=1, max_iterations=1)
        self.assertIs(None, next(result, None))

    def test_pksi(self):
        seqs = self.pksi
        results = pyhmmer.jackhmmer(seqs[-1:], seqs, cpus=1, max_iterations=1)
        result = next(results)
        self.assertEqual(result.iteration, 1)
        result.hits.sort()

        with self.table("A0A089QRB9.domtbl") as table:
            lines = iter(filter(lambda line: not line.startswith("#"), table))
            it = ((hit, domain) for hit in result.hits for domain in hit.domains)
            for line, (hit, domain) in itertools.zip_longest(lines, it):
                self.assertIsNot(line, None)
                self.assertIsNot(hit, None)
                fields = list(filter(None, line.strip().split(" ")))

                self.assertEqual(hit.name.decode(), fields[0])
                self.assertAlmostEqual(hit.score, float(fields[7]), delta=0.1)
                self.assertAlmostEqual(hit.bias, float(fields[8]), delta=0.1)
                self.assertAlmostEqual(hit.evalue, float(fields[6]), delta=0.1)

                self.assertAlmostEqual(domain.i_evalue, float(fields[12]), delta=0.1)
                self.assertAlmostEqual(domain.score, float(fields[13]), delta=0.1)

                self.assertEqual(domain.alignment.hmm_from, int(fields[15]))
                self.assertEqual(domain.alignment.hmm_to, int(fields[16]))
                self.assertEqual(domain.alignment.target_from, int(fields[17]))
                self.assertEqual(domain.alignment.target_to, int(fields[18]))
                self.assertEqual(domain.env_from, int(fields[19]))
                self.assertEqual(domain.env_to, int(fields[20]))

    def test_pksi_checkpoint(self):
        seqs = self.pksi
        # jackhmmer CLI converges in 3 iterations, 5 hits, 17 sequences in MSA
        iterations = next(
            pyhmmer.jackhmmer(
                seqs[-1],
                seqs,
                cpus=1,
                checkpoints=True,
            )
        )
        self.assertEqual(len(iterations), 3)
        self.assertTrue(iterations[-1].converged)
        self.assertEqual(len(iterations[-1].hits), 5)
        self.assertEqual(len(iterations[-1].msa.sequences), 17)

    def test_thioestherase(self):
        with self.hmm_file("Thioesterase") as hmm_file:
            hmm = hmm_file.read()
        with self.seqs_file("938293.PRJEB85.HG003687", digital=True) as seqs_file:
            seqs = seqs_file.read_block()

        result = next(
            pyhmmer.jackhmmer(
                hmm, seqs, cpus=1, max_iterations=1, incE=0.1, incdomE=0.1
            )
        )
        # unpack IterationResult
        _, hits, _, _, it = result
        self.assertEqual(it, 1)
        self.assertEqual(len(hits), 1)

        hits.sort()

        hit = hits[0]
        self.assertEqual(hit.name, b"938293.PRJEB85.HG003687_113")
        self.assertAlmostEqual(hit.score, 8.6, delta=0.1)  # printed with %6.1f
        self.assertAlmostEqual(hit.bias, 1.5, delta=0.1)  # printed with  %5.1f
        self.assertAlmostEqual(hit.evalue, 0.096, delta=0.01)  # printed with %9.2g
        self.assertEqual(len(hit.domains), 1)

        domain = hit.domains[0]
        self.assertAlmostEqual(domain.score, 8.1, delta=0.1)
        self.assertAlmostEqual(domain.bias, 1.5, delta=0.1)
        self.assertAlmostEqual(domain.i_evalue, 0.14, places=2)  # printed with %9.2g
        self.assertAlmostEqual(domain.c_evalue, 6.5e-05, places=2)  # printed with %9.2g
        self.assertEqual(domain.alignment.target_from, 115)
        self.assertEqual(domain.alignment.target_to, 129)
        self.assertEqual(domain.alignment.hmm_from, 79)
        self.assertEqual(domain.alignment.hmm_to, 93)
        self.assertEqual(domain.env_from, 115)
        self.assertEqual(domain.env_to, 129)

        last = hit.domains[-1]
        self.assertEqual(domain.alignment.hmm_name, last.alignment.hmm_name)


@unittest.skipUnless(resource_files, "importlib.resources not available")
class TestNhmmer(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = resource_files(__package__).joinpath("data", "hmms", "txt", "bmyD.hmm")
        with HMMFile(path) as hmm_file:
            cls.bmyD = next(hmm_file)
        path = resource_files(__package__).joinpath("data", "hmms", "txt", "RF00001.hmm")
        with HMMFile(path) as hmm_file:
            cls.rf00001 = next(hmm_file)

    @staticmethod
    def table(name):
        path = resource_files(__package__).joinpath("data", "tables", name)
        return path.open()

    def assertTableEqual(self, hits, table):
        lines = iter(filter(lambda line: not line.startswith("#"), table))
        reported_hits = filter(lambda hit: hit.reported, hits)
        for line, hit in itertools.zip_longest(lines, reported_hits):
            self.assertIsNot(line, None)
            self.assertIsNot(hit, None)
            fields = list(filter(None, line.strip().split(" ")))
            self.assertEqual(hit.name.decode(), fields[0])
            if fields[1] == "-":
                self.assertIs(hit.accession, None)
            else:
                self.assertEqual(hit.accession.decode(), fields[1])
            self.assertAlmostEqual(hit.best_domain.bias, float(fields[14]), delta=0.1)
            self.assertAlmostEqual(hit.best_domain.score, float(fields[13]), delta=0.1)
            self.assertAlmostEqual(
                hit.best_domain.i_evalue, float(fields[12]), delta=0.1
            )

    def test_no_queries(self):
        alphabet = Alphabet.dna()
        path = resource_files(__package__).joinpath("data", "seqs", "BGC0001090.gbk")
        with SequenceFile(path, digital=True, alphabet=alphabet) as seqs_file:
            seqs = seqs_file.read_block()
        hits = pyhmmer.nhmmer([], seqs, cpus=1)
        self.assertIs(None, next(hits, None))

    def test_bmyd_seq_bgc_block(self):
        alphabet = Alphabet.dna()

        path = resource_files(__package__).joinpath("data", "seqs", "bmyD.fna")
        with SequenceFile(path, digital=True, alphabet=alphabet) as seqs_file:
            query = next(seqs_file)

        path = resource_files(__package__).joinpath("data", "seqs", "BGC0001090.gbk")
        with SequenceFile(
            path, "genbank", digital=True, alphabet=alphabet
        ) as seqs_file:
            seqs = seqs_file.read_block()

        hits = next(pyhmmer.nhmmer(query, seqs, cpus=1))
        hits.sort()

        self.assertEqual(len(hits), 1)
        with self.table("bmyD3.tbl") as table:
            self.assertTableEqual(hits, table)

    def test_bmyd_seq_bgc_file(self):
        alphabet = Alphabet.dna()

        path = resource_files(__package__).joinpath("data", "seqs", "bmyD.fna")
        with SequenceFile(path, digital=True, alphabet=alphabet) as seqs_file:
            query = next(seqs_file)

        path = resource_files(__package__).joinpath("data", "seqs", "BGC0001090.gbk")
        with SequenceFile(path, "genbank", digital=True, alphabet=alphabet) as seqs:
            hits = list(pyhmmer.nhmmer(query, seqs, cpus=1))[0]
            hits.sort()

        self.assertEqual(len(hits.reported), 1)
        with self.table("bmyD3.tbl") as table:
            self.assertTableEqual(hits, table)

    def test_bmyd_msa_bgc_block(self):
        alphabet = Alphabet.dna()

        path = resource_files(__package__).joinpath("data", "seqs", "bmyD.fna")
        with SequenceFile(path, digital=True, alphabet=alphabet) as seqs_file:
            query = DigitalMSA(alphabet, name=b"bmyD", sequences=[next(seqs_file)])

        path = resource_files(__package__).joinpath("data", "seqs", "BGC0001090.gbk")
        with SequenceFile(
            path, "genbank", digital=True, alphabet=alphabet
        ) as seqs_file:
            seqs = seqs_file.read_block()

        hits = next(pyhmmer.nhmmer(query, seqs, cpus=1))
        self.assertEqual(len(hits), 1)

    def test_bmyd_msa_bgc_file(self):
        alphabet = Alphabet.dna()

        path = resource_files(__package__).joinpath("data", "seqs", "bmyD.fna")
        with SequenceFile(path, digital=True, alphabet=alphabet) as seqs_file:
            query = DigitalMSA(alphabet, name=b"bmyD", sequences=[next(seqs_file)])

        path = resource_files(__package__).joinpath("data", "seqs", "BGC0001090.gbk")
        with SequenceFile(path, "genbank", digital=True, alphabet=alphabet) as seqs:
            hits = list(pyhmmer.nhmmer(query, seqs, cpus=1))[0]
            self.assertEqual(len(hits.reported), 1)

    def test_bmyd_hmm_bgc_block(self):
        alphabet = Alphabet.dna()

        path = resource_files(__package__).joinpath("data", "seqs", "BGC0001090.gbk")
        with SequenceFile(
            path, "genbank", digital=True, alphabet=alphabet
        ) as seqs_file:
            seqs = seqs_file.read_block()

        hits = next(pyhmmer.nhmmer(self.bmyD, seqs, cpus=1))
        hits.sort()

        self.assertEqual(len(hits.reported), 2)
        with self.table("bmyD1.tbl") as table:
            self.assertTableEqual(hits, table)

    def test_bmyd_hmm_bgc_file(self):
        alphabet = Alphabet.dna()

        path = resource_files(__package__).joinpath("data", "seqs", "BGC0001090.gbk")
        with SequenceFile(
            path, "genbank", digital=True, alphabet=alphabet
        ) as seqs_file:
            hits = list(pyhmmer.nhmmer(self.bmyD, seqs_file, cpus=1))[0]
            hits.sort()

        self.assertEqual(len(hits.reported), 2)
        with self.table("bmyD1.tbl") as table:
            self.assertTableEqual(hits, table)

    def test_bmyd_hmm_genome_block(self):
        alphabet = Alphabet.dna()

        path = resource_files(__package__).joinpath("data", "seqs", "CP000560.2.fna")
        with SequenceFile(path, "fasta", digital=True, alphabet=alphabet) as seqs_file:
            seqs = seqs_file.read_block()

        hits = next(pyhmmer.nhmmer(self.bmyD, seqs, cpus=1))
        hits.sort()

        self.assertEqual(len(hits.reported), 6)
        with self.table("bmyD2.tbl") as table:
            self.assertTableEqual(hits, table)

    def test_bmyd_hmm_genome_file(self):
        alphabet = Alphabet.dna()

        path = resource_files(__package__).joinpath("data", "seqs", "CP000560.2.fna")
        with SequenceFile(path, "fasta", digital=True, alphabet=alphabet) as seqs_file:
            hits = list(pyhmmer.nhmmer(self.bmyD, seqs_file, cpus=1))[0]
            hits.sort()

        self.assertEqual(len(hits.reported), 6)
        with self.table("bmyD2.tbl") as table:
            self.assertTableEqual(hits, table)

    def test_rf0001_genome_file(self):
        alphabet = Alphabet.rna()
        path = resource_files(__package__).joinpath("data", "seqs", "CP000560.2.fna")
        with SequenceFile(path, "fasta", digital=True, alphabet=alphabet) as seqs_file:
            hits = list(pyhmmer.nhmmer(self.rf00001, seqs_file, cpus=1))[0]
            hits.sort()

        self.assertAlmostEqual(hits[0].evalue, 1.4e-14)
        self.assertAlmostEqual(hits[9].evalue, 1.7e-13)

        strands = collections.Counter(hit.best_domain.strand for hit in hits.reported)
        self.assertEqual(strands["+"], 9)
        self.assertEqual(strands["-"], 1)

    def test_rf0001_genome_file_wlen_3878(self):
        alphabet = Alphabet.rna()
        path = resource_files(__package__).joinpath("data", "seqs", "CP000560.2.fna")
        with SequenceFile(path, "fasta", digital=True, alphabet=alphabet) as seqs_file:
            hits = list(pyhmmer.nhmmer(self.rf00001, seqs_file, cpus=1, window_length=3878))[0]
            hits.sort()

        self.assertAlmostEqual(hits[0].evalue, 3.1e-14)
        self.assertAlmostEqual(hits[9].evalue, 3.7e-13)

        strands = collections.Counter(hit.best_domain.strand for hit in hits.reported)
        self.assertEqual(strands["+"], 9)
        self.assertEqual(strands["-"], 1)

class TestHmmalign(unittest.TestCase):
    def setUp(self):
        self.tmpout = tempfile.NamedTemporaryFile(suffix=".hmm", delete=False).name

    def tearDown(self):
        os.remove(self.tmpout)

    @unittest.skipUnless(resource_files, "importlib.resources not available")
    def test_luxc(self):
        hmm_path = resource_files(__package__).joinpath("data", "hmms", "txt", "LuxC.hmm")
        seqs_path = resource_files(__package__).joinpath("data", "seqs", "LuxC.faa")
        ref_path = resource_files(__package__).joinpath("data", "msa", "LuxC.hmmalign.sto")

        with HMMFile(hmm_path) as hmm_file:
            hmm = hmm_file.read()
        with SequenceFile(seqs_path, digital=True, alphabet=hmm.alphabet) as seqs_file:
            seqs = seqs_file.read_block()
        with MSAFile(ref_path) as ref_file:
            ref = ref_file.read()

        msa = pyhmmer.hmmalign(hmm, seqs, trim=True)
        self.assertEqual(msa, ref)


class TestHMMScan(unittest.TestCase):
    @staticmethod
    def table(name):
        path = resource_files(__package__).joinpath("data", "tables", name)
        return path.open()

    @staticmethod
    def hmm_file(name):
        path = resource_files(__package__).joinpath("data", "hmms", "db", "{}.hmm".format(name))
        return HMMFile(path)

    @unittest.skipUnless(resource_files, "importlib.resources not available")
    def test_t2pks_block(self):
        # get paths to resources
        table_path = resource_files(__package__).joinpath("data", "tables", "t2pks.scan.tbl")
        db_path = resource_files(__package__).joinpath("data", "hmms", "db", "t2pks.hmm")
        seqs_path = resource_files(__package__).joinpath("data", "seqs", "PKSI.faa")

        # load expected results from the hmmscan table
        expected = {}
        with open(table_path) as table:
            lines = filter(lambda line: not line.startswith("#"), table)
            for query_name, query_lines in itertools.groupby(
                lines, key=lambda line: line.strip().split()[2]
            ):
                expected[query_name] = list(query_lines)

        # pre-load the profile database so that `hmmscan` will create a profile block
        with HMMFile(db_path) as hmm_file:
            hmms = list(hmm_file)

        # scan with the sequences and check the hits are equal
        with SequenceFile(
            seqs_path, digital=True, alphabet=hmms[0].alphabet
        ) as seqs_file:
            for hits in pyhmmer.hmmer.hmmscan(seqs_file, hmms, cpus=1):
                expected_lines = expected.get(hits.query_name.decode())
                if expected_lines is None:
                    self.assertEqual(len(hits), 0)
                    continue
                for line, hit in itertools.zip_longest(expected_lines, hits):
                    self.assertIsNot(line, None)
                    self.assertIsNot(hit, None)
                    fields = list(filter(None, line.strip().split(" ")))
                    self.assertEqual(hit.name.decode(), fields[0])
                    self.assertAlmostEqual(hit.score, float(fields[5]), delta=0.1)
                    self.assertAlmostEqual(hit.bias, float(fields[6]), delta=0.1)
                    self.assertAlmostEqual(hit.evalue, float(fields[4]), delta=0.1)

    @unittest.skipUnless(resource_files, "importlib.resources not available")
    def test_t2pks_file(self):
        # get paths to resources
        table_path = resource_files(__package__).joinpath("data", "tables", "t2pks.scan.tbl")
        db_path = resource_files(__package__).joinpath("data", "hmms", "db", "t2pks.hmm")
        seqs_path = resource_files(__package__).joinpath("data", "seqs", "PKSI.faa")

        # load expected results from the hmmscan table
        expected = {}
        with open(table_path) as table:
            lines = filter(lambda line: not line.startswith("#"), table)
            for query_name, query_lines in itertools.groupby(
                lines, key=lambda line: line.strip().split()[2]
            ):
                expected[query_name] = list(query_lines)

        # scan with the sequences and check the hits are equal, using the file
        # to read the profiles from (the files will be rewinded before each query)
        with SequenceFile(seqs_path, digital=True) as seqs_file:
            with HMMPressedFile(db_path) as pressed_file:
                for hits in pyhmmer.hmmer.hmmscan(seqs_file, pressed_file, cpus=1):
                    expected_lines = expected.get(hits.query_name.decode())
                    if expected_lines is None:
                        self.assertEqual(len(hits), 0)
                        continue
                    for line, hit in itertools.zip_longest(expected_lines, hits):
                        self.assertIsNot(line, None)
                        self.assertIsNot(hit, None)
                        fields = list(filter(None, line.strip().split(" ")))
                        self.assertEqual(hit.name.decode(), fields[0])
                        self.assertAlmostEqual(hit.score, float(fields[5]), delta=0.1)
                        self.assertAlmostEqual(hit.bias, float(fields[6]), delta=0.1)
                        self.assertAlmostEqual(hit.evalue, float(fields[4]), delta=0.1)
