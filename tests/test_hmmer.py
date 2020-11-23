import abc
import io
import itertools
import os
import unittest
import tempfile
import threading
import pkg_resources

import pyhmmer
from pyhmmer.plan7 import Pipeline, HMMFile, TopHits
from pyhmmer.easel import Alphabet, SequenceFile, TextSequence


class _TestSearch(metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def get_hits(self, hmm, sequences):
        return NotImplemented

    @staticmethod
    def table(name):
        bin_stream = pkg_resources.resource_stream(__name__, "data/tables/{}".format(name))
        return io.TextIOWrapper(bin_stream)

    @staticmethod
    def hmm_file(name):
        path = pkg_resources.resource_filename(__name__, "data/hmm/{}.hmm".format(name))
        return HMMFile(path)

    @staticmethod
    def seqs_file(name):
        seqs_path = pkg_resources.resource_filename(__name__, "data/seqs/{}.faa".format(name))
        return SequenceFile(seqs_path)

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
            hmm = next(hmm_file)
        with self.seqs_file("938293.PRJEB85.HG003687") as seqs_file:
            seqs = [seq.digitize(hmm.alphabet) for seq in seqs_file]

        hits = self.get_hits(hmm, seqs)
        self.assertEqual(len(hits), 1)

        hits.sort()

        hit = hits[0]
        self.assertEqual(hit.name, b"938293.PRJEB85.HG003687_113")
        self.assertAlmostEqual(hit.score, 8.6, delta=0.1)      # printed with %6.1f
        self.assertAlmostEqual(hit.bias, 1.5, delta=0.1)       # printed with  %5.1f
        self.assertAlmostEqual(hit.evalue, 0.096, delta=0.01)  # printed with %9.2g
        self.assertEqual(len(hit.domains), 1)

        domain = hit.domains[0]
        self.assertAlmostEqual(domain.score, 8.1, delta=0.1)
        self.assertAlmostEqual(domain.bias, 1.5, delta=0.1)
        self.assertAlmostEqual(domain.i_evalue, 0.14, delta=0.01)  # printed with %9.2g
        self.assertAlmostEqual(domain.c_evalue, 6.5e-05, delta=0.1e-5)  # printed with %9.2g
        self.assertEqual(domain.ali_from, 115)
        self.assertEqual(domain.ali_to, 129)
        self.assertEqual(domain.env_from, 115)
        self.assertEqual(domain.env_to, 129)

    def test_pf02826(self):
        with self.hmm_file("PF02826") as hmm_file:
            hmm = next(hmm_file)
        with self.seqs_file("938293.PRJEB85.HG003687") as seqs_file:
            seqs = [seq.digitize(hmm.alphabet) for seq in seqs_file]

        pipeline = Pipeline(alphabet=hmm.alphabet)
        hits = pipeline.search(hmm, seqs)
        self.assertEqual(len(hits), 16)

        hits.sort()

        with self.table("PF02826.tbl") as table:
            lines = filter(lambda line: not line.startswith("#"), table)
            for line, hit in itertools.zip_longest(lines, hits):
                fields = list(filter(None, line.strip().split(" ")))
                self.assertIsNot(line, None)
                self.assertIsNot(hit, None)
                self.assertEqual(hit.name.decode(), fields[0])
                self.assertAlmostEqual(hit.score, float(fields[5]), delta=0.1)
                self.assertAlmostEqual(hit.bias, float(fields[6]), delta=0.1)
                self.assertAlmostEqual(hit.evalue, float(fields[4]), delta=0.1)


class TestHmmsearch(_TestSearch, unittest.TestCase):

    def get_hits(self, hmm, seqs):
        return pyhmmer.hmmsearch([hmm], seqs)[0]

    def test_background_error(self):
        # check that errors occuring in worker threads are recovered and raised
        # in the main threads (a common error is mismatching the HMM and the
        # sequence alphabets).
        seqs = [ TextSequence().digitize(Alphabet.dna()) ]
        with self.hmm_file("PF02826") as hmm_file:
            hmm = next(hmm_file)
        self.assertRaises(ValueError, self.get_hits, hmm, seqs)

class TestHmmsearchSingle(TestHmmsearch, unittest.TestCase):

    def get_hits(self, hmm, seqs):
        return pyhmmer.hmmsearch([hmm], seqs, cpus=1)[0]


class TestPipelinesearch(_TestSearch, unittest.TestCase):

    def get_hits(self, hmm, seqs):
        pipeline = Pipeline(alphabet=hmm.alphabet)
        hits = pipeline.search(hmm, seqs)
        return hits
