import os
import unittest
import tempfile
import pkg_resources

from pyhmmer.plan7 import Pipeline, HMMFile, TopHits
from pyhmmer.easel import Alphabet, SequenceFile


class TestHmmsearch(unittest.TestCase):

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

        hmm_path = pkg_resources.resource_filename(__name__, "data/hmm/Thioesterase.hmm")
        with HMMFile(hmm_path) as hmm_file:
            hmm = next(hmm_file)

        seqs_path = pkg_resources.resource_filename(__name__, "data/seqs/938293.PRJEB85.HG003687.faa")
        with SequenceFile(seqs_path) as seqs_file:
            seqs = list(seqs_file)

        pipeline = Pipeline(hmm.alphabet)
        hits = pipeline.search(hmm, seqs)
        self.assertEqual(len(hits), 1)

        hit = hits[0]
        self.assertEqual(hit.name, b"938293.PRJEB85.HG003687_113")
        self.assertAlmostEqual(hit.score, 8.6, delta=0.1)      # printed with %6.1f
        self.assertAlmostEqual(hit.bias, 1.5, delta=0.1)       # printed with  %5.1f
        # self.assertAlmostEqual(hit.i_evalue, 0.096, delta=0.01)  # printed with %9.2g
        self.assertEqual(len(hit.domains), 1)

        domain = hit.domains[0]
        self.assertAlmostEqual(domain.score, 8.1, delta=0.1)
        self.assertAlmostEqual(domain.bias, 1.5, delta=0.1)
        # self.assertAlmostEqual(hit.i_evalue, 0.14, delta=0.01)  # printed with %9.2g
        self.assertAlmostEqual(domain.c_evalue, 6.5e-05, delta=0.01)  # printed with %9.2g
        self.assertEqual(domain.ali_from, 115)
        self.assertEqual(domain.ali_to, 129)
        self.assertEqual(domain.env_from, 115)
        self.assertEqual(domain.env_to, 129)
