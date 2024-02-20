import io
import itertools
import os
import pickle
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
class TestTopHits(unittest.TestCase):

    def assertAlignmentEqual(self, a1, a2):
        for attr in (
            "hmm_from",
            "hmm_to",
            "hmm_name",
            "hmm_accession",
            "hmm_sequence",
            "posterior_probabilities",
            "target_from",
            "target_to",
            "target_name",
            "target_sequence",
            "identity_sequence",
        ):
            self.assertEqual(
                getattr(a1, attr),
                getattr(a2, attr),
                "attribute {!r} differs".format(attr)
            )

    def assertDomainEqual(self, d1, d2):
        for attr in (
            "env_from",
            "env_to",
            "score",
            "bias",
            "correction",
            "envelope_score",
            "c_evalue",
            "i_evalue",
            "pvalue",
        ):
            self.assertEqual(
                getattr(d1, attr),
                getattr(d2, attr),
                "attribute {!r} differs".format(attr)
            )
        self.assertAlignmentEqual(d1.alignment, d2.alignment)

    def assertHitEqual(self, h1, h2):
        for attr in (
            "name",
            "accession",
            "description",
            "score",
            "pre_score",
            "sum_score",
            "bias",
            "evalue",
            "pvalue",
        ):
            self.assertEqual(
                getattr(h1, attr),
                getattr(h2, attr),
                "attribute {!r} differs".format(attr)
            )
        for d1, d2 in zip(h1.domains, h2.domains):
            self.assertDomainEqual(d1, d2)
        self.assertEqual(len(h1.domains), len(h2.domains))
        self.assertEqual(len(h1.domains.reported), len(h2.domains.reported))
        self.assertEqual(len(h1.domains.included), len(h2.domains.included))
        self.assertEqual(sum(d.reported for d in h1.domains), sum(d.reported for d in h2.domains))
        self.assertEqual(sum(d.included for d in h1.domains), sum(d.included for d in h2.domains))

    def assertHitsEqual(self, hits1, hits2):
        self.assertEqual(hits1.query_name, hits2.query_name)
        self.assertEqual(hits1.query_accession, hits2.query_accession)
        self.assertEqual(hits1.query_length, hits2.query_length)
        self.assertEqual(len(hits1), len(hits2))
        self.assertEqual(len(hits1.included), len(hits2.included))
        self.assertEqual(len(hits1.reported), len(hits2.reported))
        for h1, h2 in zip(hits1, hits2):
            self.assertHitEqual(h1, h2)

    @classmethod
    def setUpClass(cls):
        hmm_path = resource_files(__package__).joinpath("data", "hmms", "txt", "PF02826.hmm")
        seqs_path = resource_files(__package__).joinpath("data", "seqs", "938293.PRJEB85.HG003687.faa")

        with HMMFile(hmm_path) as hmm_file:
            cls.hmm = hmm_file.read()
        with SequenceFile(seqs_path, digital=True, alphabet=cls.hmm.alphabet) as seqs_file:
            cls.seqs = seqs_file.read_block()

        cls.pipeline = Pipeline(alphabet=cls.hmm.alphabet)
        cls._hits = cls.pipeline.search_hmm(cls.hmm, cls.seqs)

    def setUp(self):
        self.hits = self._hits.copy()

    def test_mode(self):
        pipeline = Pipeline(alphabet=self.hmm.alphabet)
        search_hits = pipeline.search_hmm(self.hmm, self.seqs)
        self.assertEqual(search_hits.mode, "search")

    def test_bool(self):
        self.assertFalse(pyhmmer.plan7.TopHits())
        self.assertTrue(self.hits)

    def test_index_error(self):
        # check out of bounds indexing is caught
        with self.assertRaises(IndexError):
            self.hits[len(self.hits) + 1]
        # check out of bounds negative indexing works
        with self.assertRaises(IndexError):
            self.hits[-len(self.hits) - 1]

    def test_index_negative(self):
        dom = self.hits[len(self.hits) - 1]
        dom_last = self.hits[-1]
        self.assertEqual(dom.name, dom_last.name)

    def test_Z(self):
        empty = TopHits()
        self.assertEqual(empty.Z, 0)
        self.assertEqual(self.hits.Z, len(self.seqs))

    def test_strand(self):
        empty = TopHits()
        self.assertIs(empty.strand, None)

    def test_searched_sequences(self):
        empty = TopHits()
        self.assertEqual(empty.searched_sequences, 0)
        self.assertEqual(self.hits.searched_sequences, len(self.seqs))

    def test_searched_nodes(self):
        empty = TopHits()
        self.assertEqual(empty.searched_nodes, 0)
        self.assertEqual(self.hits.searched_nodes, self.hmm.M)

    def test_searched_residues(self):
        empty = TopHits()
        self.assertEqual(empty.searched_residues, 0)
        self.assertEqual(self.hits.searched_residues, sum(map(len, self.seqs)))

    def test_searched_models(self):
        empty = TopHits()
        self.assertEqual(empty.searched_sequences, 0)
        self.assertEqual(self.hits.searched_models, 1)

    def test_merge_empty(self):
        empty = TopHits()
        self.assertFalse(empty.long_targets)
        self.assertEqual(empty.Z, 0.0)
        self.assertEqual(empty.domZ, 0.0)

        empty2 = empty.merge()
        self.assertFalse(empty2.long_targets)
        self.assertEqual(empty2.Z, 0.0)
        self.assertEqual(empty2.domZ, 0.0)

        merged_empty = empty.merge(TopHits())
        self.assertHitsEqual(merged_empty, empty)
        self.assertEqual(merged_empty.searched_residues, 0)
        self.assertEqual(merged_empty.searched_sequences, 0)
        self.assertFalse(merged_empty.long_targets)
        self.assertFalse(self.hits.long_targets)
        self.assertEqual(merged_empty.Z, 0.0)
        self.assertEqual(merged_empty.domZ, 0.0)

        merged = empty.merge(self.hits)
        self.assertEqual(merged.query_name, self.hits.query_name)
        self.assertEqual(merged.query_length, self.hits.query_length)
        self.assertEqual(merged.query_accession, self.hits.query_accession)
        self.assertEqual(merged.E, self.hits.E)
        self.assertHitsEqual(merged, self.hits)

    def test_merge_pipeline(self):
        pipeline = Pipeline(alphabet=self.hmm.alphabet)

        hits1 = pipeline.search_hmm(self.hmm, self.seqs[:1000])
        self.assertEqual(hits1.Z, 1000)
        hits2 = pipeline.search_hmm(self.hmm, self.seqs[1000:2000])
        self.assertEqual(hits2.Z, 1000)
        hits3 = pipeline.search_hmm(self.hmm, self.seqs[2000:])
        self.assertEqual(hits3.Z, len(self.seqs) - 2000)

        merged = hits1.merge(hits2, hits3)
        self.assertEqual(merged.mode, "search")

        hits = pipeline.search_hmm(self.hmm, self.seqs)
        self.assertEqual(len(hits1) + len(hits2) + len(hits3), len(hits))
        self.assertEqual(len(merged), len(hits))

        self.assertEqual(merged.searched_sequences, hits.searched_sequences)
        self.assertEqual(merged.Z, hits.Z)
        self.assertEqual(merged.domZ, hits.domZ)
        self.assertEqual(merged.query_name, hits.query_name)
        self.assertEqual(merged.query_length, hits.query_length)
        self.assertEqual(merged.query_accession, hits.query_accession)
        self.assertEqual(merged.E, hits.E)
        self.assertEqual(merged.domE, hits.domE)

        for hit in merged:
            nincluded = sum(domain.included for domain in hit.domains)
            self.assertEqual(nincluded, len(hit.domains.included))
            nreported = sum(domain.reported for domain in hit.domains)
            self.assertEqual(nreported, len(hit.domains.reported))

        self.assertHitsEqual(merged, hits)

        # number of reported & included hits in the final merged hits may be
        # different of the sum here, because we used E-value thresholds, which
        # are affected by how many hits were effectively reported

    def test_merge_pipeline_byscore(self):
        pipeline = Pipeline(alphabet=self.hmm.alphabet, T=10.0, domT=10.0, incT=10.0, incdomT=10.0)

        hits1 = pipeline.search_hmm(self.hmm, self.seqs[:1000])
        self.assertEqual(hits1.Z, 1000)
        hits2 = pipeline.search_hmm(self.hmm, self.seqs[1000:2000])
        self.assertEqual(hits2.Z, 1000)
        hits3 = pipeline.search_hmm(self.hmm, self.seqs[2000:])
        self.assertEqual(hits3.Z, len(self.seqs) - 2000)

        merged = hits1.merge(hits2, hits3)
        self.assertEqual(merged.mode, "search")

        hits = pipeline.search_hmm(self.hmm, self.seqs)
        self.assertEqual(len(hits1) + len(hits2) + len(hits3), len(hits))
        self.assertEqual(len(merged), len(hits))

        self.assertEqual(merged.searched_sequences, hits.searched_sequences)
        self.assertEqual(merged.Z, hits.Z)
        self.assertEqual(merged.domZ, hits.domZ)
        self.assertEqual(merged.query_name, hits.query_name)
        self.assertEqual(merged.query_length, hits.query_length)
        self.assertEqual(merged.query_accession, hits.query_accession)
        self.assertEqual(merged.E, hits.E)
        self.assertEqual(merged.domE, hits.domE)

        for hit in merged:
            nincluded = sum(domain.included for domain in hit.domains)
            self.assertEqual(nincluded, len(hit.domains.included))
            nreported = sum(domain.reported for domain in hit.domains)
            self.assertEqual(nreported, len(hit.domains.reported))

        self.assertHitsEqual(merged, hits)

        # number of reported & included hits in the final merged hits should
        # be equal to the sum of reported & included hits in each sub hits,
        # because we used bitscore thresholds, which are independent of
        # the number of hits

        reported = lambda hits: sum(len(hit.domains.reported) for hit in hits)
        included = lambda hits: sum(len(hit.domains.reported) for hit in hits)
        hits_nreported = reported(hits1) + reported(hits2) + reported(hits3)
        hits_nincluded = included(hits1) + included(hits2) + included(hits3)
        self.assertEqual(reported(merged), hits_nreported)
        self.assertEqual(included(merged), hits_nincluded)

    def test_merged_pipeline_fixed_Z(self):
        pipeline = Pipeline(alphabet=self.hmm.alphabet, Z=200.0)
        hits1 = pipeline.search_hmm(self.hmm, self.seqs[:1000])
        hits2 = pipeline.search_hmm(self.hmm, self.seqs[1000:2000])
        hits3 = pipeline.search_hmm(self.hmm, self.seqs[2000:])

        self.assertEqual(hits1.Z, 200.0)
        self.assertEqual(hits2.Z, 200.0)
        self.assertEqual(hits3.Z, 200.0)

        merged = hits1.merge(hits2, hits3)
        self.assertEqual(merged.Z, 200.0)

        hits = pipeline.search_hmm(self.hmm, self.seqs)
        self.assertHitsEqual(merged, hits)

    def test_copy(self):
        copy = self.hits.copy()
        self.assertHitsEqual(copy, self.hits)

    def test_sort(self):
        # check the hits are sorted by default
        self.assertTrue(self.hits.is_sorted())

        # check sorting
        self.hits.sort()
        self.assertTrue(self.hits.is_sorted())
        self.assertFalse(self.hits.is_sorted(by="seqidx"))

        # check sorting using a different key
        self.hits.sort(by="seqidx")
        self.assertTrue(self.hits.is_sorted(by="seqidx"))
        self.assertFalse(self.hits.is_sorted(by="key"))

        # check error otherwise
        self.assertRaises(ValueError, self.hits.sort, by="nonsense")
        self.assertRaises(ValueError, self.hits.is_sorted, by="nonsense")

    def test_sort_while_hit(self):
        # check sorting while a Hit instance exists does not crash, and that
        # the Hit still points to the right data
        hit = self.hits[-1]
        self.assertEqual(hit.name, b"938293.PRJEB85.HG003687_187")
        self.hits.sort(by="seqidx")
        if self.hits[-1].name != b"938293.PRJEB85.HG003687_187":
            self.assertEqual(hit.name, b"938293.PRJEB85.HG003687_187")

    def test_threshold(self):
        # after thresholding, one of the hits will not be included.
        self.assertEqual(len(self.hits.included), 15)
        self.assertEqual(len(self.hits.reported), 22)

    def test_to_msa(self):
        msa = self.hits.to_msa(self.hmm.alphabet)
        self.assertIsInstance(msa, TextMSA)
        unique_names = { seq.name.split(b"/")[0] for seq in msa.sequences }
        self.assertEqual(len(unique_names), len(self.hits.included))

        msa_d = self.hits.to_msa(
            self.hmm.alphabet,
            trim=True,
            digitize=True,
            all_consensus_cols=True
        )
        self.assertIsInstance(msa_d, DigitalMSA)

    def test_pickle(self):
        pickled = pickle.loads(pickle.dumps(self.hits))
        self.assertHitsEqual(pickled, self.hits)

    def test_query_name(self):
        self.assertEqual(self.hits.query_name, self.hmm.name)

    def test_query_accession(self):
        self.assertEqual(self.hits.query_accession, self.hmm.accession)

    def test_query_length(self):
        self.assertEqual(self.hits.query_length, self.hmm.M)

    def test_write_target(self):
        buffer = io.BytesIO()
        self.hits.write(buffer, format="targets")
        lines = buffer.getvalue().decode().splitlines()

        with resource_files(__package__).joinpath("data", "tables", "PF02826.tbl").open() as f:
            expected = f.read().splitlines()
        while expected[-1].startswith("#"):
            expected.pop()

        self.assertMultiLineEqual("\n".join(lines), "\n".join(expected))

    def test_write_domains(self):
        buffer = io.BytesIO()
        self.hits.write(buffer, format="domains")
        lines = buffer.getvalue().decode().splitlines()

        with resource_files(__package__).joinpath("data", "tables", "PF02826.domtbl").open() as f:
            expected = f.read().splitlines()
        while expected[-1].startswith("#"):
            expected.pop()

        self.assertMultiLineEqual("\n".join(lines), "\n".join(expected))

    def test_manual_report(self):
        self.assertEqual(len(self.hits.reported), 22)
        self.assertTrue(self.hits[0].reported)

        self.hits[0].reported = False
        self.assertEqual(len(self.hits.reported), 21)
        self.assertFalse(self.hits[0].reported)

        self.hits[0].reported = True
        self.assertEqual(len(self.hits.reported), 22)
        self.assertTrue(self.hits[0].reported)

    def test_manual_inclusion(self):
        self.assertEqual(len(self.hits.included), 15)
        self.assertEqual(len(self.hits.reported), 22)
        self.assertTrue(self.hits[0].included)
        self.assertTrue(self.hits[0].reported)

        self.hits[0].included = False
        self.assertEqual(len(self.hits.included), 14)
        self.assertEqual(len(self.hits.reported), 22)
        self.assertFalse(self.hits[0].included)
        self.assertTrue(self.hits[0].reported)

        self.hits[0].included = True
        self.assertEqual(len(self.hits.included), 15)
        self.assertEqual(len(self.hits.reported), 22)
        self.assertTrue(self.hits[0].included)
        self.assertTrue(self.hits[0].reported)

    def test_manual_drop(self):
        self.assertEqual(len(self.hits.included), 15)
        self.assertEqual(len(self.hits.reported), 22)
        self.assertTrue(self.hits[0].included)
        self.assertTrue(self.hits[0].reported)

        self.hits[0].dropped = True
        self.assertEqual(len(self.hits.included), 14)
        self.assertEqual(len(self.hits.reported), 22)
        self.assertTrue(self.hits[0].dropped)
        self.assertFalse(self.hits[0].included)
        self.assertTrue(self.hits[0].reported)

        self.hits[0].dropped = False
        self.assertEqual(len(self.hits.included), 14)
        self.assertEqual(len(self.hits.reported), 22)
        self.assertFalse(self.hits[0].dropped)
        self.assertFalse(self.hits[0].included)
        self.assertTrue(self.hits[0].reported)

    def test_manual_duplicate(self):
        self.assertEqual(len(self.hits.included), 15)
        self.assertEqual(len(self.hits.reported), 22)
        self.assertTrue(self.hits[0].included)
        self.assertTrue(self.hits[0].reported)

        self.hits[0].duplicate = True
        self.assertEqual(len(self.hits.included), 14)
        self.assertEqual(len(self.hits.reported), 21)
        self.assertFalse(self.hits[0].included)
        self.assertFalse(self.hits[0].reported)
        self.assertTrue(self.hits[0].duplicate)

        self.hits[0].duplicate = False
        self.assertEqual(len(self.hits.included), 14)
        self.assertEqual(len(self.hits.reported), 21)
        self.assertFalse(self.hits[0].included)
        self.assertFalse(self.hits[0].reported)
        self.assertFalse(self.hits[0].duplicate)
