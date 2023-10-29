import abc
import pickle
import unittest

from pyhmmer.errors import AlphabetMismatch
from pyhmmer.easel import (
    Alphabet,
    GeneticCode,
    SequenceBlock,
    TextSequenceBlock,
    DigitalSequenceBlock,
    Sequence,
    TextSequence,
    DigitalSequence,
)


class _TestSequenceBlock(abc.ABC):

    @abc.abstractmethod
    def _new_sequence(self, name, seq):
        pass

    @abc.abstractmethod
    def _new_block(self, sequences=()):
        pass

    def test_truth(self):
        self.assertTrue(not self._new_block())
        seq1 = self._new_sequence(b"seq1", "ATGC")
        block = self._new_block([seq1])
        self.assertTrue(self._new_block([seq1]))

    def test_identity(self):
        self.assertIsNot(self._new_block(), self._new_block())

    def test_len(self):
        seq1 = self._new_sequence(b"seq1", "ATGC")
        seq2 = self._new_sequence(b"seq2", "ATGCA")
        seq3 = self._new_sequence(b"seq3", "TTGA")

        self.assertEqual(len(self._new_block([])), 0)
        self.assertEqual(len(self._new_block([seq1])), 1)
        self.assertEqual(len(self._new_block([seq1, seq2, seq3])), 3)

    def test_append(self):
        block = self._new_block()
        self.assertEqual(len(block), 0)
        seq1 = self._new_sequence(b"seq1", "ATGC")
        block.append(seq1)
        self.assertEqual(len(block), 1)
        self.assertEqual(block[0], seq1)
        seq2 = self._new_sequence(b"seq2", "ATGCA")
        block.append(seq2)
        self.assertEqual(len(block), 2)
        self.assertEqual(block[1], seq2)

    def test_iter(self):
        seq1 = self._new_sequence(b"seq1", "ATGC")
        seq2 = self._new_sequence(b"seq2", "ATGCA")
        seq3 = self._new_sequence(b"seq3", "TTGA")

        block = self._new_block([seq1, seq2, seq3])
        self.assertEqual(len(block), 3)

        iterator = iter(block)
        self.assertIs(next(iterator), seq1)
        self.assertIs(next(iterator), seq2)
        self.assertIs(next(iterator), seq3)
        self.assertRaises(StopIteration, next, iterator)

    def test_clear(self):
        seq1 = self._new_sequence(b"seq1", "ATGC")
        seq2 = self._new_sequence(b"seq2", "ATGCA")

        block = self._new_block([seq1, seq2])
        self.assertEqual(len(block), 2)
        block.clear()
        self.assertEqual(len(block), 0)
        block.clear()
        self.assertEqual(len(block), 0)

        block = self._new_block()
        self.assertEqual(len(block), 0)
        block.clear()
        self.assertEqual(len(block), 0)

    def test_getitem(self):
        seq1 = self._new_sequence(b"seq1", "ATGC")
        seq2 = self._new_sequence(b"seq2", "ATGCA")

        block = self._new_block([seq1, seq2])
        self.assertEqual(len(block), 2)

        self.assertIs(block[0], seq1)
        self.assertIs(block[1], seq2)
        self.assertIs(block[-1], seq2)
        self.assertIs(block[-2], seq1)

        with self.assertRaises(IndexError):
            block[3]
        with self.assertRaises(IndexError):
            block[-5]

    def test_setitem(self):
        seq1 = self._new_sequence(b"seq1", "ATGC")
        seq2 = self._new_sequence(b"seq2", "ATGCA")
        seq3 = self._new_sequence(b"seq3", "TTGA")

        block = self._new_block([seq1, seq2])
        self.assertEqual(len(block), 2)

        block[0] = seq3
        self.assertIs(block[0], seq3)
        self.assertIs(block[1], seq2)

    def test_remove(self):
        seq1 = self._new_sequence(b"seq1", "ATGC")
        seq2 = self._new_sequence(b"seq2", "ATGCA")

        block = self._new_block([seq1, seq2])
        block.remove(seq1)
        self.assertEqual(len(block), 1)
        self.assertIs(block[0], seq2)

    def test_index(self):
        seq1 = self._new_sequence(b"seq1", "ATGC")
        seq2 = self._new_sequence(b"seq2", "ATGCA")
        seq3 = self._new_sequence(b"seq3", "TTGA")

        block = self._new_block([seq1, seq2])
        self.assertEqual(len(block), 2)

        self.assertEqual(block.index(seq1), 0)
        self.assertEqual(block.index(seq2), 1)

        self.assertRaises(ValueError, block.index, seq3)
        self.assertRaises(ValueError, block.index, seq1, start=1)

    def test_pop(self):
        seq1 = self._new_sequence(b"seq1", "ATGC")
        seq2 = self._new_sequence(b"seq2", "ATGCA")
        seq3 = self._new_sequence(b"seq3", "TTGA")

        block = self._new_block([seq1, seq2, seq3])
        self.assertEqual(len(block), 3)

        self.assertIs(block.pop(), seq3)
        self.assertIs(block.pop(0), seq1)

        self.assertEqual(len(block), 1)
        self.assertIs(block[0], seq2)

        self.assertIs(block.pop(-1), seq2)
        self.assertEqual(len(block), 0)
        self.assertRaises(IndexError, block.pop)

        empty = self._new_block()
        self.assertRaises(IndexError, block.pop)

    def test_setitem_slice(self):
        seq1 = self._new_sequence(b"seq1", "ATGC")
        seq2 = self._new_sequence(b"seq2", "ATGCA")
        seq3 = self._new_sequence(b"seq3", "TTGA")

        block = self._new_block([seq1, seq2])
        self.assertEqual(len(block), 2)
        self.assertIs(block[0], seq1)
        self.assertIs(block[1], seq2)

        block[:] = [seq1, seq2, seq3]
        self.assertEqual(len(block), 3)
        self.assertIs(block[0], seq1)
        self.assertIs(block[1], seq2)
        self.assertIs(block[2], seq3)

        block[2:5] = [seq3, seq3, seq3]
        self.assertEqual(len(block), 5)
        self.assertIs(block[3], seq3)
        self.assertIs(block[4], seq3)

    def test_delitem_slice(self):
        seq1 = self._new_sequence(b"seq1", "ATGC")
        seq2 = self._new_sequence(b"seq2", "ATGCA")
        seq3 = self._new_sequence(b"seq3", "TTGA")

        block = self._new_block([seq1, seq2, seq3])
        self.assertEqual(len(block), 3)
        del block[:]
        self.assertEqual(len(block), 0)

        block = self._new_block([seq1, seq2, seq3])
        self.assertEqual(len(block), 3)
        del block[::2]
        self.assertEqual(len(block), 1)
        self.assertIs(block[0], seq2)

    def test_contains(self):
        seq1 = self._new_sequence(b"seq1", "ATGC")
        seq2 = self._new_sequence(b"seq2", "ATGCA")
        seq3 = self._new_sequence(b"seq3", "TTGA")

        block = self._new_block([seq1, seq2])
        self.assertEqual(len(block), 2)
        self.assertIn(seq1, block)
        self.assertIn(seq2, block)
        self.assertNotIn(seq3, block)

        self.assertNotIn(42, block)
        self.assertNotIn(object(), block)

    def test_largest(self):
        seq1 = self._new_sequence(b"seq1", "ATGC")
        seq2 = self._new_sequence(b"seq2", "ATGCA")
        seq3 = self._new_sequence(b"seq3", "TTGA")

        block = self._new_block([seq1, seq2])
        self.assertIs(block.largest(), seq2)

        block.pop(1)
        self.assertIs(block.largest(), seq1)

        with self.assertRaises(ValueError):
            block = self._new_block([])
            block.largest()

    def test_copy(self):
        seq1 = self._new_sequence(b"seq1", "ATGC")
        seq2 = self._new_sequence(b"seq2", "ATGCA")
        seq3 = self._new_sequence(b"seq3", "TTGA")

        block = self._new_block([seq1, seq2, seq3])
        block_copy = block.copy()

        self.assertEqual(len(block), len(block_copy))
        self.assertEqual(block, block_copy)

        self.assertIs(block[0], block_copy[0])
        self.assertIs(block[1], block_copy[1])
        self.assertIs(block[2], block_copy[2])

        block.remove(seq3)
        self.assertEqual(len(block), 2)
        self.assertEqual(len(block_copy), 3)

    def test_pickle(self):
        seq1 = self._new_sequence(b"seq1", "ATGC")
        seq2 = self._new_sequence(b"seq2", "ATGCA")
        seq3 = self._new_sequence(b"seq3", "TTGA")

        block = self._new_block([seq1, seq2, seq3])
        block_pickled = pickle.loads(pickle.dumps(block))
        
        self.assertEqual(len(block), len(block_pickled))
        self.assertEqual(block[0], block_pickled[0])
        self.assertEqual(block[1], block_pickled[1])
        self.assertEqual(block[2], block_pickled[2])


class TestTextSequenceBlock(_TestSequenceBlock, unittest.TestCase):

    def _new_sequence(self, name, seq):
        return TextSequence(name=name, sequence=seq)

    def _new_block(self, sequences=()):
        return TextSequenceBlock(sequences)

    def test_digitize(self):
        alphabet = Alphabet.dna()

        seq1 = self._new_sequence(b"seq1", "ATGC")
        seq2 = self._new_sequence(b"seq2", "ATGCA")

        block = self._new_block([seq1, seq2])
        self.assertEqual(len(block), 2)

        dblock = block.digitize(alphabet)
        self.assertEqual(len(dblock), 2)
        self.assertEqual(dblock[0], seq1.digitize(alphabet))
        self.assertEqual(dblock[1], seq2.digitize(alphabet))


class TestDigitalSequenceBlock(_TestSequenceBlock, unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.alphabet = Alphabet.dna()

    def _new_sequence(
        self,
        name,
        seq,
        description=None,
        source=None,
    ):
        return TextSequence(
            name=name,
            sequence=seq,
            description=description,
            source=source,
        ).digitize(self.alphabet)

    def _new_block(self, sequences=()):
        return DigitalSequenceBlock(self.alphabet, sequences)

    def test_textize(self):
        seq1 = self._new_sequence(b"seq1", "ATGC")
        seq2 = self._new_sequence(b"seq2", "ATGCA")

        block = self._new_block([seq1, seq2])
        self.assertEqual(len(block), 2)

        tblock = block.textize()
        self.assertEqual(len(tblock), 2)
        self.assertEqual(tblock[0], seq1.textize())
        self.assertEqual(tblock[1], seq2.textize())

    def test_translate(self):
        seq1 = self._new_sequence(
            b"seq1",
            "ATGCTG",
            description=b"one",
            source=b"some_source"
        )
        seq2 = self._new_sequence(
            b"seq2",
            "ATGCCC",
            description=b"two",
            source=b"other_source"
        )
        seq3 = self._new_sequence(
            b"seq3",
            "ATTCTGATGATA",
        )

        block = self._new_block([seq1, seq2, seq3])
        prots = block.translate().textize()

        self.assertEqual(len(prots), 3)
        # test sequences
        self.assertEqual(prots[0].sequence, "ML")
        self.assertEqual(prots[1].sequence, "MP")
        self.assertEqual(prots[2].sequence, "ILMI")
        # test names
        self.assertEqual(prots[0].name, b"seq1")
        self.assertEqual(prots[1].name, b"seq2")
        self.assertEqual(prots[2].name, b"seq3")
        # test other metadata fields
        self.assertEqual(prots[0].description, b"one")
        self.assertEqual(prots[1].description, b"two")
        self.assertEqual(prots[2].description, b"")
        self.assertEqual(prots[0].source, b"some_source")
        self.assertEqual(prots[1].source, b"other_source")
        self.assertEqual(prots[2].source, b"")

    def test_translate_empty(self):
        gencode = GeneticCode()
        block = DigitalSequenceBlock(gencode.nucleotide_alphabet)
        prots = block.translate(gencode)
        self.assertEqual(len(prots), 0)

    def test_translate_alphabet_mismatch(self):
        gencode = GeneticCode()
        block = DigitalSequenceBlock(gencode.amino_alphabet)
        with self.assertRaises(AlphabetMismatch):
            prots = block.translate(gencode)