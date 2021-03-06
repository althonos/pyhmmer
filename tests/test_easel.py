import copy
import gc
import os
import unittest
import tempfile

from pyhmmer import easel


class TestAlphabet(unittest.TestCase):

    def test_eq(self):
        self.assertEqual(easel.Alphabet.dna(), easel.Alphabet.dna())
        self.assertEqual(easel.Alphabet.rna(), easel.Alphabet.rna())
        self.assertEqual(easel.Alphabet.amino(), easel.Alphabet.amino())
        self.assertNotEqual(easel.Alphabet.amino(), easel.Alphabet.dna())
        self.assertNotEqual(easel.Alphabet.amino(), object())


class TestMSAFile(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.formats_folder = os.path.realpath(os.path.join(
            __file__, os.pardir, os.pardir, "vendor", "easel", "formats"
        ))
        cls.testsuite_folder = os.path.realpath(os.path.join(
            __file__, os.pardir, os.pardir, "vendor", "easel", "testsuite"
        ))

    def test_author(self):
        trna_5 = os.path.join(self.testsuite_folder, "trna-5.stk")
        with easel.MSAFile(trna_5) as f:
            msa = f.read()

        self.assertEqual(msa.author, b"Infernal 0.1")


    def test_readformat_stockholm(self):
        stockholm = os.path.join(self.formats_folder, "stockholm.1")

        # check reading with specified format works
        with easel.MSAFile(stockholm, "stockholm") as f:
            msa = f.read()
            self.assertNotEqual(msa.name, b"")

        # check reading without specifying the format works too
        with easel.MSAFile(stockholm) as f:
            msa2 = f.read()
            self.assertNotEqual(msa2.name, b"")

        # check reading while giving another file format fails
        with easel.MSAFile(stockholm, "clustal") as f:
            self.assertRaises(ValueError, f.read)


class TestSequenceFile(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.formats_folder = os.path.realpath(os.path.join(
            __file__, os.pardir, os.pardir, "vendor", "easel", "formats"
        ))

    def test_guess_alphabet(self):
        embl = os.path.join(self.formats_folder, "embl")
        with easel.SequenceFile(embl, "embl") as f:
            self.assertEqual(f.guess_alphabet(), easel.Alphabet.dna())

    def test_init_error_unknownformat(self):
        with self.assertRaises(ValueError):
            _file = easel.SequenceFile("file.x", format="nonsense")

    def test_init_error_empty(self):
        with tempfile.NamedTemporaryFile() as empty:
            # without a format argument, we can't determine the file format
            self.assertRaises(ValueError, easel.SequenceFile, empty.name)
            # with a format argument, we don't expect the file to be empty
            self.assertRaises(EOFError, easel.SequenceFile, empty.name, "fasta")

    def test_init_error_filenotfound(self):
        self.assertRaises(FileNotFoundError, easel.SequenceFile, "path/to/missing/file.fa")

    def test_iter(self):
        fasta = os.path.join(self.formats_folder, "fasta")

        with easel.SequenceFile(fasta) as f:
            seq1, seq2 = list(f)

        self.assertEqual(seq1.name, b"SNRPA_DROME")
        self.assertEqual(seq2.name, b"SNRPA_HUMAN")

    def test_parse(self):
        seq = easel.SequenceFile.parse(b"> EcoRI\nGAATTC\n", format="fasta")
        self.assertEqual(seq.name, b"EcoRI")

    def test_parse_error_unknownformat(self):
        with self.assertRaises(ValueError):
            _file = easel.SequenceFile.parse(b"> EcoRI\nGAATTC\n", format="nonsense")

    def test_read_skip_info(self):
        fasta = os.path.join(self.formats_folder, "fasta")

        with easel.SequenceFile(fasta) as f:
            seq = f.read()
            self.assertEqual(seq.name, b"SNRPA_DROME")
            self.assertNotEqual(seq.sequence, "")

        with easel.SequenceFile(fasta) as f:
            seq = f.read(skip_info=True)
            self.assertEqual(seq.name, b"")
            self.assertNotEqual(seq.sequence, "")

    def test_read_skip_sequence(self):
        fasta = os.path.join(self.formats_folder, "fasta")

        with easel.SequenceFile(fasta) as f:
            seq = f.read()
            self.assertEqual(seq.name, b"SNRPA_DROME")
            self.assertNotEqual(seq.sequence, "")

        with easel.SequenceFile(fasta) as f:
            seq = f.read(skip_sequence=True)
            self.assertEqual(seq.name, b"SNRPA_DROME")
            self.assertEqual(seq.sequence, "")

    def test_readformat_fasta(self):
        for filename in ["fasta", "fasta.2", "fasta.odd.1"]:
            fasta = os.path.join(self.formats_folder, filename)

            # check reading a FASTA file works
            with easel.SequenceFile(fasta, "fasta") as f:
                seq = f.read()
                self.assertNotEqual(seq.name, b"")

            # check reading a FASTA file without specifying the format works too
            with easel.SequenceFile(fasta) as f:
                seq2 = f.read()
                self.assertNotEqual(seq2.name, b"")

            # check reading a FASTA file while giving another format fails
            with easel.SequenceFile(fasta, "embl") as f:
                self.assertRaises(ValueError, f.read)

    def test_readformat_fasta_error(self):
        fasta = os.path.join(self.formats_folder, "fasta.bad.1")
        with easel.SequenceFile(fasta, "fasta") as f:
            self.assertRaises(ValueError, f.read)

    def test_readformat_embl(self):
        embl = os.path.join(self.formats_folder, "embl")

        # check reading an EMBL file works
        with easel.SequenceFile(embl, "embl") as f:
            seq = f.read()
            self.assertNotEqual(seq.accession, b"")

        # check reading an EMBL file without specifying the format works too
        with easel.SequenceFile(embl) as f:
            seq2 = f.read()
            self.assertNotEqual(seq2.accession, b"")

        # check reading an EMBL file while giving another format fails
        with easel.SequenceFile(embl, "fasta") as f:
            self.assertRaises(ValueError, f.read)

    def test_readformat_genbank(self):
        gbk = os.path.join(self.formats_folder, "genbank")

        # check reading an EMBL file works
        with easel.SequenceFile(gbk, "genbank") as f:
            seq = f.read()
            self.assertNotEqual(seq.accession, b"")

        # check reading an GenBank file without specifying the format works too
        with easel.SequenceFile(gbk) as f:
            seq2 = f.read()
            self.assertNotEqual(seq2.accession, b"")

        # check reading an GenBank file while giving another format fails
        with easel.SequenceFile(gbk, "fasta") as f:
            self.assertRaises(ValueError, f.read)


class TestDigitalSequence(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.abc = easel.Alphabet.dna()

    def test_init_empty(self):
        seq = easel.DigitalSequence(self.abc)
        self.assertEqual(seq.name, b"")
        self.assertEqual(seq.sequence, bytearray())
        self.assertEqual(len(seq), 0)

    def test_init_kwargs(self):
        arr = bytearray([0, 1, 2, 3])
        seq = easel.DigitalSequence(self.abc, name=b"TEST", sequence=arr)
        self.assertEqual(seq.name, b"TEST")
        self.assertEqual(bytearray(seq.sequence), arr)
        self.assertEqual(len(seq), 4)

    def test_setter_name(self):
        seq = easel.DigitalSequence(self.abc)
        self.assertEqual(seq.name, b"")
        seq.name = b"OTHER"
        self.assertEqual(seq.name, b"OTHER")

    def test_setter_description(self):
        seq = easel.DigitalSequence(self.abc)
        self.assertIs(seq.description, b"")
        seq.description = b"a test sequence"
        self.assertEqual(seq.description, b"a test sequence")

    def test_copy(self):
        arr = bytearray([0, 1, 2, 3])
        seq = easel.DigitalSequence(self.abc, name=b"TEST", sequence=arr)

        # copy the sequence
        cpy = seq.copy()

        # check the copy is equal and attributes are equal
        self.assertEqual(seq, cpy)
        self.assertEqual(seq.name, cpy.name)
        self.assertEqual(seq.checksum(), cpy.checksum())

        # check attributes can be read even after the original object
        # is (hopefully) deallocated, to make sure internal strings were copied
        del seq
        gc.collect()
        self.assertEqual(cpy.name, b"TEST")

        # check `copy.copy` works too
        cpy2 = copy.copy(cpy)
        self.assertEqual(cpy, cpy2)
        self.assertEqual(cpy.name, cpy2.name)

    def test_eq(self):
        arr = bytearray([0, 1, 2, 3])
        seq1 = easel.DigitalSequence(self.abc, name=b"TEST", sequence=arr)
        self.assertEqual(seq1, seq1)
        self.assertNotEqual(seq1, object())

        seq2 = easel.DigitalSequence(self.abc, name=b"TEST", sequence=arr)
        self.assertEqual(seq1, seq2)

        seq3 = easel.DigitalSequence(self.abc, name=b"OTHER", sequence=arr)
        self.assertNotEqual(seq1, seq3)

        arr2 = bytearray([0, 1, 2, 0])
        seq4 = easel.DigitalSequence(self.abc, name=b"TEST", sequence=arr2)
        self.assertNotEqual(seq1, seq4)

    def test_digitize_roundtrip(self):
        arr = bytearray([0, 1, 2, 3])
        dsq1 = easel.DigitalSequence(self.abc, name=b"TEST", sequence=arr)
        tsq1 = dsq1.textize()
        self.assertEqual(tsq1.name, dsq1.name)

        dsq2 = tsq1.digitize(self.abc)
        self.assertEqual(dsq1.name, dsq2.name)
        self.assertEqual(list(dsq1.sequence), list(dsq2.sequence))
        self.assertEqual(dsq1, dsq2)


class TestTextSequence(unittest.TestCase):

    def test_init_empty(self):
        seq = easel.TextSequence()
        self.assertEqual(seq.name, b"")
        self.assertEqual(seq.sequence, "")
        self.assertEqual(len(seq), 0)

    def test_init_kwargs(self):
        seq = easel.TextSequence(name=b"TEST", sequence="ATGC")
        self.assertEqual(seq.name, b"TEST")
        self.assertEqual(seq.sequence, "ATGC")
        self.assertEqual(len(seq), 4)

    def test_setter_name(self):
        seq = easel.TextSequence()
        self.assertEqual(seq.name, b"")
        seq.name = b"OTHER"
        self.assertEqual(seq.name, b"OTHER")

    def test_setter_description(self):
        seq = easel.TextSequence()
        self.assertIs(seq.description, b"")
        seq.description = b"a test sequence"
        self.assertEqual(seq.description, b"a test sequence")

    def test_copy(self):
        seq = easel.TextSequence(name=b"TEST", sequence="ATGC")

        # copy the sequence
        cpy = seq.copy()

        # check the copy is equal and attributes are equal
        self.assertEqual(seq, cpy)
        self.assertEqual(seq.name, cpy.name)
        self.assertEqual(seq.checksum(), cpy.checksum())

        # check attributes can be read even after the original object
        # is (hopefully) deallocated, to make sure internal strings were copied
        del seq
        gc.collect()
        self.assertEqual(cpy.name, b"TEST")

        # check `copy.copy` works too
        cpy2 = copy.copy(cpy)
        self.assertEqual(cpy, cpy2)
        self.assertEqual(cpy.name, cpy2.name)

    def test_eq(self):
        seq1 = easel.TextSequence(name=b"TEST", sequence="ATGC")
        self.assertEqual(seq1, seq1)
        self.assertNotEqual(seq1, object())

        seq2 = easel.TextSequence(name=b"TEST", sequence="ATGC")
        self.assertEqual(seq1, seq2)

        seq3 = easel.TextSequence(name=b"OTHER", sequence="ATGC")
        self.assertNotEqual(seq1, seq3)

        seq4 = easel.TextSequence(name=b"TEST", sequence="ATGT")
        self.assertNotEqual(seq1, seq4)

    def test_digitize_roundtrip(self):
        seq1 = easel.TextSequence(name=b"TEST", sequence="ATGC")
        dsq1 = seq1.digitize(easel.Alphabet.dna())
        self.assertEqual(seq1.name, dsq1.name)

        seq2 = dsq1.textize()
        self.assertEqual(seq1.name, seq2.name)
        self.assertEqual(seq1, seq2)


class TestSSIReader(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.h3i = os.path.realpath(os.path.join(
            __file__, os.pardir, "data", "hmms", "db", "Thioesterase.hmm.h3i"
        ))

    def test_init_error_filenotfound(self):
       self.assertRaises(FileNotFoundError, easel.SSIReader, "path/to/missing/file.ssi")

    def test_init_error_wrongformat(self):
       self.assertRaises(ValueError, easel.SSIReader, __file__)

    def test_find_name(self):
        with easel.SSIReader(self.h3i) as reader:
            entry = reader.find_name(b"Thioesterase")
            self.assertEqual(entry.fd, 0)
            self.assertEqual(entry.record_offset, 0)

    def test_find_name_closed(self):
        with easel.SSIReader(self.h3i) as reader:
            reader.close()
            self.assertRaises(ValueError, reader.find_name, b"Thioesterase")

    def test_find_name_error(self):
        with easel.SSIReader(self.h3i) as reader:
            self.assertRaises(KeyError, reader.find_name, b"Ketosynthase")

    def test_file_info(self):
        with easel.SSIReader(self.h3i) as reader:
            info = reader.file_info(0)
            self.assertEqual(info.format, 0)
            self.assertEqual(info.name, "Thioesterase.hmm")

    def test_file_info_error(self):
        with easel.SSIReader(self.h3i) as reader:
            self.assertRaises(IndexError, reader.file_info, 1)

    def test_file_info_closed(self):
        with easel.SSIReader(self.h3i) as reader:
            reader.close()
            self.assertRaises(ValueError, reader.file_info, 0)


class TestSSIWriter(unittest.TestCase):

    def test_init_error_fileexists(self):
        with tempfile.NamedTemporaryFile() as tmp:
            self.assertRaises(FileExistsError, easel.SSIWriter, tmp.name, exclusive=True)
