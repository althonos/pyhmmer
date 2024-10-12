import copy
import gc
import io
import os
import unittest
import tempfile

from pyhmmer import easel

from .. import __name__ as __package__
from .utils import resource_files


@unittest.skipUnless(resource_files, "importlib.resources.files not available")
class TestSSIReader(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.h3i = resource_files(__package__).joinpath("data", "hmms", "db", "Thioesterase.hmm.h3i")

    def test_init_error_filenotfound(self):
        self.assertRaises(
            FileNotFoundError, easel.SSIReader, "path/to/missing/file.ssi"
        )

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
            self.assertRaises(
                FileExistsError, easel.SSIWriter, tmp.name, exclusive=True
            )
