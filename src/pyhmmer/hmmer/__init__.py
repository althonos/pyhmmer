# coding: utf-8
"""Reimplementation of HMMER binaries with the PyHMMER API.

Note:
    Functions of this module handle parallelization using threads to run
    searches in parallel for the different queries. If less queries are
    given, the number of threads will be reduced to avoid spawning idle
    threads.

"""

from ._hmmalign import hmmalign
from ._hmmpress import hmmpress
from ._hmmscan import hmmscan
from ._hmmsearch import hmmsearch
from ._jackhmmer import jackhmmer
from ._nhmmer import nhmmer
from ._phmmer import phmmer

__all__ = [
    "hmmalign",
    "hmmpress",
    "hmmscan",
    "hmmsearch",
    "jackhmmer",
    "nhmmer",
    "phmmer",
]