"""Cython bindings and Python interface to HMMER3.
"""

import collections.abc as _collections_abc
import contextlib as _contextlib

from . import errors
from . import easel
from . import plan7

from .hmmer import hmmsearch


__author__ = "Martin Larralde <martin.larralde@embl.de>"
__license__ = "MIT"
__version__ = "0.1.0-a2"
__all__ = [errors.__name__, easel.__name__, plan7.__name__, hmmsearch.__name__]


# Register collections using the `collections.abc` module (this is probably
# not required with later versions of Python)
_collections_abc.Iterator.register(easel.SequenceFile)
_collections_abc.Iterator.register(plan7.HMMFile)
_collections_abc.Sized.register(plan7.Alignment)
_collections_abc.Sequence.register(easel.Bitfield)
_collections_abc.Sequence.register(plan7.Domains)
_collections_abc.Sequence.register(plan7.TopHits)

if hasattr(_contextlib, "AbstractContextManager"):
    _contextlib.AbstractContextManager.register(easel.SequenceFile)
    _contextlib.AbstractContextManager.register(plan7.HMMFile)
