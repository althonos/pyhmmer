# coding: utf-8
"""Cython bindings and Python interface to HMMER3.

HMMER is a biological sequence analysis tool that uses profile hidden Markov
models to search for sequence homologs. HMMER3 is maintained by members of the
the `Eddy/Rivas Laboratory <http://eddylab.org/>`_ at Harvard University.

``pyhmmer`` is a module, implemented using the `Cython <https://cython.org/>`_
language, that provides bindings to HMMER3. It directly interacts with the
HMMER internals, which has several advantages over CLI wrappers like
`hmmer-py <https://pypi.org/project/hmmer/>`_.

"""

import collections.abc as _collections_abc
import contextlib as _contextlib
import os as _os

from . import errors
from . import easel
from . import plan7
from . import daemon

from .hmmer import hmmalign, hmmsearch, hmmpress, nhmmer, hmmscan, phmmer


__author__ = "Martin Larralde <martin.larralde@embl.de>"
__license__ = "MIT"
__version__ = "0.7.0-rc3"
__all__ = [
    "errors",
    "easel",
    "plan7",
    "daemon",
    "hmmalign",
    "hmmsearch",
    "hmmscan",
    "hmmpress",
    "phmmer",
    "nhmmer",
]

# Small addition to the docstring: we want to show a link redirecting to the
# rendered version of the documentation, but this can only work when Python
# is running with docstrings enabled
if __doc__ is not None:
    __doc__ += """See Also:
    An online rendered version of the documentation for this version of the
    library on `Read The Docs <https://pyhmmer.readthedocs.io/en/v{}/>`_.

    """.format(__version__)

# Register collections using the `collections.abc` module (this is probably
# not required with later versions of Python)
_collections_abc.Iterator.register(easel.SequenceFile)
_collections_abc.Iterator.register(plan7.HMMFile)
_collections_abc.Mapping.register(easel.KeyHash)
_collections_abc.Sized.register(plan7.Alignment)
_collections_abc.Sequence.register(easel.Bitfield)
_collections_abc.Sequence.register(plan7.Domains)
_collections_abc.Sequence.register(plan7.TopHits)
_collections_abc.Sequence.register(easel.SequenceBlock)
_collections_abc.Sequence.register(plan7.OptimizedProfileBlock)

if hasattr(_contextlib, "AbstractContextManager"):
    _contextlib.AbstractContextManager.register(easel.SequenceFile)
    _contextlib.AbstractContextManager.register(easel.SSIReader)
    _contextlib.AbstractContextManager.register(easel.SSIWriter)
    _contextlib.AbstractContextManager.register(plan7.HMMFile)
