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

from . import errors
from . import easel
from . import plan7

from .hmmer import hmmsearch


__author__ = "Martin Larralde <martin.larralde@embl.de>"
__license__ = "MIT"
__version__ = "0.1.0-a3"
__all__ = [errors.__name__, easel.__name__, plan7.__name__, hmmsearch.__name__]
if __doc__ is not None:
    __doc__ += """See Also:
    An online rendered version of the documentation for this version of the
    library on `Read The Docs <https://pyhmmer.readthedocs.io/en/v{}/>`_.

    """.format(__version__)

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
