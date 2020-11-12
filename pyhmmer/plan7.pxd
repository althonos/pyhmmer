# coding: utf-8
# cython: language_level=3, linetrace=True

# --- C imports --------------------------------------------------------------

from libhmmer.p7_alidisplay cimport P7_ALIDISPLAY
from libhmmer.p7_bg cimport P7_BG
from libhmmer.p7_domain cimport P7_DOMAIN
from libhmmer.p7_hmm cimport P7_HMM
from libhmmer.p7_hmmfile cimport P7_HMMFILE
from libhmmer.p7_pipeline cimport P7_PIPELINE
from libhmmer.p7_profile cimport P7_PROFILE
from libhmmer.p7_tophits cimport P7_TOPHITS, P7_HIT

from pyhmmer.easel cimport Alphabet

# --- Cython classes ---------------------------------------------------------


cdef class Alignment:
    cdef readonly Domain domain
    cdef P7_ALIDISPLAY* _ad


cdef class Domain:
    cdef readonly Alignment alignment
    cdef readonly Hit hit
    cdef P7_DOMAIN* _dom


cdef class Domains:
    cdef readonly Hit hit


cdef class Hit:
    # a reference to the TopHits that owns the wrapped P7_HIT, kept so that
    # the internal data is never deallocated before the Python class.
    cdef readonly TopHits hits
    cdef P7_HIT* _hit


cdef class HMM:
    # a reference to the Alphabet Python object to avoid deallocation of the
    # inner ESL_ALPHABET; the Python object provides reference counting for free
    cdef readonly Alphabet alphabet
    cdef P7_HMM* _hmm

    cpdef void zero(self)


cdef class HMMFile:
    cdef P7_HMMFILE* _hfp
    cdef Alphabet _alphabet

    cpdef void close(self)


cdef class Profile:
    cdef readonly Alphabet alphabet
    cdef P7_PROFILE* _gm

    cpdef void clear(self)
    cpdef Profile copy(self)
    cpdef bint is_local(self)
    cpdef bint is_multihit(self)


cdef class TopHits:
    cdef P7_TOPHITS* _th

    cpdef void threshold(self, Pipeline pipeline)
    cpdef void clear(self)
    cpdef void sort(self, str by=*)
    cpdef bint is_sorted(self, str by=*)


cdef class Pipeline:
    cdef readonly Alphabet alphabet
    cdef P7_PIPELINE* _pli
    cdef P7_BG* _bg

    cpdef TopHits search(self, HMM hmm, object seqs, TopHits hits=?)
