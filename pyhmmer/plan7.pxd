# coding: utf-8
# cython: language_level=3, linetrace=True

# --- C imports --------------------------------------------------------------

from libc.stdint cimport uint32_t
from posix.types cimport off_t

from libhmmer.p7_alidisplay cimport P7_ALIDISPLAY
from libhmmer.p7_bg cimport P7_BG
from libhmmer.p7_domain cimport P7_DOMAIN
from libhmmer.p7_hmm cimport P7_HMM
from libhmmer.p7_hmmfile cimport P7_HMMFILE
from libhmmer.p7_pipeline cimport P7_PIPELINE, p7_pipemodes_e
from libhmmer.p7_profile cimport P7_PROFILE
from libhmmer.p7_tophits cimport P7_TOPHITS, P7_HIT

IF HMMER_IMPL == "VMX":
    from libhmmer.impl_vmx.p7_oprofile cimport P7_OPROFILE
ELIF HMMER_IMPL == "SSE":
    from libhmmer.impl_sse.p7_oprofile cimport P7_OPROFILE

from .easel cimport Alphabet


cdef extern from "hmmer.h" nogil:
    DEF p7_NOFFSETS = 3

# --- Cython classes ---------------------------------------------------------


cdef class Alignment:
    cdef readonly Domain domain
    cdef P7_ALIDISPLAY* _ad


cdef class Background:
    cdef readonly Alphabet alphabet
    cdef P7_BG* _bg

    cpdef Background copy(self)


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

    cpdef bint is_included(self)
    cpdef bint is_reported(self)
    cpdef bint is_new(self)
    cpdef bint is_dropped(self)
    cpdef bint is_duplicate(self)


cdef class HMM:
    # a reference to the Alphabet Python object to avoid deallocation of the
    # inner ESL_ALPHABET; the Python object provides reference counting for free
    cdef readonly Alphabet alphabet
    cdef P7_HMM* _hmm

    cpdef void write(self, object fh, bint binary=*) except *
    cpdef void zero(self)


cdef class HMMFile:
    cdef P7_HMMFILE* _hfp
    cdef Alphabet _alphabet

    cpdef void close(self)

    cdef P7_HMMFILE* _open_fileobj(self, object fh) except *


cdef class OptimizedProfile:
    cdef readonly Alphabet alphabet
    cdef P7_OPROFILE* _om

    cpdef OptimizedProfile copy(self)
    cpdef bint is_local(self)
    cpdef void write(self, object fh_filter, object fh_profile)


cdef class _Offsets:
    cdef OptimizedProfile    opt
    cdef off_t[p7_NOFFSETS]* _offs


cdef class Pipeline:
    cdef public   uint32_t   seed
    cdef public   bint       null2
    cdef public   bint       bias_filter
    cdef public   float      report_e
    cdef          object     _Z
    cdef          object     _domZ

    cdef readonly Alphabet   alphabet
    cdef readonly Background background
    cdef readonly Profile    profile

    cdef OptimizedProfile _optimized
    cdef P7_PIPELINE* _pli

    cpdef void clear(self)
    cpdef TopHits search(self, HMM hmm, object seqs, TopHits hits=?)


cdef class Profile:
    cdef readonly Alphabet alphabet
    cdef P7_PROFILE* _gm

    cpdef void clear(self)
    cpdef void configure(self, HMM hmm, Background background, int L, bint multihit=*, bint local=*)
    cpdef Profile copy(self)
    cpdef bint is_local(self)
    cpdef bint is_multihit(self)
    cpdef OptimizedProfile optimized(self)


cdef class TopHits:
    cdef public float Z
    cdef public float domZ
    cdef public bint  long_targets

    cdef P7_TOPHITS* _th

    cdef void threshold(self, Pipeline pipeline)
    cpdef void clear(self)
    cpdef void sort(self, str by=*)
    cpdef bint is_sorted(self, str by=*)
