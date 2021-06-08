# coding: utf-8
# cython: language_level=3, linetrace=True

# --- C imports --------------------------------------------------------------

from libc.stdint cimport uint32_t
from posix.types cimport off_t

from libeasel.sq cimport ESL_SQ
from libhmmer.p7_alidisplay cimport P7_ALIDISPLAY
from libhmmer.p7_bg cimport P7_BG
from libhmmer.p7_builder cimport P7_BUILDER
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

from .easel cimport Alphabet, DigitalSequence, DigitalMSA, MSA, Randomness, VectorF


cdef extern from "hmmer.h" nogil:
    DEF p7_NOFFSETS = 3

# --- Cython classes ---------------------------------------------------------


cdef class Alignment:
    cdef readonly Domain domain
    cdef P7_ALIDISPLAY* _ad


cdef class Background:
    cdef readonly bint     uniform
    cdef readonly Alphabet alphabet
    cdef readonly VectorF  residue_frequencies
    cdef          P7_BG*   _bg
    cdef          int      _L

    cpdef Background copy(self)


cdef class Builder:
    cdef readonly object prior_scheme
    cdef readonly object effective_number
    cdef readonly object architecture
    cdef readonly object weighting
    cdef readonly float  popen
    cdef readonly float  pextend
    cdef readonly Alphabet alphabet
    cdef readonly Randomness randomness

    cdef uint32_t    _seed    # the seed passed at builder initialization
    cdef P7_BUILDER* _bld

    cpdef tuple build(self, DigitalSequence sequence, Background background)
    cpdef tuple build_msa(self, DigitalMSA msa, Background background)
    cpdef Builder copy(self)


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

    cpdef HMM copy(self)
    cpdef void write(self, object fh, bint binary=*) except *
    cpdef void zero(self)
    cpdef void renormalize(self)
    cpdef void scale(self, double scale, bint exponential=?)
    cpdef void set_composition(self)


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

    cpdef object ssv_filter(self, DigitalSequence seq)


cdef class _Offsets:
    cdef object              _owner
    cdef off_t[p7_NOFFSETS]* _offs


cdef class Pipeline:
    cdef          object     _Z       # either `Z` as an int, or `None`
    cdef          object     _domZ    # either `domZ` as an int, or `None`
    cdef          uint32_t   _seed    # the seed passed at pipeline initialization

    cdef readonly Alphabet   alphabet
    cdef readonly Background background
    cdef readonly Profile    profile
    cdef readonly Randomness randomness

    cdef OptimizedProfile _optimized
    cdef P7_PIPELINE* _pli

    cpdef void    clear(self)
    cpdef TopHits search_hmm(self, HMM query, object seqs)
    cpdef TopHits search_msa(self, DigitalMSA query, object seqs, Builder builder = ?)
    cpdef TopHits search_seq(self, DigitalSequence query, object seqs, Builder builder = ?)
    cdef  int    _search_loop(
        self,
        P7_PIPELINE* pli,
        P7_OPROFILE* om,
        P7_BG*       bg,
        ESL_SQ*      sq,
        P7_TOPHITS*  th,
        object       seqs_iter,
        Alphabet     seq_alphabet,
    ) except 1
    cpdef TopHits scan_seq(self, DigitalSequence query, object hmms)
    cdef int _scan_loop(
        self,
        P7_PIPELINE* pli,
        ESL_SQ*      sq,
        P7_BG*       bg,
        P7_HMM*      hm,
        P7_TOPHITS*  th,
        object       hmm_iter,
        Alphabet     hmm_alphabet
    ) except 1


cdef class Profile:
    cdef readonly Alphabet alphabet
    cdef P7_PROFILE* _gm

    cdef int _clear(self) except 1
    cdef int _configure(self, HMM hmm, Background background, int L, bint multihit=*, bint local=*) except 1
    cpdef Profile copy(self)
    cpdef bint is_local(self)
    cpdef bint is_multihit(self)
    cpdef OptimizedProfile optimized(self)


cdef class TopHits:
    cdef public float Z
    cdef public float domZ
    cdef public bint  long_targets

    cdef P7_TOPHITS* _th

    cdef int _threshold(self, Pipeline pipeline) except 1
    cdef int _sort(self, str by=*) except 1
    cpdef MSA to_msa(self, Alphabet alphabet, bint trim=*, bint digitize=?, bint all_consensus_cols=?)
