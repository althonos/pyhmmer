# coding: utf-8
# cython: language_level=3, linetrace=True

# --- C imports --------------------------------------------------------------

from libc.stdint cimport uint8_t, uint32_t
from posix.types cimport off_t
from cpython.pythread cimport PyThread_type_lock

from libeasel import eslINFINITY
from libeasel.sq cimport ESL_SQ
from libeasel.sqio cimport ESL_SQFILE
from libhmmer.p7_alidisplay cimport P7_ALIDISPLAY
from libhmmer.p7_bg cimport P7_BG
from libhmmer.p7_builder cimport P7_BUILDER
from libhmmer.p7_domain cimport P7_DOMAIN
from libhmmer.p7_hit cimport P7_HIT
from libhmmer.p7_hmm cimport P7_HMM, p7_NOFFSETS, p7_NEVPARAM, p7_NCUTOFFS
from libhmmer.p7_hmmfile cimport P7_HMMFILE
from libhmmer.p7_pipeline cimport P7_PIPELINE, p7_pipemodes_e, p7_zsetby_e
from libhmmer.p7_profile cimport P7_PROFILE
from libhmmer.p7_scoredata cimport P7_SCOREDATA
from libhmmer.p7_tophits cimport P7_TOPHITS
from libhmmer.p7_trace cimport P7_TRACE
from libhmmer.nhmmer cimport ID_LENGTH_LIST

if HMMER_IMPL == "VMX":
    from libhmmer.impl_vmx.p7_omx cimport P7_OM_BLOCK
    from libhmmer.impl_vmx.p7_oprofile cimport P7_OPROFILE
elif HMMER_IMPL == "SSE":
    from libhmmer.impl_sse.p7_omx cimport P7_OM_BLOCK
    from libhmmer.impl_sse.p7_oprofile cimport P7_OPROFILE
elif HMMER_IMPL == "NEON":
    from libhmmer.impl_neon.p7_omx cimport P7_OM_BLOCK
    from libhmmer.impl_neon.p7_oprofile cimport P7_OPROFILE

from .easel cimport (
    Alphabet,
    DigitalSequence,
    DigitalSequenceBlock,
    DigitalMSA,
    KeyHash,
    MSA,
    Randomness,
    VectorF,
    VectorU8,
    SequenceFile,
)


# --- Fused types ------------------------------------------------------------

ctypedef fused ScanTargets:
    HMMPressedFile
    OptimizedProfileBlock

ctypedef fused SearchTargets:
    SequenceFile
    DigitalSequenceBlock


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
    cdef public   str    score_matrix
    cdef public   double popen
    cdef public   double pextend
    cdef readonly Alphabet alphabet
    cdef readonly Randomness randomness

    cdef uint32_t    _seed    # the seed passed at builder initialization
    cdef P7_BUILDER* _bld

    cpdef tuple build(self, DigitalSequence sequence, Background background)
    cpdef tuple build_msa(self, DigitalMSA msa, Background background)
    cpdef Builder copy(self)


cdef class Cutoffs:
    cdef object              _owner
    cdef int*                _flags
    cdef bint                _is_profile
    cdef float[p7_NCUTOFFS]* _cutoffs

    cpdef VectorF as_vector(self)
    cpdef bint gathering_available(self)
    cpdef bint trusted_available(self)
    cpdef bint noise_available(self)


cdef class Domain:
    cdef readonly Alignment alignment
    cdef readonly Hit hit
    cdef P7_DOMAIN* _dom

cdef class Domains:
    cdef readonly Hit hit


cdef class EvalueParameters:
    cdef object              _owner
    cdef float[p7_NEVPARAM]* _evparams

    cpdef VectorF as_vector(self)


cdef class Hit:
    # a reference to the TopHits that owns the wrapped P7_HIT, kept so that
    # the internal data is never deallocated before the Python class.
    cdef readonly TopHits hits
    cdef P7_HIT* _hit


cdef class HMM:
    cdef P7_HMM* _hmm
    # a reference to the Alphabet Python object to avoid deallocation of the
    # inner ESL_ALPHABET; the Python object provides reference counting for free
    cdef readonly Alphabet          alphabet

    cdef void _initialize(self) noexcept nogil

    cpdef HMM copy(self)
    cpdef VectorF match_occupancy(self)
    cpdef double mean_match_entropy(self) except *
    cpdef double mean_match_information(self, Background background) except *
    cpdef double mean_match_relative_entropy(self, Background background) except *
    cpdef void renormalize(self)
    cpdef void scale(self, double scale, bint exponential=?)
    cpdef void set_composition(self) except *
    cpdef void set_consensus(self, DigitalSequence sequence=?) except *
    cpdef Profile to_profile(self,  Background background=?,  int L=*,  bint multihit =*,  bint local=*)
    cpdef void validate(self, float tolerance=*) except *
    cpdef void write(self, object fh, bint binary=*) except *
    cpdef void zero(self) noexcept


cdef class HMMFile:
    cdef str         _name
    cdef P7_HMMFILE* _hfp
    cdef Alphabet    _alphabet
    cdef object      _file

    cpdef void close(self) except *
    cpdef void rewind(self) except *
    cpdef HMM read(self)
    cpdef bint is_pressed(self) except *
    cpdef HMMPressedFile optimized_profiles(self)

    @staticmethod
    cdef P7_HMMFILE* _open_fileobj(object fh) except *


cdef class HMMPressedFile:
    cdef P7_HMMFILE* _hfp
    cdef Alphabet    _alphabet
    cdef HMMFile     _hmmfile

    cpdef void close(self) except *
    cpdef void rewind(self) except *
    cpdef OptimizedProfile read(self)


cdef class IterationResult:
    cdef readonly TopHits    hits
    cdef readonly DigitalMSA msa
    cdef readonly HMM        hmm
    cdef readonly bint       converged
    cdef readonly size_t     iteration


cdef class IterativeSearch:
    cdef readonly object               query
    cdef readonly Pipeline             pipeline
    cdef readonly Background           background
    cdef readonly Builder              builder
    cdef readonly bint                 converged
    cdef readonly DigitalSequenceBlock targets
    cdef readonly KeyHash              ranking
    cdef readonly size_t               iteration
    cdef          DigitalMSA           msa
    cdef          object               select_hits

    cpdef TopHits _search_hmm(self, HMM hmm)


cdef class OptimizedProfile:
    cdef P7_OPROFILE* _om
    cdef readonly Alphabet alphabet

    cpdef OptimizedProfile copy(self)
    cpdef void convert(self, Profile profile) except *
    cpdef void write(self, object fh_filter, object fh_profile) except *
    cpdef object ssv_filter(self, DigitalSequence seq)


cdef class OptimizedProfileBlock:
    cdef          PyThread_type_lock* _locks
    cdef          P7_OM_BLOCK*        _block
    cdef          list                _storage
    cdef readonly Alphabet            alphabet

    cdef void _allocate(self, size_t n) except *

    cpdef void append(self, OptimizedProfile optimized_profile) except *
    cpdef void clear(self) except *
    cpdef void extend(self, object iterable) except *
    cpdef size_t index(self, OptimizedProfile optimized_profile, ssize_t start=*, ssize_t stop=*) except *
    cpdef void insert(self, ssize_t index, OptimizedProfile optimized_profile) except *
    cpdef OptimizedProfile pop(self, ssize_t index=*)
    cpdef void remove(self, OptimizedProfile optimized_profile) except *
    cpdef OptimizedProfileBlock copy(self)


cdef class Offsets:
    cdef object              _owner
    cdef off_t[p7_NOFFSETS]* _offs


cdef class Pipeline:
    cdef          object     _Z            # either `Z` as an int, or `None`
    cdef          object     _domZ         # either `domZ` as an int, or `None`
    cdef          uint32_t   _seed         # the seed passed at pipeline initialization
    cdef          dict       _cutoff_save  # a local save of the reporting parameters

    cdef readonly Alphabet         alphabet
    cdef readonly Background       background
    cdef readonly Profile          profile
    cdef readonly OptimizedProfile opt
    cdef readonly Randomness       randomness

    cdef OptimizedProfile _optimized
    cdef P7_PIPELINE* _pli

    cdef int _save_cutoff_parameters(self) except 1
    cdef int _restore_cutoff_parameters(self) except 1
    
    @staticmethod
    cdef int _missing_cutoffs(const P7_PIPELINE* pli, const P7_OPROFILE* om) except 1 nogil
    cdef P7_OPROFILE* _get_om_from_query(self, object query, int L = *) except NULL

    cpdef list    arguments(self)
    cpdef void    clear(self)

    cpdef TopHits search_hmm(
        self,
        object query,
        SearchTargets sequences
    )
    cpdef TopHits search_msa(
        self,
        DigitalMSA query,
        object sequences,
        Builder builder = ?
    )
    cpdef TopHits search_seq(
        self,
        DigitalSequence query,
        object sequences,
        Builder builder = ?
    )
    @staticmethod
    cdef  int  _search_loop(
              P7_PIPELINE* pli,
              P7_OPROFILE* om,
              P7_BG*       bg,
        const ESL_SQ**     sq,
        const size_t       n_targets,
              P7_TOPHITS*  th,
    ) except 1 nogil
    @staticmethod
    cdef  int  _search_loop_file(
              P7_PIPELINE* pli,
              P7_OPROFILE* om,
              P7_BG*       bg,
              ESL_SQFILE*  sqfp,
              P7_TOPHITS*  th,
    ) except 1 nogil

    cpdef TopHits scan_seq(
        self,
        DigitalSequence query,
        ScanTargets targets
    )
    @staticmethod
    cdef int _scan_loop(
              P7_PIPELINE*        pli,
        const ESL_SQ*             sq,
              P7_BG*              bg,
              P7_OPROFILE**       om,
        const size_t              n_targets,
              P7_TOPHITS*         th,
              PyThread_type_lock* locks,
    ) except 1 nogil
    @staticmethod
    cdef int _scan_loop_file(
              P7_PIPELINE*     pli,
        const ESL_SQ*          sq,
              P7_BG*           bg,
              P7_HMMFILE*      hfp,
              P7_TOPHITS*      th,
    ) except 1 nogil

    cpdef IterativeSearch iterate_hmm(
        self,
        HMM query,
        DigitalSequenceBlock sequences,
        Builder builder = ?,
        object select_hits = ?,
    )
    cpdef IterativeSearch iterate_seq(
        self,
        DigitalSequence query,
        DigitalSequenceBlock sequences,
        Builder builder = ?,
        object select_hits = ?,
    )


cdef class LongTargetsPipeline(Pipeline):
    cdef int             _window_length
    cdef double          _window_beta
    cdef ID_LENGTH_LIST* _idlens

    cpdef TopHits search_hmm(
        self,
        object query,
        SearchTargets sequences
    )
    cpdef TopHits search_msa(
        self,
        DigitalMSA query,
        object sequences,
        Builder builder = ?
    )
    cpdef TopHits search_seq(
        self,
        DigitalSequence query,
        object sequences,
        Builder builder = ?
    )
    @staticmethod
    cdef int _search_loop_longtargets(
              P7_PIPELINE*    pli,
              P7_OPROFILE*    om,
              P7_BG*          bg,
        const ESL_SQ**        sq,
        const size_t          n_targets,
              P7_TOPHITS*     th,
              P7_SCOREDATA*   scoredata,
              ID_LENGTH_LIST* idlens
    ) except 1 nogil
    @staticmethod
    cdef int _search_loop_longtargets_file(
        P7_PIPELINE*  pli,
        P7_OPROFILE*  om,
        P7_BG*        bg,
        ESL_SQFILE*   sqfp,
        P7_TOPHITS*   th,
        P7_SCOREDATA* scoredata,
        ID_LENGTH_LIST* idlens
    ) except 1 nogil
    cpdef TopHits scan_seq(
        self,
        DigitalSequence query,
        ScanTargets targets
    )

cdef class Profile:
    cdef          P7_PROFILE* _gm
    cdef readonly Alphabet    alphabet

    cpdef void clear(self) except *
    cpdef void configure(self, HMM hmm, Background background, int L=?, bint multihit=*, bint local=*) except *
    cpdef Profile copy(self)
    cpdef OptimizedProfile to_optimized(self)


cdef class ScoreData:
    cdef P7_SCOREDATA* _sd
    cdef readonly int  Kp

    cpdef ScoreData copy(self)


cdef class TopHits:

    # NOTE(@althonos): this is not a full pipeline, but a local copy of the
    #                  accounting parameters so that the e-value can be
    #                  computed and thresholding can be done correctly.
    cdef P7_PIPELINE _pli
    cdef P7_TOPHITS* _th
    cdef bytes       _qname
    cdef bytes       _qacc
    cdef int         _qlen

    cdef int _threshold(self, Pipeline pipeline) except 1 nogil
    cdef int _sort_by_key(self) except 1 nogil
    cdef int _sort_by_seqidx(self) except 1 nogil
    cdef void _check_threshold_parameters(self, const P7_PIPELINE* other) except *

    cpdef TopHits copy(self)
    cpdef int compare_ranking(self, KeyHash) except -1
    cpdef bint is_sorted(self, str by=*) except *
    cpdef void sort(self, str by=*) except *
    cpdef MSA to_msa(self, Alphabet alphabet, list sequences=?, list traces=?, bint trim=*, bint digitize=?, bint all_consensus_cols=?)
    cpdef void write(self, object fh, str format=*, bint header=*) except *


cdef class Trace:
    cdef readonly Traces traces
    cdef P7_TRACE* _tr

    cpdef float expected_accuracy(self)
    cpdef float score(self, DigitalSequence sequence, Profile profile) except *


cdef class Traces:
    cdef P7_TRACE** _traces
    cdef size_t     _ntraces


cdef class TraceAligner:
    cdef ESL_SQ** _seqs
    cdef size_t   _nseq

    cpdef Traces compute_traces(self, HMM hmm, DigitalSequenceBlock sequences)
    cpdef MSA align_traces(
        self,
        HMM hmm,
        DigitalSequenceBlock sequences,
        Traces traces,
        bint digitize=*,
        bint trim=*,
        bint all_consensus_cols=*
    )
