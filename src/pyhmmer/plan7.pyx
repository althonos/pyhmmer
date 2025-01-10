# coding: utf-8
# cython: language_level=3
"""High-level interface to the Plan7 data model.

Plan7 is the model architecture used by HMMER since HMMER2.

See Also:
    Details about the Plan 7 architecture in the `HMMER documentation
    <http://www.csb.yale.edu/userguides/seq/hmmer/docs/node11.html>`_.

"""

# --- C imports --------------------------------------------------------------

cimport cython
from cpython.bytes cimport PyBytes_FromStringAndSize, PyBytes_FromString
from cpython.list cimport PyList_New, PyList_SET_ITEM
from cpython.ref cimport PyObject
from cpython.exc cimport PyErr_Clear
from cpython.unicode cimport PyUnicode_DecodeASCII
from libc.math cimport exp, ceil
from libc.stddef cimport ptrdiff_t
from libc.stdio cimport printf, rewind
from libc.stdlib cimport calloc, malloc, realloc, free, llabs
from libc.stdint cimport uint8_t, uint32_t, uint64_t, int64_t
from libc.stdio cimport fprintf, FILE, stdout, fclose
from libc.string cimport memset, memcpy, memmove, strdup, strndup, strlen, strcmp, strncpy
from libc.time cimport ctime, strftime, time, time_t, tm, localtime_r
from cpython.unicode cimport (
    PyUnicode_DATA,
    PyUnicode_KIND,
    PyUnicode_READ,
    PyUnicode_GET_LENGTH,
    PyUnicode_FromStringAndSize,
)
from cpython.pythread cimport (
    PyThread_type_lock,
    PyThread_allocate_lock,
    PyThread_free_lock,
    PyThread_acquire_lock,
    PyThread_release_lock,
    PyLockStatus,
    WAIT_LOCK,
)

cimport libeasel
cimport libeasel.sq
cimport libeasel.sqio
cimport libeasel.alphabet
cimport libeasel.dmatrix
cimport libeasel.fileparser
cimport libeasel.random
cimport libeasel.scorematrix
cimport libeasel.getopts
cimport libeasel.vec
cimport libhmmer
cimport libhmmer.modelconfig
cimport libhmmer.modelstats
cimport libhmmer.p7_alidisplay
cimport libhmmer.p7_hmm
cimport libhmmer.p7_builder
cimport libhmmer.p7_bg
cimport libhmmer.p7_domain
cimport libhmmer.p7_domaindef
cimport libhmmer.p7_hit
cimport libhmmer.p7_hmmfile
cimport libhmmer.p7_pipeline
cimport libhmmer.p7_prior
cimport libhmmer.p7_profile
cimport libhmmer.p7_scoredata
cimport libhmmer.p7_tophits
cimport libhmmer.p7_trace
cimport libhmmer.tracealign
cimport libhmmer.nhmmer
from libeasel cimport eslERRBUFSIZE, eslCONST_LOG2R, eslINFINITY
from libeasel.alphabet cimport ESL_ALPHABET, esl_alphabet_Create, esl_abc_ValidateType
from libeasel.getopts cimport ESL_GETOPTS, ESL_OPTIONS
from libeasel.sq cimport ESL_SQ
from libeasel.keyhash cimport ESL_KEYHASH
from libeasel.fileparser cimport ESL_FILEPARSER
from libhmmer cimport (
    p7_MAXABET,
    p7_LOCAL,
    p7_EVPARAM_UNSET,
    p7_CUTOFF_UNSET,
    p7_NEVPARAM,
    p7_NCUTOFFS,
    p7_offsets_e,
    p7_cutoffs_e,
    p7_evparams_e,
)
from libhmmer.logsum cimport p7_FLogsumInit
from libhmmer.p7_builder cimport P7_BUILDER, p7_archchoice_e, p7_wgtchoice_e, p7_effnchoice_e
from libhmmer.p7_hmm cimport p7H_NTRANSITIONS, p7H_TC, p7H_GA, p7H_NC, p7H_MAP, p7h_transitions_e
from libhmmer.p7_hmmfile cimport p7_hmmfile_formats_e
from libhmmer.p7_hit cimport p7_hitflags_e, P7_HIT
from libhmmer.p7_alidisplay cimport P7_ALIDISPLAY
from libhmmer.p7_pipeline cimport P7_PIPELINE, p7_pipemodes_e, p7_zsetby_e, p7_strands_e, p7_complementarity_e
from libhmmer.p7_profile cimport p7_LOCAL, p7_GLOCAL, p7_UNILOCAL, p7_UNIGLOCAL
from libhmmer.p7_trace cimport P7_TRACE, p7t_statetype_e
from libhmmer.p7_prior cimport P7_PRIOR
from libhmmer.nhmmer cimport ID_LENGTH_LIST
from capacity cimport new_capacity

if HMMER_IMPL == "VMX":
    from libhmmer.impl_vmx cimport p7_oprofile, p7_omx, impl_Init, p7O_EXTRA_SB
    from libhmmer.impl_vmx.io cimport p7_oprofile_Write, p7_oprofile_ReadMSV, p7_oprofile_ReadRest
    from libhmmer.impl_vmx.p7_omx cimport (
        P7_OM_BLOCK,
        p7_oprofile_CreateBlock,
        p7_oprofile_DestroyBlock,
    )
    from libhmmer.impl_vmx.p7_oprofile cimport (
        P7_OPROFILE,
        p7O_NXSTATES,
        p7O_NXTRANS,
        p7O_NQB,
        p7O_NQF,
        p7_oprofile_Compare,
        p7_oprofile_Dump,
        p7_oprofile_Sizeof,
        p7_oprofile_Destroy,
    )
elif HMMER_IMPL == "SSE":
    from libhmmer.impl_sse cimport p7_oprofile, p7_omx, impl_Init, p7_SSVFilter, p7O_EXTRA_SB
    from libhmmer.impl_sse.io cimport p7_oprofile_Write, p7_oprofile_ReadMSV, p7_oprofile_ReadRest
    from libhmmer.impl_sse.p7_omx cimport (
        P7_OM_BLOCK,
        p7_oprofile_CreateBlock,
        p7_oprofile_DestroyBlock,
    )
    from libhmmer.impl_sse.p7_oprofile cimport (
        P7_OPROFILE,
        p7O_NXSTATES,
        p7O_NXTRANS,
        p7O_NQB,
        p7O_NQF,
        p7_oprofile_Compare,
        p7_oprofile_Dump,
        p7_oprofile_Sizeof,
        p7_oprofile_Destroy
    )
elif HMMER_IMPL == "NEON":
    from libhmmer.impl_neon cimport p7_oprofile, p7_omx, impl_Init, p7O_EXTRA_SB
    from libhmmer.impl_neon.io cimport p7_oprofile_Write, p7_oprofile_ReadMSV, p7_oprofile_ReadRest
    from libhmmer.impl_neon.p7_omx cimport (
        P7_OM_BLOCK,
        p7_oprofile_CreateBlock,
        p7_oprofile_DestroyBlock,
    )
    from libhmmer.impl_neon.p7_oprofile cimport (
        P7_OPROFILE,
        p7O_NXSTATES,
        p7O_NXTRANS,
        p7O_NQB,
        p7O_NQF,
        p7_oprofile_Compare,
        p7_oprofile_Dump,
        p7_oprofile_Sizeof,
        p7_oprofile_Destroy,
    )

from .easel cimport (
    Alphabet,
    Sequence,
    DigitalSequence,
    DigitalSequenceBlock,
    KeyHash,
    MSA,
    TextMSA,
    DigitalMSA,
    VectorF,
    VectorU8,
    MatrixF,
    MatrixU8,
    Randomness,
)
from .reexports.p7_tophits cimport p7_tophits_Reuse
from .reexports.p7_hmmfile cimport (
    read_asc20hmm,
    read_asc30hmm,
    read_bin30hmm,
    v3a_magic,
    v3b_magic,
    v3c_magic,
    v3d_magic,
    v3e_magic,
    v3f_magic
)

if TARGET_SYSTEM == "Linux":
    from .fileobj.linux cimport fileobj_linux_open as fopen_obj
elif TARGET_SYSTEM == "Darwin" or TARGET_SYSTEM.endswith("BSD"):
    from .fileobj.bsd cimport fileobj_bsd_open as fopen_obj

include "exceptions.pxi"
include "_getid.pxi"


# --- Python imports ---------------------------------------------------------

import array
import collections.abc
import copy
import datetime
import enum
import errno
import math
import io
import itertools
import operator
import os
import sys
import warnings

from .utils import peekable, SizedIterator
from .errors import (
    AllocationError,
    UnexpectedError,
    AlphabetMismatch,
    MissingCutoffs,
    InvalidParameter,
    InvalidHMM,
)


# --- Constants --------------------------------------------------------------

cdef dict BUILDER_ARCHITECTURE_STRATEGY = {
    "fast": p7_archchoice_e.p7_ARCH_FAST,
    "hand": p7_archchoice_e.p7_ARCH_HAND,
}

cdef dict BUILDER_WEIGHTING_STRATEGY = {
    "pb": p7_wgtchoice_e.p7_WGT_PB,
    "gsc": p7_wgtchoice_e.p7_WGT_GSC,
    "blosum": p7_wgtchoice_e.p7_WGT_BLOSUM,
    "none": p7_wgtchoice_e.p7_WGT_NONE,
    "given": p7_wgtchoice_e.p7_WGT_GIVEN,
}

cdef dict BUILDER_EFFECTIVE_STRATEGY = {
    "entropy": p7_effnchoice_e.p7_EFFN_ENTROPY,
    "exp": p7_effnchoice_e.p7_EFFN_ENTROPY_EXP,
    "clust": p7_effnchoice_e.p7_EFFN_CLUST,
    "none": p7_effnchoice_e.p7_EFFN_NONE,
}

cdef dict HMM_FILE_FORMATS = {
    "2.0": p7_hmmfile_formats_e.p7_HMMFILE_20,
    "3/a": p7_hmmfile_formats_e.p7_HMMFILE_3a,
    "3/b": p7_hmmfile_formats_e.p7_HMMFILE_3b,
    "3/c": p7_hmmfile_formats_e.p7_HMMFILE_3c,
    "3/d": p7_hmmfile_formats_e.p7_HMMFILE_3d,
    "3/e": p7_hmmfile_formats_e.p7_HMMFILE_3e,
    "3/f": p7_hmmfile_formats_e.p7_HMMFILE_3f,
}

cdef dict HMM_FILE_MAGIC = {
    v3a_magic: p7_hmmfile_formats_e.p7_HMMFILE_3a,
    v3b_magic: p7_hmmfile_formats_e.p7_HMMFILE_3b,
    v3c_magic: p7_hmmfile_formats_e.p7_HMMFILE_3c,
    v3d_magic: p7_hmmfile_formats_e.p7_HMMFILE_3d,
    v3e_magic: p7_hmmfile_formats_e.p7_HMMFILE_3e,
    v3f_magic: p7_hmmfile_formats_e.p7_HMMFILE_3f,
}

cdef dict PIPELINE_BIT_CUTOFFS = {
    "gathering": libhmmer.p7_hmm.p7H_GA,
    "noise": libhmmer.p7_hmm.p7H_NC,
    "trusted": libhmmer.p7_hmm.p7H_TC,
}

# --- Cython classes ---------------------------------------------------------

cdef class Alignment:
    """An alignment of a sequence to a profile.

    Attributes:
        domain (`Domain`): The domain this alignment corresponds to.

    .. versionadded:: 0.6.1
       `pickle` protocol support.

    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self, Domain domain):
        self.domain = domain
        self._ad = domain._dom.ad

    def __len__(self):
        return self.hmm_to - self.hmm_from

    def __str__(self):
        assert self._ad != NULL

        cdef int    status
        cdef object buffer = io.BytesIO()
        cdef FILE*  fp     = fopen_obj(buffer, "w")

        try:
            status = libhmmer.p7_alidisplay.p7_nontranslated_alidisplay_Print(
                fp,
                self._ad,
                0,
                -1,
                False,
            )
            if status == libeasel.eslEWRITE:
                raise OSError("Failed to write alignment")
            elif status != libeasel.eslOK:
                raise UnexpectedError(status, "p7_alidisplay_Print")
        finally:
            fclose(fp)

        return buffer.getvalue().decode("ascii")

    def __sizeof__(self):
        assert self._ad != NULL
        return sizeof(self) + libhmmer.p7_alidisplay.p7_alidisplay_Sizeof(self._ad)

    def __getstate__(self):
        cdef int       status
        cdef uint32_t  n      = 0
        cdef uint32_t  nalloc = 0
        cdef uint8_t*  buffer = NULL
        cdef VectorU8  vec

        with nogil:
            status = libhmmer.p7_alidisplay.p7_alidisplay_Serialize(self._ad, &buffer, &n, &nalloc)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_alidisplay_Serialize")

        vec = VectorU8.__new__(VectorU8)
        vec._owner = None
        vec._n = vec._shape[0] = nalloc
        vec._data = <void*> buffer
        return vec

    def __setstate__(self, uint8_t[::1] state):
        cdef int      status
        cdef uint32_t offset = 0

        if self._ad == NULL:
            self.domain._dom.ad = self._ad = libhmmer.p7_alidisplay.p7_alidisplay_Create_empty()
            if self._ad == NULL:
                raise AllocationError("P7_ALIDISPLAY", sizeof(P7_ALIDISPLAY))

        with nogil:
            status = libhmmer.p7_alidisplay.p7_alidisplay_Deserialize(&state[0], &offset, self._ad)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_alidisplay_Deserialize")

    # --- Properties ---------------------------------------------------------

    @property
    def hmm_from(self):
        """`int`: The start coordinate of the alignment in the query HMM.
        """
        assert self._ad != NULL
        return self._ad.hmmfrom

    @property
    def hmm_to(self):
        """`int`: The end coordinate of the alignment in the query HMM.
        """
        assert self._ad != NULL
        return self._ad.hmmto

    @property
    def hmm_name(self):
        """`bytes`: The name of the query HMM.
        """
        assert self._ad != NULL
        assert self._ad.hmmname != NULL
        return <bytes> self._ad.hmmname

    @property
    def hmm_accession(self):
        """`bytes`: The accession of the query, or its name if it has none.

        .. versionadded:: 0.1.4

        """
        assert self._ad != NULL
        assert self._ad.hmmacc != NULL
        return <bytes> self._ad.hmmacc

    @property
    def hmm_sequence(self):
        """`str`: The sequence of the query HMM in the alignment.
        """
        assert self._ad != NULL
        return self._ad.model.decode('ascii')

    @property
    def hmm_length(self):
        """`int`: The length of the query HMM in the alignment.

        .. versionadded:: 0.10.5

        """
        return self._ad.M

    @property
    def posterior_probabilities(self):
        """`str`: Posterior probability annotation of the alignment.

        .. versionadded:: 0.10.5

        """
        assert self._ad != NULL
        return self._ad.ppline.decode('ascii')


    @property
    def target_from(self):
        """`int`: The start coordinate of the alignment in the target sequence.
        """
        assert self._ad != NULL
        return self._ad.sqfrom

    @property
    def target_name(self):
        """`bytes`: The name of the target sequence.
        """
        assert self._ad != NULL
        assert self._ad.sqname != NULL
        return <bytes> self._ad.sqname

    @property
    def target_sequence(self):
        """`str`: The sequence of the target sequence in the alignment.
        """
        assert self._ad != NULL
        return self._ad.aseq.decode('ascii')

    @property
    def target_to(self):
        """`int`: The end coordinate of the alignment in the target sequence.
        """
        assert self._ad != NULL
        return self._ad.sqto

    @property
    def target_length(self):
        """`int`: The length of the target sequence in the alignment.

        .. versionadded:: 0.10.5

        """
        return self._ad.L

    @property
    def identity_sequence(self):
        """`str`: The identity sequence between the query and the target.
        """
        assert self._ad != NULL
        return self._ad.mline.decode('ascii')


cdef class Background:
    """The null background model of HMMER.

    Attributes:
        alphabet (`~pyhmmer.easel.Alphabet`): The alphabet of the
            background model.
        uniform (`bool`): Whether or not the null model has been created
            with uniform frequencies.
        residue_frequencies (`~pyhmmer.easel.VectorF`): The null1 background
            residue frequencies.

    .. versionchanged:: 0.4.0
       Added the `residue_frequencies` attribute.

    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._bg = NULL
        self._L  = 350
        self.alphabet = None
        self.residue_frequencies = None
        self.uniform = False

    def __init__(self, Alphabet alphabet, bint uniform=False):
        """__init__(self, alphabet, uniform=False)\n--\n

        Create a new background model for the given ``alphabet``.

        Arguments:
          alphabet (`~pyhmmer.easel.Alphabet`): The alphabet to create
              the background model with.
          uniform (`bool`): Whether or not to create the null model with
              uniform frequencies. Defaults to `False`.

        """
        # store the alphabet so it's not deallocated
        self.alphabet = alphabet
        # store whether or not the null model has uniform frequencies
        self.uniform = uniform
        # create the background profile
        with nogil:
            if uniform:
                self._bg = libhmmer.p7_bg.p7_bg_CreateUniform(alphabet._abc)
            else:
                self._bg = libhmmer.p7_bg.p7_bg_Create(alphabet._abc)
        if self._bg == NULL:
            raise AllocationError("P7_BG", sizeof(P7_BG))
        # expose the residue frequencies as the `residue_frequencies` attribute
        self.residue_frequencies = VectorF.__new__(VectorF)
        self.residue_frequencies._data = &(self._bg.f[0])
        self.residue_frequencies._owner = self
        self.residue_frequencies._n = self.residue_frequencies._shape[0] = self.alphabet.K

    def __dealloc__(self):
        libhmmer.p7_bg.p7_bg_Destroy(self._bg)

    def __copy__(self):
        return self.copy()

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"{ty}({self.alphabet!r}, uniform={self.uniform!r})"

    def __reduce__(self):
        return type(self), (self.alphabet, self.uniform)

    # --- Properties ---------------------------------------------------------

    @property
    def L(self):
        """`int`: The mean of the null model length distribution, in residues.
        """
        assert self._bg != NULL
        return self._L

    @L.setter
    def L(self, int L):
        assert self._bg != NULL

        cdef int    status
        cdef P7_BG* bg     = self._bg

        with nogil:
            status = libhmmer.p7_bg.p7_bg_SetLength(bg, L)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_bg_SetLength")
        self._L = L

    @property
    def transition_probability(self):
        r"""`float`: The null1 transition probability (:math:`\frac{L}{L+1}`).

        .. versionadded:: 0.4.0

        """
        assert self._bg != NULL
        return self._bg.p1

    @property
    def omega(self):
        """`float`: The *prior* on *null2*/*null3*.

        .. versionadded:: 0.4.0

        """
        assert self._bg != NULL
        return self._bg.omega

    @omega.setter
    def omega(self, float omega):
        assert self._bg != NULL
        self._bg.omega = omega

    # --- Methods ------------------------------------------------------------

    cpdef Background copy(self):
        """Create a copy of the null model with the same parameters.
        """
        assert self._bg != NULL

        cdef Background new = Background.__new__(Background)
        with nogil:
            new._bg = libhmmer.p7_bg.p7_bg_Clone(self._bg)
        if new._bg == NULL:
            raise AllocationError("P7_BG", sizeof(P7_BG))

        new.alphabet = self.alphabet
        new.uniform = self.uniform
        new.residue_frequencies = VectorF.__new__(VectorF)
        new.residue_frequencies._data = &(new._bg.f[0])
        new.residue_frequencies._owner = new
        new.residue_frequencies._n = new.residue_frequencies._shape[0] = new.alphabet.K
        return new


cdef class Builder:
    """A factory for constructing new HMMs from raw sequences.

    Attributes:
        alphabet (`~pyhmmer.easel.Alphabet`): The alphabet the builder is
            configured to use to convert sequences to HMMs.
        randomness (`~pyhmmer.easel.Randomness`): The random number generator
            being used by the builder.
        score_matrix (`str`): The name of the substitution matrix used to
            build HMMs for single sequences.
        popen (`float`): The *gap open* probability to use when building
            HMMs from single sequences.
        pextend (`float`): The *gap extend* probability to use when building
            HMMs from single sequences.

    .. versionadded:: 0.2.0

    .. versionchanged:: 0.4.2
       Added the `~Builder.randomness` attribute.

    """

    _ARCHITECTURE_STRATEGY = dict(BUILDER_ARCHITECTURE_STRATEGY)
    _WEIGHTING_STRATEGY = dict(BUILDER_WEIGHTING_STRATEGY)
    _EFFECTIVE_STRATEGY = dict(BUILDER_EFFECTIVE_STRATEGY)

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._bld = NULL
        self.alphabet = None

    def __init__(
        self,
        Alphabet alphabet,
        *,
        str architecture="fast",
        str weighting="pb",
        object effective_number="entropy",
        str prior_scheme="alphabet",
        float symfrac=0.5,
        float fragthresh=0.5,
        double wid=0.62,
        double esigma=45.0,
        double eid=0.62,
        int EmL=200,
        int EmN=200,
        int EvL=200,
        int EvN=200,
        int EfL=100,
        int EfN=200,
        double Eft=0.04,
        uint32_t seed=42,
        object ere=None,
        object popen=None,
        object pextend=None,
        object score_matrix=None,
        object window_length=None,
        object window_beta=None,
    ):
        """__init__(self, alphabet, *, architecture="fast", weighting="pb", effective_number="entropy", prior_scheme="alphabet", symfrac=0.5, fragthres=0.5, wid=0.62, esigma=45.0, eid=0.62, EmL=200, EmN=200, EvL=200, EvN=200, EfL=100, EfN=200, Eft=0.04, seed=42, ere=None, popen=None, pextend=None, score_matrix=None, window_length=None, window_beta=None)\n--\n

        Create a new sequence builder with the given configuration.

        Arguments:
            alphabet (`~pyhmmer.easel.Alphabet`): The alphabet the builder
                expects the sequences to be in.

        Keyword Arguments:
            architecture (`str`): The algorithm to use to determine the
                model architecture, either ``"fast"`` (the default), or
                ``"hand"``.
            weighting (`str`): The algorithm to use for relative sequence
                weighting, either ``"pb"`` (the default), ``"gsc"``,
                ``"blosum"``, ``"none"``, or ``"given"``.
            effective_number (`str`, `int`, or `float`): The algorithm to
                use to determine the effective sequence number, either
                ``"entropy"`` (the default), ``"exp"``, ``"clust"``, ``"none"``.
                A number can also be given directly to set the effective
                sequence number manually.
            prior_scheme (`str`): The choice of mixture Dirichlet prior when
                parameterizing  from counts, either ``"laplace"``
                or ``"alphabet"`` (the default).
            symfrac (`float`): The residue occurrence threshold for fast
                architecture determination.
            fragthresh (`float`): A threshold such that a sequence is called
                a fragment when :math:`L \\le fragthresh \times alen`.
            wid (`double`): The percent identity threshold for BLOSUM relative
                weighting.
            esigma (`float`): The minimum total relative entropy parameter
                for effective number entropy weights.
            eid (`float`): The percent identity threshold for effective
                number clustering.
            EmL (`int`): The length of sequences generated for MSV fitting.
            EmN (`int`): The number of sequences generated for MSV fitting.
            EvL (`int`): The lenght of sequences generated for Viterbi fitting.
            EvN (`int`): The number of sequences generated for Viterbi fitting.
            EfL (`int`): The lenght of sequences generated for Forward fitting.
            EfN (`int`): The number of sequences generated for Forward fitting.
            Eft (`float`): The tail mass used for Forward fitting.
            seed (`int`): The seed to use to initialize the internal random
                number generator. If ``0`` is given, an arbitrary seed will
                be chosen based on the current time.
            ere (`double`, optional): The relative entropy target for effective
                number weighting, or `None`.
            popen (`float`, optional): The *gap open* probability to use
                when building HMMs from single sequences. The default value
                depends on the alphabet: *0.02* for proteins,
                *0.03125* for nucleotides.
            pextend (`float`, optional): The *gap extend* probability to use
                when building HMMs from single sequences. Default depends on
                the alphabet: *0.4* for proteins, *0.75* for nucleotides.
            score_matrix (`str`, optional): The name of the score matrix to
                use when building HMMs from single sequences. The only
                allowed value for nucleotide alphabets is *DNA1*. For
                proteins, *PAM30*, *PAM70*, *PAM120*, *PAM240*, *BLOSUM45*,
                *BLOSUM50*, *BLOSUM62* (the default), *BLOSUM80* or
                *BLOSUM90* can be used.
            window_length (`float`, optional): The window length for
                nucleotide sequences, essentially the max expected hit length.
                *If given, takes precedence over* ``window_beta``.
            window_beta (`float`, optional): The tail mass at which window
                length is determined for nucleotide sequences.

        """
        # extract alphabet and create raw P7_BUILDER object
        self.alphabet = alphabet
        abcty = alphabet._abc.type
        with nogil:
            self._bld = libhmmer.p7_builder.p7_builder_Create(NULL, alphabet._abc)
        if self._bld == NULL:
            raise AllocationError("P7_BG", sizeof(P7_BG))

        # create a Randomness object exposing the internal RNG
        self.randomness = Randomness.__new__(Randomness)
        self.randomness._owner = self
        self.randomness._rng = self._bld.r

        # store the seed given at initialization and reseed the RNG
        self.seed = seed

        # set numeric parameters
        self._bld.symfrac = symfrac
        self._bld.fragthresh = fragthresh
        self._bld.wid = wid
        self._bld.esigma = esigma
        self._bld.eid = eid
        self._bld.EmL = EmL
        self._bld.EmN = EmN
        self._bld.EvL = EvL
        self._bld.EvN = EvN
        self._bld.EfL = EfL
        self._bld.EfN = EfN
        self._bld.Eft = Eft

        # set the architecture strategy
        self.architecture = architecture
        _arch = BUILDER_ARCHITECTURE_STRATEGY.get(architecture)
        if _arch is not None:
            self._bld.arch_strategy = _arch
        else:
            raise InvalidParameter(
                "architecture",
                architecture,
                choices=list(BUILDER_ARCHITECTURE_STRATEGY)
            )

        # set the weighting strategy
        self.weighting = weighting
        _weighting = BUILDER_WEIGHTING_STRATEGY.get(weighting)
        if _weighting is not None:
            self._bld.wgt_strategy = _weighting
        else:
            raise InvalidParameter(
                "weighting",
                weighting,
                choices=list(BUILDER_WEIGHTING_STRATEGY)
            )

        # set the effective sequence number strategy
        self.effective_number = effective_number
        if isinstance(effective_number, (int, float)):
            self._bld.effn_strategy = p7_effnchoice_e.p7_EFFN_SET
            self._bld.eset = effective_number
        elif isinstance(effective_number, str):
            effn = BUILDER_EFFECTIVE_STRATEGY.get(effective_number)
            if effn is None:
                raise InvalidParameter(
                    "effective_number",
                    effective_number,
                    choices=list(BUILDER_EFFECTIVE_STRATEGY) + [int, float]
                )
            self._bld.effn_strategy = BUILDER_EFFECTIVE_STRATEGY[effective_number]
        else:
            ty = type(effective_number).__name__
            raise TypeError(f"Expected str, int or float, found {ty}")

        # set the RE target if given one, otherwise the default
        # is alphabet dependent
        if ere is not None:
            self._bld.re_target = ere
        elif alphabet.is_nucleotide():
            self._bld.re_target = libhmmer.p7_ETARGET_DNA
        elif alphabet.is_amino():
            self._bld.re_target = libhmmer.p7_ETARGET_AMINO
        else:
            self._bld.re_target = libhmmer.p7_ETARGET_OTHER

        # set the prior scheme (this involves re-allocating a new prior)
        self.prior_scheme = prior_scheme
        libhmmer.p7_prior.p7_prior_Destroy(self._bld.prior)
        self._bld.prior = NULL
        if prior_scheme is not None:
            if prior_scheme == "laplace":
                self._bld.prior = libhmmer.p7_prior.p7_prior_CreateLaplace(self.alphabet._abc)
            elif prior_scheme == "alphabet":
                if alphabet.is_amino():
                    self._bld.prior = libhmmer.p7_prior.p7_prior_CreateAmino()
                elif alphabet.is_nucleotide():
                    self._bld.prior = libhmmer.p7_prior.p7_prior_CreateNucleic()
                else:
                    self._bld.prior = libhmmer.p7_prior.p7_prior_CreateLaplace(self.alphabet._abc)
            else:
                raise InvalidParameter("prior_scheme", prior_scheme, choices=["laplace", "alphabet", None])
            if self._bld.prior == NULL:
                raise AllocationError("P7_PRIOR", sizeof(P7_PRIOR))

        # set the gap probabilities and score matrix using alphabet-specific
        # defaults (this is only used when building a HMM for a single sequence)
        if alphabet.is_nucleotide():
            self.score_matrix = "DNA1" if score_matrix is None else score_matrix
            self.popen = 0.03125 if popen is None else popen
            self.pextend = 0.75 if pextend is None else pextend
        else:
            self.score_matrix = "BLOSUM62" if score_matrix is None else score_matrix
            self.popen = 0.02 if popen is None else popen
            self.pextend = 0.4 if pextend is None else pextend

        # set the window options
        self.window_length = window_length
        if window_beta is None:
            self.window_beta = libhmmer.p7_builder.p7_DEFAULT_WINDOW_BETA
        else:
            self.window_beta = window_beta

    def __dealloc__(self):
        libhmmer.p7_builder.p7_builder_Destroy(self._bld)

    def __copy__(self):
        return self.copy()

    # --- Properties ---------------------------------------------------------

    @property
    def seed(self):
        """`int`: The seed given at builder initialization.

        Setting this attribute to a different value will cause the internal
        random number generator to be reseeded immediately.

        .. versionchanged:: 0.4.2
           Avoid shadowing initial null seeds given on builder initialization.

        """
        return self._seed

    @seed.setter
    def seed(self, uint32_t seed):
        self._seed = seed
        self._bld.do_reseeding = seed != 0
        self.randomness.seed(seed)

    @property
    def window_length(self):
        """`int` or `None`: The window length for nucleotide sequences.
        """
        assert self._bld != NULL
        return None if self._bld.w_len == -1 else self._bld.w_len

    @window_length.setter
    def window_length(self, object window_length):
        assert self._bld != NULL
        if window_length is None:
            self._bld.w_len = -1
        elif window_length > 3:
            self._bld.w_len = window_length
        else:
            raise InvalidParameter("window_length", window_length, hint="integer greater than 3 or None")

    @property
    def window_beta(self):
        """`float`: The tail mass at which window length is determined.
        """
        assert self._bld != NULL
        return self._bld.w_beta

    @window_beta.setter
    def window_beta(self, double window_beta):
        assert self._bld != NULL
        if window_beta > 1 or window_beta < 0:
            raise InvalidParameter("window_beta", window_beta, hint="real number between 0 and 1")
        self._bld.w_beta = window_beta

    # --- Methods ------------------------------------------------------------

    cpdef tuple build(
        self,
        DigitalSequence sequence,
        Background background,
    ):
        """Build a new HMM from ``sequence`` using the builder configuration.

        Arguments:
            sequence (`~pyhmmer.easel.DigitalSequence`): A single biological
                sequence in digital mode to build a HMM with.
            background (`~pyhmmer.plan7.Background`): The background model
                to use to create the HMM.

        Returns:
            (`HMM`, `Profile`, `OptimizedProfile`): A tuple containing the
            new HMM as well as profiles to be used directly in a `Pipeline`.

        Raises:
            `~pyhmmer.errors.AlphabetMismatch`: When either ``sequence`` or
                ``background`` have the wrong alphabet for this builder.
            `ValueError`: When the HMM cannot be created successfully
                from the input; the error message contains more details.

        Hint:
            The score matrix and the gap probabilities used here can be set
            when initializing the builder, or changed by setting a new value
            to the right attribute::

                >>> alphabet = easel.Alphabet.amino()
                >>> background = plan7.Background(alphabet)
                >>> builder = plan7.Builder(alphabet, popen=0.05)
                >>> builder.score_matrix = "BLOSUM90"
                >>> hmm, _, _ = builder.build(proteins[0], background)

        .. versionchanged:: 0.4.6
           Sets the `HMM.creation_time` attribute with the current time.

        .. versionchanged:: 0.4.10
           The score system is now cached between calls to `Builder.build`.

        """
        assert self._bld != NULL

        cdef int              status
        cdef HMM              hmm     = HMM.__new__(HMM)
        cdef Profile          profile = Profile.__new__(Profile)
        cdef OptimizedProfile opti    = OptimizedProfile.__new__(OptimizedProfile)
        cdef bytes            mx      = self.score_matrix.encode('utf-8')
        cdef char*            mx_ptr  = <char*> mx
        cdef str              msg

        # reseed RNG used by the builder if needed
        if self._bld.do_reseeding:
            self.randomness.seed(self.seed)

        # check alphabet and assign it to newly created objects
        hmm.alphabet = profile.alphabet = opti.alphabet = self.alphabet
        if not self.alphabet._eq(background.alphabet):
            raise AlphabetMismatch(self.alphabet, background.alphabet)
        if not self.alphabet._eq(sequence.alphabet):
            raise AlphabetMismatch(self.alphabet, sequence.alphabet)

        # only load the score system if it hasn't been loaded already,
        # or if a different score system is in use
        if self._bld.S == NULL or strcmp(self._bld.S.name, mx_ptr) != 0:
            with nogil:
                status = libhmmer.p7_builder.p7_builder_LoadScoreSystem(
                    self._bld,
                    mx_ptr,
                    self.popen,
                    self.pextend,
                    background._bg
                )
            if status == libeasel.eslEINVAL:
                msg = self._bld.errbuf.decode('utf-8', 'ignore')
                raise RuntimeError(f"failed to set single sequence score system: {msg}")
            elif status != libeasel.eslOK:
                raise UnexpectedError(status, "p7_builder_LoadScoreSystem")
        else:
            self._bld.popen = self.popen
            self._bld.pextend = self.pextend

        # build HMM and profiles
        with nogil:
            status = libhmmer.p7_builder.p7_SingleBuilder(
                self._bld,
                sequence._sq,
                background._bg,
                &hmm._hmm, # HMM
                NULL, # traceback
                &profile._gm, # profile,
                &opti._om, # optimized profile
            )
        if status == libeasel.eslOK:
            hmm.command_line = " ".join(sys.argv)
        elif status == libeasel.eslEINVAL:
            msg = self._bld.errbuf.decode("utf-8", "replace")
            raise ValueError(f"Could not build HMM: {msg}")
        else:
            raise UnexpectedError(status, "p7_SingleBuilder")

        # update max length (used in `nhmmer`) which is only set on `hmm`
        profile._gm.max_length = opti._om.max_length = hmm._hmm.max_length

        # return newly built HMM, profile and optimized profile
        return hmm, profile, opti

    cpdef tuple build_msa(
        self,
        DigitalMSA msa,
        Background background,
    ):
        """Build a new HMM from ``msa`` using the builder configuration.

        Arguments:
            msa (`~pyhmmer.easel.DigitalMSA`): A multiple sequence
                alignment in digital mode to build a HMM with.
            background (`~pyhmmer.plan7.Background`): The background model
                to use to create the HMM.

        Returns:
            (`HMM`, `Profile`, `OptimizedProfile`): A tuple containing the
            new HMM as well as profiles to be used directly in a `Pipeline`.

        Raises:
            `~pyhmmer.errors.AlphabetMismatch`: When either ``msa`` or
                ``background`` have the wrong alphabet for this builder.
            `ValueError`: When the HMM cannot be created successfully
                from the input; the error message contains more details.

        Caution:
            HMMER requires that every HMM has a name, so the `Builder` will
            attempt to use the name of the sequence to name the HMM. Passing
            an MSA without a name will result in an error::

                >>> alphabet = easel.Alphabet.amino()
                >>> msa = easel.TextMSA(sequences=[
                ...   easel.TextSequence(name=b"ustiA", sequence="YAIG"),
                ...   easel.TextSequence(name=b"ustiB", sequence="YVIG")
                ... ])
                >>> builder = plan7.Builder(alphabet)
                >>> background = plan7.Background(alphabet)
                >>> builder.build_msa(msa.digitize(alphabet), background)
                Traceback (most recent call last):
                  ...
                ValueError: Could not build HMM: Unable to name the HMM.

        .. versionadded:: 0.3.0

        .. versionchanged:: 0.4.6
           Sets the `HMM.creation_time` attribute with the current time.

        """
        assert self._bld != NULL

        cdef int              status
        cdef HMM              hmm     = HMM.__new__(HMM)
        cdef Profile          profile = Profile.__new__(Profile)
        cdef OptimizedProfile opti    = OptimizedProfile.__new__(OptimizedProfile)
        cdef str              msg
        cdef size_t           k

        # unset builder probabilities and score system if they were set
        # by a call to `build`, since we won't be needing them here
        self._bld.popen = self._bld.pextend = -1
        if self._bld.S != NULL:
            libeasel.scorematrix.esl_scorematrix_Destroy(self._bld.S)
            self._bld.S = NULL
        if self._bld.Q != NULL:
            libeasel.dmatrix.esl_dmatrix_Destroy(self._bld.Q)
            self._bld.Q = NULL

        # reseed RNG used by the builder if needed
        if self._bld.do_reseeding:
            self.randomness.seed(self.seed)

        # check alphabet and assign it to newly created objects
        hmm.alphabet = profile.alphabet = opti.alphabet = self.alphabet
        if not self.alphabet._eq(background.alphabet):
            raise AlphabetMismatch(self.alphabet, background.alphabet)
        if not self.alphabet._eq(msa.alphabet):
            raise AlphabetMismatch(self.alphabet, msa.alphabet)

        # build HMM
        with nogil:
            # build HMM and profiles
            status = libhmmer.p7_builder.p7_Builder(
                self._bld,
                msa._msa,
                background._bg,
                &hmm._hmm, # HMM
                NULL, # traceback
                &profile._gm, # profile,
                &opti._om, # optimized profile
                NULL # modified msa
            )
        if status == libeasel.eslOK:
            hmm.command_line = " ".join(sys.argv)
        elif status == libeasel.eslEINVAL:
            msg = self._bld.errbuf.decode("utf-8", "replace")
            raise ValueError(f"Could not build HMM: {msg}")
        else:
            raise UnexpectedError(status, "p7_Builder")

        # update max length (used in `nhmmer`) which is only set on `hmm`
        profile._gm.max_length = opti._om.max_length = hmm._hmm.max_length

        # return newly built HMM, profile and optimized profile
        return hmm, profile, opti

    cpdef Builder copy(self):
        """Create a duplicate `Builder` instance with the same arguments.
        """
        assert self._bld != NULL
        return Builder(
            self.alphabet,
            architecture=self.architecture,
            weighting=self.weighting,
            effective_number=self.effective_number,
            prior_scheme=self.prior_scheme,
            symfrac=self._bld.symfrac,
            fragthresh=self._bld.fragthresh,
            wid=self._bld.wid,
            esigma=self._bld.esigma,
            eid=self._bld.eid,
            EmL=self._bld.EmL,
            EmN=self._bld.EmN,
            EvL=self._bld.EvL,
            EvN=self._bld.EvN,
            EfL=self._bld.EfL,
            EfN=self._bld.EfN,
            Eft=self._bld.Eft,
            seed=self.seed,
            ere=self._bld.re_target,
            popen=self.popen,
            pextend=self.pextend,
            window_length=self.window_length,
            window_beta=self.window_beta,
        )


@cython.freelist(8)
@cython.no_gc_clear
cdef class Cutoffs:
    """A mutable view over the score cutoffs of a `HMM` or a `Profile`.

    .. versionadded:: 0.4.6

    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._owner = None
        self._flags = NULL
        self._cutoffs = NULL
        self._is_profile = True

    def __init__(self, object owner):
        """__init__(self, owner)\n--\n
        """
        cdef str ty
        if isinstance(owner, Profile):
            self._cutoffs = &(<Profile> owner)._gm.cutoff
            self._flags = NULL
            self._owner = owner
            self._is_profile = True
        elif isinstance(owner, HMM):
            self._cutoffs = &(<HMM> owner)._hmm.cutoff
            self._flags = &(<HMM> owner)._hmm.flags
            self._owner = owner
            self._is_profile = False
        else:
            ty = type(owner).__name__
            raise TypeError(f"expected Profile or HMM, found {ty}")

    def __copy__(self):
        cdef Cutoffs c
        c = Cutoffs.__new__(Cutoffs)
        c._owner = self._owner
        c._cutoffs = self._cutoffs
        c._flags = self._flags
        c._owner = self._owner
        return c

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"<{ty} of {self._owner!r}>"

    def __eq__(self, object other):
        if isinstance(other, Cutoffs):
            return self.as_vector() == other.as_vector()
        return NotImplemented

    # --- Properties ---------------------------------------------------------

    @property
    def gathering(self):
        """`tuple` of `float`, or `None`: The gathering thresholds, if any.

        Example:
            This property can be used to set the gathering cutoffs by
            passing it an iterable of two `float`::

                >>> thioesterase.cutoffs.gathering = (180.0, 120.0)
                >>> thioesterase.cutoffs.gathering_available()
                True
                >>> thioesterase.cutoffs.gathering
                (180.0, 120.0)

            Set the attribute to `None` or delete it with ``del`` to clear
            the gathering thresholds::

                >>> thioesterase.cutoffs.gathering = None
                >>> thioesterase.cutoffs.gathering_available()
                False

        .. versionadded:: 0.4.8

        """
        if self.gathering_available():
            return (
                self._cutoffs[0][<int> p7_cutoffs_e.p7_GA1],
                self._cutoffs[0][<int> p7_cutoffs_e.p7_GA2]
            )
        return None

    @gathering.setter
    def gathering(self, object gathering):
        assert self._cutoffs != NULL

        if gathering is None:
            del self.gathering
            return

        gathering1, gathering2 = gathering
        self._cutoffs[0][<int> p7_cutoffs_e.p7_GA1] = gathering1
        self._cutoffs[0][<int> p7_cutoffs_e.p7_GA2] = gathering2

        if not self._is_profile:
            self._flags[0] = self._flags[0] | p7H_GA

    @gathering.deleter
    def gathering(self):
        assert self._cutoffs != NULL

        self._cutoffs[0][<int> p7_cutoffs_e.p7_GA1] = p7_CUTOFF_UNSET
        self._cutoffs[0][<int> p7_cutoffs_e.p7_GA2] = p7_CUTOFF_UNSET

        if not self._is_profile:
            self._flags[0] = self._flags[0] & (~p7H_GA)

    @property
    def gathering1(self):
        """`float` or `None`: The first gathering threshold, if any.
        """
        if self.gathering_available():
            return self._cutoffs[0][<int> p7_cutoffs_e.p7_GA1]
        return None

    @property
    def gathering2(self):
        """`float` or `None`: The second gathering threshold, if any.
        """
        if self.gathering_available():
            return self._cutoffs[0][<int> p7_cutoffs_e.p7_GA2]
        return None

    @property
    def trusted(self):
        """`tuple` of `float`, or `None`: The trusted cutoffs, if available.

        .. versionadded:: 0.4.8

        """
        if self.trusted_available():
            return (
                self._cutoffs[0][<int> p7_cutoffs_e.p7_TC1],
                self._cutoffs[0][<int> p7_cutoffs_e.p7_TC2]
            )
        return None

    @trusted.setter
    def trusted(self, object trusted):
        assert self._cutoffs != NULL

        if trusted is None:
            del self.trusted
            return

        trusted1, trusted2 = trusted
        self._cutoffs[0][<int> p7_cutoffs_e.p7_TC1] = trusted1
        self._cutoffs[0][<int> p7_cutoffs_e.p7_TC2] = trusted2

        if not self._is_profile:
            self._flags[0] = self._flags[0] | p7H_TC

    @trusted.deleter
    def trusted(self):
        assert self._cutoffs != NULL

        self._cutoffs[0][<int> p7_cutoffs_e.p7_TC1] = p7_CUTOFF_UNSET
        self._cutoffs[0][<int> p7_cutoffs_e.p7_TC2] = p7_CUTOFF_UNSET

        if not self._is_profile:
            self._flags[0] = self._flags[0] & (~p7H_TC)

    @property
    def trusted1(self):
        """`float` or `None`: The first trusted score cutoff, if any.
        """
        if self.trusted_available():
            return self._cutoffs[0][<int> p7_cutoffs_e.p7_TC1]
        return None

    @property
    def trusted2(self):
        """`float` or `None`: The second trusted score cutoff, if any.
        """
        if self.trusted_available():
            return self._cutoffs[0][<int> p7_cutoffs_e.p7_TC2]
        return None

    @property
    def noise(self):
        """`tuple` of `float`, or `None`: The noise cutoffs, if available.

        .. versionadded:: 0.4.8

        """
        if self.noise_available():
            return (
                self._cutoffs[0][<int> p7_cutoffs_e.p7_NC1],
                self._cutoffs[0][<int> p7_cutoffs_e.p7_NC2]
            )
        return None

    @noise.setter
    def noise(self, object noise):
        assert self._cutoffs != NULL

        if noise is None:
            del self.noise
            return

        noise1, noise2 = noise
        self._cutoffs[0][<int> p7_cutoffs_e.p7_NC1] = noise1
        self._cutoffs[0][<int> p7_cutoffs_e.p7_NC2] = noise2

        if not self._is_profile:
            self._flags[0] = self._flags[0] | p7H_NC

    @noise.deleter
    def noise(self):
        assert self._cutoffs != NULL

        self._cutoffs[0][<int> p7_cutoffs_e.p7_NC1] = p7_CUTOFF_UNSET
        self._cutoffs[0][<int> p7_cutoffs_e.p7_NC2] = p7_CUTOFF_UNSET

        if not self._is_profile:
            self._flags[0] = self._flags[0] & (~p7H_NC)

    @property
    def noise1(self):
        """`float` or `None`: The first noise cutoff, if any.
        """
        if self.noise_available():
            return self._cutoffs[0][<int> p7_cutoffs_e.p7_NC1]
        return None

    @property
    def noise2(self):
        """`float` or `None`: The second noise cutoff, if any.
        """
        if self.noise_available():
            return self._cutoffs[0][<int> p7_cutoffs_e.p7_NC2]
        return None

    # --- Methods ------------------------------------------------------------

    cpdef bint gathering_available(self):
        """Check whether the gathering thresholds are available.
        """
        assert self._cutoffs != NULL
        if self._is_profile:
            return (
                    self._cutoffs[0][<int> p7_cutoffs_e.p7_GA1] != p7_CUTOFF_UNSET
                and self._cutoffs[0][<int> p7_cutoffs_e.p7_GA2] != p7_CUTOFF_UNSET
            )
        else:
            return (self._flags[0] & p7H_GA) != 0

    cpdef bint trusted_available(self):
        """Check whether the trusted cutoffs are available.
        """
        assert self._cutoffs != NULL
        if self._is_profile:
            return (
                    self._cutoffs[0][<int> p7_cutoffs_e.p7_TC1] != p7_CUTOFF_UNSET
                and self._cutoffs[0][<int> p7_cutoffs_e.p7_TC2] != p7_CUTOFF_UNSET
            )
        else:
            return (self._flags[0] & p7H_TC) != 0

    cpdef bint noise_available(self):
        """Check whether the noise cutoffs are available.
        """
        assert self._cutoffs != NULL
        if self._is_profile:
            return (
                    self._cutoffs[0][<int> p7_cutoffs_e.p7_NC1] != p7_CUTOFF_UNSET
                and self._cutoffs[0][<int> p7_cutoffs_e.p7_NC2] != p7_CUTOFF_UNSET
            )
        else:
            return (self._flags[0] & p7H_NC) != 0

    cpdef VectorF as_vector(self):
        """Return a view over the score cutoffs as a `~pyhmmer.easel.VectorF`.
        """
        assert self._cutoffs != NULL

        cdef VectorF new

        new = VectorF.__new__(VectorF)
        new._owner = self
        new._n = new._shape[0] = p7_NCUTOFFS
        new._data = <float*> self._cutoffs
        return new


cdef class Domain:
    """A single domain in a query `~pyhmmer.plan7.Hit`.

    Attributes:
        hit (`~pyhmmer.plan7.Hit`): The hit this domains is part of.
        alignment (`~pyhmmer.plan7.Alignment`): The alignment of this domain
            to a target sequence.

    .. versionadded:: 0.6.1
       `pickle` protocol support.

    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self, Hit hit, size_t index):
        self._dom = &hit._hit.dcl[index]
        self.hit = hit
        self.alignment = Alignment(self)

    def __getstate__(self):
        assert self._dom != NULL

        cdef int       status
        cdef uint32_t  n      = 0
        cdef uint32_t  nalloc = 0
        cdef uint8_t*  buffer = NULL
        cdef VectorU8  vec

        with nogil:
            status = libhmmer.p7_domain.p7_domain_Serialize(self._dom, &buffer, &n, &nalloc)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_domain_Serialize")

        vec = VectorU8.__new__(VectorU8)
        vec._owner = None
        vec._n = vec._shape[0] = nalloc
        vec._data = <void*> buffer
        return vec

    def __setstate__(self, uint8_t[::1] state):
        cdef int      status
        cdef uint32_t offset = 0

        if self._dom == NULL:
            self._dom = libhmmer.p7_domain.p7_domain_Create_empty()
            if self._dom == NULL:
                raise AllocationError("P7_DOMAIN", sizeof(P7_DOMAIN))

        with nogil:
            status = libhmmer.p7_domain.p7_domain_Deserialize(&state[0], &offset, self._dom)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_domain_Deserialize")

    # --- Properties ---------------------------------------------------------

    @property
    def env_from(self):
        """`int`: The start coordinate of the domain envelope.
        """
        assert self._dom != NULL
        return self._dom.ienv

    @property
    def env_to(self):
        """`int`: The end coordinate of the domain envelope.
        """
        assert self._dom != NULL
        return self._dom.jenv

    @property
    def strand(self):
        """`str` or `None`: The strand where the domain is located.

        When running a search with the `LongTargetsPipeline`, both strands
        of each target sequence are processed (unless disabled), so the
        domain may be located on either strand, either ``+`` or ``-``.
        For default `Pipeline` searches, this is always `None`.

        .. versionadded:: 0.10.8

        """
        assert self._dom != NULL
        if not self.hit.hits._pli.long_targets:
            return None
        return "+" if self._dom.iali < self._dom.jali else "-"

    @property
    def score(self):
        """`float`: The overall score in bits, *null2*-corrected.
        """
        assert self._dom != NULL
        return self._dom.bitscore

    @property
    def bias(self):
        """`float`: The *null2* score contribution to the domain score.
        """
        assert self._dom != NULL
        return self._dom.dombias * libeasel.eslCONST_LOG2R

    @property
    def correction(self):
        """`float`: The *null2* score when calculating a per-domain score.
        """
        assert self._dom != NULL
        return self._dom.domcorrection * libeasel.eslCONST_LOG2R

    @property
    def envelope_score(self):
        """`float`: The forward score in the envelope, without *null2* correction.
        """
        assert self._dom != NULL
        return self._dom.envsc * libeasel.eslCONST_LOG2R

    @property
    def c_evalue(self):
        """`float`: The conditional e-value for the domain.
        """
        assert self._dom != NULL
        if self.hit.hits.long_targets:
            return exp(self._dom.lnP)
        else:
            return exp(self._dom.lnP) * self.hit.hits._pli.domZ

    @property
    def i_evalue(self):
        """`float`: The independent e-value for the domain.
        """
        assert self._dom != NULL
        if self.hit.hits.long_targets:
            return exp(self._dom.lnP)
        else:
            return exp(self._dom.lnP) * self.hit.hits._pli.Z

    @property
    def pvalue(self):
        """`float`: The p-value of the domain bitscore.
        """
        assert self._dom != NULL
        return exp(self._dom.lnP)

    @property
    def included(self):
        """`bool`: Whether this domain is marked as *included*.

        .. versionadded:: 0.7.0

        """
        assert self._dom != NULL
        return self._dom.is_included

    @included.setter
    def included(self, bint included):
        assert self._dom != NULL
        assert self.hit._hit != NULL
        if included:
            if not self._dom.is_included:
                self.hit._hit.nincluded += 1
            if not self._dom.is_reported:
                self.hit._hit.nreported += 1
            self._dom.is_included = self._dom.is_reported = True
        else:
            if self._dom.is_included:
                self.hit._hit.nincluded -= 1
            self._dom.is_included = False

    @property
    def reported(self):
        """`bool`: Whether this domain is marked as *reported*.

        .. versionadded:: 0.7.0

        """
        assert self._dom != NULL
        return self._dom.is_reported

    @reported.setter
    def reported(self, bint reported):
        assert self._dom != NULL
        assert self.hit._hit != NULL
        if reported:
            if not self._dom.is_reported:
                self.hit._hit.nreported += 1
            self._dom.is_reported = True
        else:
            if self._dom.is_reported:
                self.hit._hit.nreported -= 1
            if self._dom.is_included:
                self.hit._hit.nincluded -= 1
            self._dom.is_reported = self._dom.is_included = False


cdef class Domains:
    """A read-only view over the domains of a single `~pyhmmer.plan7.Hit`.

    Attributes:
        hit (`~pyhmmer.plan7.Hit`): The target hit these domains belong hit.

    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self, Hit hit):
        self.hit = hit

    def __len__(self):
        assert self.hit._hit != NULL
        return self.hit._hit.ndom

    def __getitem__(self, int index):
        assert self.hit._hit != NULL
        if index < 0:
            index += self.hit._hit.ndom
        if index >= self.hit._hit.ndom or index < 0:
            raise IndexError("list index out of range")
        return Domain(self.hit, <size_t> index)

    # --- Properties ---------------------------------------------------------

    @property
    def included(self):
        """iterator of `Domain`: An iterator over *included* domains only.

        .. versionadded:: 0.7.0

        """
        assert self.hit._hit != NULL
        return SizedIterator(
            self.hit._hit.nincluded,
            filter(operator.attrgetter("included"), self)
        )

    @property
    def reported(self):
        """iterator of `Domain`: An iterator over *reported* domains only.

        .. versionadded:: 0.7.0

        """
        assert self.hit._hit != NULL
        return SizedIterator(
            self.hit._hit.nreported,
            filter(operator.attrgetter("reported"), self)
        )

@cython.freelist(8)
@cython.no_gc_clear
cdef class EvalueParameters:
    """A mutable view over the e-value parameters of a `HMM` or a `Profile`.

    The e-value for each filter is estimated based off a maximum likelihood
    distribution fitted for each profile HMM, either a Gumbel distribution
    for the MSV and Viterbi filters, or an exponential distribution for the
    Forward filter.

    .. versionadded:: 0.4.6

    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._owner = None
        self._evparams = NULL

    def __init__(self, object owner):
        """__init__(self, owner)\n--\n
        """
        cdef str ty
        if isinstance(owner, Profile):
            self._evparams = &(<Profile> owner)._gm.evparam
            self._owner = owner
        elif isinstance(owner, HMM):
            self._evparams = &(<HMM> owner)._hmm.evparam
            self._owner = owner
        else:
            ty = type(owner).__name__
            raise TypeError(f"expected Profile or HMM, found {ty}")

    def __copy__(self):
        cdef EvalueParameters ev = EvalueParameters.__new__(EvalueParameters)
        ev._owner = self._owner
        ev._evparams = self._evparams
        return ev

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"<{ty} of {self._owner!r}>"

    def __eq__(self, object other):
        if isinstance(other, EvalueParameters):
            return self.as_vector() == other.as_vector()
        return NotImplemented

    # --- Properties ---------------------------------------------------------

    @property
    def m_mu(self):
        r"""`float` or `None`: The :math:`\mu` parameter for the MSV filter distribution.
        """
        assert self._evparams != NULL
        cdef float m = self._evparam[0][<int> p7_evparams_e.p7_MMU]
        return None if m == p7_EVPARAM_UNSET else m

    @m_mu.setter
    def m_mu(self, object m):
        assert self._evparams != NULL
        if m is None:
            self._evparams[0][<int> p7_evparams_e.p7_MMU] = p7_EVPARAM_UNSET
        else:
            self._evparams[0][<int> p7_evparams_e.p7_MMU] = m

    @property
    def m_lambda(self):
        r"""`float` or `None`: The :math:`\lambda` parameter for the MSV filter distribution.
        """
        assert self._evparams != NULL
        cdef float l = self._evparams[0][<int> p7_evparams_e.p7_MLAMBDA]
        return None if l == p7_EVPARAM_UNSET else l

    @m_lambda.setter
    def m_lambda(self, object l):
        assert self._evparams != NULL
        if l is None:
            self._evparams[0][<int> p7_evparams_e.p7_MLAMBDA] = p7_EVPARAM_UNSET
        else:
            self._evparams[0][<int> p7_evparams_e.p7_MLAMBDA] = l

    @property
    def v_mu(self):
        r"""`float` or `None`: The :math:`\mu` parameter for the Viterbi filter distribution.
        """
        assert self._evparams != NULL
        cdef float m = self._evparams[0][<int> p7_evparams_e.p7_VMU]
        return None if m == p7_EVPARAM_UNSET else m

    @v_mu.setter
    def v_mu(self, object m):
        assert self._evparams != NULL
        if m is None:
            self._evparams[0][<int> p7_evparams_e.p7_VMU] = p7_EVPARAM_UNSET
        else:
            self._evparams[0][<int> p7_evparams_e.p7_VMU] = m

    @property
    def v_lambda(self):
        r"""`float` or `None`: The :math:`\lambda` parameter for the Viterbi filter distribution.
        """
        assert self._evparams != NULL
        cdef float l = self._evparams[0][<int> p7_evparams_e.p7_VLAMBDA]
        return None if l == p7_EVPARAM_UNSET else l

    @v_lambda.setter
    def v_lambda(self, object l):
        assert self._evparams != NULL
        if l is None:
            self._evparams[0][<int> p7_evparams_e.p7_VLAMBDA] = p7_EVPARAM_UNSET
        else:
            self._evparams[0][<int> p7_evparams_e.p7_VLAMBDA] = l

    @property
    def f_tau(self):
        r"""`float` or `None`: The :math:`\tau` parameter for the Forward filter distribution.
        """
        assert self._evparams != NULL
        cdef float t = self._evparams[0][<int> p7_evparams_e.p7_FTAU]
        return None if t == p7_EVPARAM_UNSET else t

    @f_tau.setter
    def f_tau(self, object t):
        assert self._evparams != NULL
        if t is None:
            self._evparams[0][<int> p7_evparams_e.p7_FTAU] = p7_EVPARAM_UNSET
        else:
            self._evparams[0][<int> p7_evparams_e.p7_FTAU] = t

    @property
    def f_lambda(self):
        r"""`float` or `None`: The :math:`\lambda` parameter for the Forward filter distribution.
        """
        assert self._evparams != NULL
        cdef float l = self._evparams[0][<int> p7_evparams_e.p7_FLAMBDA]
        return None if l == p7_EVPARAM_UNSET else l

    @f_lambda.setter
    def f_lambda(self, object l):
        assert self._evparams != NULL
        if l is None:
            self._evparams[0][<int> p7_evparams_e.p7_FLAMBDA] = p7_EVPARAM_UNSET
        else:
            self._evparams[0][<int> p7_evparams_e.p7_FLAMBDA] = l

    # --- Methods ------------------------------------------------------------

    cpdef VectorF as_vector(self):
        """Return a view over the e-value parameters as a `~pyhmmer.easel.VectorF`.
        """
        assert self._evparams != NULL

        cdef VectorF new

        new = VectorF.__new__(VectorF)
        new._owner = self
        new._n = new._shape[0] = p7_NEVPARAM
        new._data = <float*> self._evparams
        return new


cdef class Hit:
    r"""A high-scoring database hit found by the comparison pipeline.

    A hit is obtained in HMMER for every target where one or more
    significant domain alignment was found by a `Pipeline`. A `Hit` comes
    with a *score*, which is obtained after correcting of the individual
    bit scores of all its domains; a *P-value*, which is computed by
    testing the likelihood to obtain the same alignment using a random
    background model; and an *E-value*, which is obtained after Bonferonni
    correction of the *p-value*, taking into account the total number of
    targets in the target database.

    Hits also store several information as flags. `Hit.included` and
    `Hit.reported` show whether a `Hit` is considered for inclusion
    (resp. reporting) with respects to the thresholds defined on the
    original `~pyhmmer.plan7.Pipeline`. These flags can be modified
    manually to force inclusion or exclusion of certains hits independently
    of their score or E-value. The `~pyhmmer.plan7.TopHits.write` method
    of `~pyhmmer.plan7.TopHits` objects will only write a line for hits
    marked as *reported*. Included hits are necessarily reported:

    .. math::

        \text{included} \implies \text{reported}

    When used during an iterative search, hits can also be marked as dropped
    by setting the `Hit.dropped` flag to `False`. Dropped hits will not be
    used for building HMMs during the next iteration. Hits newly found in an
    iteration will be marked as *new* with the `Hit.new` flag. `Hit.dropped`
    and `Hit.included` are mutually exclusive, and setting one will unset
    the other. Dropped hits can be reported, but are not included:

    .. math::

        \text{dropped} \implies \neg \text{included}

    When running a long target pipeline, some hits may appear as
    duplicates if they were found across multiple windows. These hits
    will be marked as duplicates with the `Hit.duplicate` flag. Duplicate
    hits are neither reported nor included:

    .. math::

        \text{duplicate} \implies \neg \text{reported}

    .. versionadded:: 0.6.1
       `pickle` protocol support.

    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self, TopHits hits, size_t index):
        assert hits._th != NULL
        assert index < hits._th.N
        self.hits = hits
        self._hit = hits._th.hit[index]

    def __getstate__(self):
        assert self._hit != NULL

        cdef int       status
        cdef uint32_t  n      = 0
        cdef uint32_t  nalloc = 0
        cdef uint8_t*  buffer = NULL
        cdef VectorU8  vec

        with nogil:
            status = libhmmer.p7_hit.p7_hit_Serialize(self._hit, &buffer, &n, &nalloc)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_hit_Serialize")

        vec = VectorU8.__new__(VectorU8)
        vec._owner = None
        vec._n = vec._shape[0] = nalloc
        vec._data = <void*> buffer
        return vec

    def __setstate__(self, uint8_t[::1] state):
        cdef int      status
        cdef uint32_t offset = 0

        if self._hit == NULL:
            self._hit = libhmmer.p7_hit.p7_hit_Create_empty()
            if self._hit == NULL:
                raise AllocationError("P7_HIT", sizeof(P7_HIT))

        with nogil:
            status = libhmmer.p7_hit.p7_hit_Deserialize(&state[0], &offset, self._hit)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_hit_Deserialize")


    # --- Properties ---------------------------------------------------------

    @property
    def name(self):
        """`bytes`: The name of the database hit.
        """
        assert self._hit != NULL
        assert self._hit.name != NULL
        return <bytes> self._hit.name

    @name.setter
    def name(self, bytes name not None):
        assert self._hit != NULL
        free(self._hit.name)
        self._hit.name = strdup(<const char*> name)
        if self._hit.name == NULL:
            raise AllocationError("char", sizeof(char), strlen(name))

    @property
    def accession(self):
        """`bytes` or `None`: The accession of the database hit, if any.
        """
        assert self._hit != NULL
        if self._hit.acc == NULL:
            return None
        return <bytes> self._hit.acc

    @accession.setter
    def accession(self, bytes accession):
        assert self._hit != NULL
        free(self._hit.acc)
        if accession is None:
            self._hit.acc = NULL
        else:
            self._hit.acc = strdup(<const char*> accession)
            if self._hit.acc == NULL:
                raise AllocationError("char", sizeof(char), strlen(accession))

    @property
    def description(self):
        """`bytes` or `None`: The description of the database hit, if any.
        """
        assert self._hit != NULL
        if self._hit.desc == NULL:
            return None
        return <bytes> self._hit.desc

    @description.setter
    def description(self, bytes description):
        assert self._hit != NULL
        free(self._hit.desc)
        if description is None:
            self._hit.desc = NULL
        else:
            self._hit.desc = strdup(<const char*> description)
            if self._hit.desc == NULL:
                raise AllocationError("char", sizeof(char), strlen(description))

    @property
    def length(self):
        """`int`: The length of the database hit.

        .. versionadded:: 0.10.5

        """
        cdef P7_DOMAIN* domain = &self._hit.dcl[self._hit.best_domain]
        if self.hits._pli.mode == p7_pipemodes_e.p7_SEARCH_SEQS:
            return domain.ad.L
        else:
            return domain.ad.M

    @property
    def score(self):
        """`float`: Bit score of the sequence with all domains after correction.
        """
        assert self._hit != NULL
        return self._hit.score

    @property
    def pre_score(self):
        """`float`: Bit score of the sequence before *null2* correction.
        """
        assert self._hit != NULL
        return self._hit.pre_score

    @property
    def sum_score(self):
        """`float`: Bit score reconstructed from the sum of domain envelopes.

        .. versionadded:: 0.4.6

        """
        assert self._hit != NULL
        return self._hit.sum_score

    @property
    def bias(self):
        """`float`: The *null2* contribution to the uncorrected score.
        """
        assert self._hit != NULL
        return self._hit.pre_score - self._hit.score

    @property
    def domains(self):
        """`~pyhmmer.plan7.Domains`: The list of domains aligned to this hit.
        """
        assert self._hit != NULL
        return Domains(self)

    @property
    def best_domain(self):
        """`~pyhmmer.plan7.Domain`: The best scoring domain in this hit.

        .. versionadded:: 0.4.2

        """
        assert self._hit != NULL
        return Domain(self, self._hit.best_domain)

    @property
    def evalue(self):
        """`float`: The e-value of the hit.
        """
        assert self._hit != NULL
        if self.hits._pli.long_targets:
            return exp(self._hit.lnP)
        else:
            return exp(self._hit.lnP) * self.hits._pli.Z

    @property
    def pvalue(self):
        """`float`: The p-value of the bitscore.

        .. versionadded:: 0.4.2

        """
        assert self._hit != NULL
        return exp(self._hit.lnP)

    @property
    def included(self):
        """`bool`: Whether this hit is marked as *included*.

        .. versionadded:: 0.7.0

        """
        assert self._hit != NULL
        return self._hit.flags & p7_hitflags_e.p7_IS_INCLUDED != 0

    @included.setter
    def included(self, bint included):
        assert self._hit != NULL
        if included:
            if (self._hit.flags & p7_hitflags_e.p7_IS_INCLUDED) == 0:
                self.hits._th.nincluded += 1
            if (self._hit.flags & p7_hitflags_e.p7_IS_REPORTED) == 0:
                self.hits._th.nreported += 1
            self._hit.flags |= p7_hitflags_e.p7_IS_INCLUDED | p7_hitflags_e.p7_IS_REPORTED
            self._hit.flags &= ~(p7_hitflags_e.p7_IS_DROPPED | p7_hitflags_e.p7_IS_DUPLICATE)
        else:
            if (self._hit.flags & p7_hitflags_e.p7_IS_INCLUDED) != 0:
                self.hits._th.nincluded -= 1
            self._hit.flags &= (~p7_hitflags_e.p7_IS_INCLUDED)

    @property
    def reported(self):
        """`bool`: Whether this hit is marked as *reported*.

        .. versionadded:: 0.7.0

        """
        assert self._hit != NULL
        return self._hit.flags & p7_hitflags_e.p7_IS_REPORTED != 0

    @reported.setter
    def reported(self, bint reported):
        assert self._hit != NULL
        if reported:
            if (self._hit.flags & p7_hitflags_e.p7_IS_REPORTED) == 0:
                self.hits._th.nreported += 1
            self._hit.flags |= p7_hitflags_e.p7_IS_REPORTED
        else:
            if (self._hit.flags & p7_hitflags_e.p7_IS_INCLUDED) != 0:
                self.hits._th.nincluded -= 1
            if (self._hit.flags & p7_hitflags_e.p7_IS_REPORTED) != 0:
                self.hits._th.nreported -= 1
            self._hit.flags &= ~(p7_hitflags_e.p7_IS_REPORTED | p7_hitflags_e.p7_IS_INCLUDED)

    @property
    def new(self):
        """`bool`: Whether this hit is marked as *new*.

        .. versionadded:: 0.7.0

        """
        assert self._hit != NULL
        return self._hit.flags & p7_hitflags_e.p7_IS_NEW != 0

    @new.setter
    def new(self, bint new):
        assert self._hit != NULL
        if new:
            self._hit.flags |= p7_hitflags_e.p7_IS_NEW
        else:
            self._hit.flags &= ~p7_hitflags_e.p7_IS_NEW

    @property
    def dropped(self):
        """`bool`: Whether this hit is marked as *dropped*.

        .. versionadded:: 0.7.0

        """
        assert self._hit != NULL
        return self._hit.flags & p7_hitflags_e.p7_IS_DROPPED != 0

    @dropped.setter
    def dropped(self, bint dropped):
        assert self._hit != NULL
        if dropped:
            if (self._hit.flags & p7_hitflags_e.p7_IS_INCLUDED) != 0:
                self.hits._th.nincluded -= 1
            self._hit.flags |= p7_hitflags_e.p7_IS_DROPPED
            self._hit.flags &= ~p7_hitflags_e.p7_IS_INCLUDED
        else:
            self._hit.flags &= ~p7_hitflags_e.p7_IS_DROPPED

    @property
    def duplicate(self):
        """`bool`: Whether this hit is marked as *duplicate*.

        .. versionadded:: 0.7.0

        """
        assert self._hit != NULL
        return self._hit.flags & p7_hitflags_e.p7_IS_DUPLICATE != 0

    @duplicate.setter
    def duplicate(self, bint duplicate):
        assert self._hit != NULL
        if duplicate:
            if (self._hit.flags & p7_hitflags_e.p7_IS_REPORTED) != 0:
                self.hits._th.nreported -= 1
            if (self._hit.flags & p7_hitflags_e.p7_IS_INCLUDED) != 0:
                self.hits._th.nincluded -= 1
            self._hit.flags |= p7_hitflags_e.p7_IS_DUPLICATE
            self._hit.flags &= ~(p7_hitflags_e.p7_IS_INCLUDED | p7_hitflags_e.p7_IS_REPORTED)
        else:
            self._hit.flags &= ~p7_hitflags_e.p7_IS_DUPLICATE


@cython.freelist(8)
cdef class HMM:
    """A data structure storing the Plan7 Hidden Markov Model.

    Attributes:
        alphabet (`~pyhmmer.easel.Alphabet`): The alphabet of the model.

    .. versionchanged:: 0.4.6
       Added the `~HMM.evalue_parameters` and `~HMM.cutoffs` attributes.

    """

    @classmethod
    def sample(
        cls,
        Alphabet alphabet not None,
        int M,
        Randomness randomness not None,
        bint ungapped=False,
        bint enumerable=False,
    ):
        """Sample an HMM of length ``M`` at random.

        Arguments:
            alphabet (`~pyhmmer.easel.Alphabet`): The alphabet of the model.
            M (`int`): The length of the model to generate (i.e. the
                number of nodes).
            randomness (`~pyhmmer.easel.Randomness`): The random number
                generator to use for sampling.
            ungapped (`bool`): Set to `True` to build an ungapped HMM, i.e.
                an HMM where the :math:`M_n \to M_{n+1}` are all one and the
                remaining transitions are zero. Ignored when ``enumerable``
                is `True`.
            enumerable (`bool`): Set to `True` to build a random HMM with no
                nonzero insertion transitions.

        Returns:
            `~pyhmmer.plan7.HMM`: A new HMM generated at random.

        .. versionadded:: 0.7.0

        """
        cdef str fname
        cdef int status
        cdef HMM hmm    = cls.__new__(cls)

        if enumerable:
            fname = "p7_hmm_SampleEnumerable"
            status = libhmmer.p7_hmm.p7_hmm_SampleEnumerable(
                randomness._rng,
                M,
                alphabet._abc,
                &hmm._hmm
            )
        elif ungapped:
            fname = "p7_hmm_SampleUngapped"
            status = libhmmer.p7_hmm.p7_hmm_SampleUngapped(
                randomness._rng,
                M,
                alphabet._abc,
                &hmm._hmm
            )
        else:
            fname = "p7_hmm_Sample"
            status = libhmmer.p7_hmm.p7_hmm_Sample(
                randomness._rng,
                M,
                alphabet._abc,
                &hmm._hmm
            )
        if status != libeasel.eslOK:
            raise UnexpectedError(status, fname)

        hmm.alphabet = alphabet
        hmm._hmm.max_length = M*4 if alphabet.is_nucleotide() else -1
        return hmm


    # --- Magic methods ------------------------------------------------------

    cdef void _initialize(self) noexcept nogil:
        cdef int i
        cdef int j
        cdef int K = self.alphabet._abc.K
        for i in range(self._hmm.M+1):
            libeasel.vec.esl_vec_FSet(self._hmm.mat[i], K, 0.0)
            libeasel.vec.esl_vec_FSet(self._hmm.ins[i], K, 0.0)
            libeasel.vec.esl_vec_FSet(self._hmm.t[i], p7H_NTRANSITIONS, 0.0)
            self._hmm.mat[i][0] = 1.0
            self._hmm.ins[i][0] = 1.0
            self._hmm.t[i][<int> p7h_transitions_e.p7H_MM] = 1.0
            self._hmm.t[i][<int> p7h_transitions_e.p7H_IM] = 1.0
            self._hmm.t[i][<int> p7h_transitions_e.p7H_DM] = 1.0

    cdef void _set_annotation_line(self, str line, char** ptr, int flag) except *:
        assert self._hmm != NULL

        cdef bytes _line
        cdef const unsigned char[::1] view

        if line is None:
            self._hmm.flags &= (~flag)
            free(ptr[0])
            ptr[0] = NULL
        else:
            if len(line) != self._hmm.M:
                raise ValueError(f"Annotation line must be of length {self._hmm.M}, got {len(line)}")
            if ptr[0] == NULL:
                ptr[0] = <char*> calloc(self._hmm.M + 2, sizeof(char))
                if ptr[0] == NULL:
                    raise AllocationError("char", sizeof(char), self._hmm.M + 2)
                ptr[0][0] = ord(' ')
                ptr[0][self._hmm.M + 1] = 0
            self._hmm.flags |= flag
            _line = line.encode()
            view = _line
            memcpy(&ptr[0][1], &view[0], self._hmm.M * sizeof(char))

    def __cinit__(self):
        self.alphabet = None
        self._hmm = NULL

    def __init__(self, Alphabet alphabet not None, int M, bytes name not None):
        """__init__(self, alphabet, M, name)\n--\n

        Create a new HMM from scratch.

        Arguments:
            alphabet (`~pyhmmer.easel.Alphabet`): The alphabet of the model.
            M (`int`): The length of the model (i.e. the number of nodes).
            name (`bytes`): The name of the model.

        """
        # store the alphabet so it's not deallocated
        self.alphabet = alphabet
        # create a new HMM suitable for at least M nodes
        with nogil:
            self._hmm = libhmmer.p7_hmm.p7_hmm_Create(M, alphabet._abc)
        if not self._hmm:
            raise AllocationError("P7_HMM", sizeof(P7_HMM))

        # initialize HMM with arbitrary probabilities so that the HMM
        # is valid and doesn't crash if used as-is (see #36).
        with nogil:
            self._initialize()
        self.name = name

    def __dealloc__(self):
        libhmmer.p7_hmm.p7_hmm_Destroy(self._hmm)

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"<{ty} alphabet={self.alphabet!r} M={self.M!r} name={self.name!r}>"

    def __eq__(self, object other):
        assert self._hmm != NULL

        cdef HMM other_
        cdef int status

        if not isinstance(other, HMM):
            return NotImplemented

        other_ = <HMM> other
        with nogil:
            status = libhmmer.p7_hmm.p7_hmm_Compare(self._hmm, other_._hmm, 0.0)
        if status == libeasel.eslOK:
            return True
        elif status == libeasel.eslFAIL:
            return False
        else:
            raise UnexpectedError(status, "p7_profile_Compare")

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        cdef HMM new
        if id(self) not in memo:
            new = memo[id(self)] = self.copy()
            new.alphabet = copy.deepcopy(self.alphabet, memo=memo)
            new._hmm.abc = new.alphabet._abc
        return memo[id(self)]

    def __sizeof__(self):
        assert self._hmm != NULL
        assert self.alphabet is not None
        assert self.alphabet._abc != NULL

        cdef size_t s = 0
        # level 1
        s += (self._hmm.M + 1) * sizeof(float*)  # hmm->t
        s += (self._hmm.M + 1) * sizeof(float*)  # hmm->mat
        s += (self._hmm.M + 1) * sizeof(float*)  # hmm->inst
        # level 2
        s += p7H_NTRANSITIONS * (self._hmm.M + 1) * sizeof(float)
        s += self.alphabet._abc.K * (self._hmm.M + 1) * sizeof(float)
        s += self.alphabet._abc.K * (self._hmm.M + 1) * sizeof(float)
        # optional fields
        if self._hmm.flags & libhmmer.p7_hmm.p7H_RF:
            s += (self._hmm.M + 2) * sizeof(char)
        if self._hmm.flags & libhmmer.p7_hmm.p7H_MMASK:
            s += (self._hmm.M + 2) * sizeof(char)
        if self._hmm.flags & libhmmer.p7_hmm.p7H_CONS:
            s += (self._hmm.M + 2) * sizeof(char)
        if self._hmm.flags & libhmmer.p7_hmm.p7H_CS:
            s += (self._hmm.M + 2) * sizeof(char)
        if self._hmm.flags & libhmmer.p7_hmm.p7H_CA:
            s += (self._hmm.M + 2) * sizeof(char)
        if self._hmm.flags & libhmmer.p7_hmm.p7H_MAP:
            s += (self._hmm.M + 1) * sizeof(int)
        # annotations
        if self._hmm.name != NULL:
            s += strlen(self._hmm.name) * sizeof(char)
        if self._hmm.acc != NULL:
            s += strlen(self._hmm.acc) * sizeof(char)
        if self._hmm.desc != NULL:
            s += strlen(self._hmm.desc) * sizeof(char)
        if self._hmm.comlog != NULL:
            s += strlen(self._hmm.comlog) * sizeof(char)
        if self._hmm.ctime != NULL:
            s += strlen(self._hmm.ctime) * sizeof(char)
        # struct size
        s += sizeof(P7_HMM)
        s += sizeof(self)
        return s

    def __reduce__(self):
        return HMM, (self.alphabet, self.M, self.name), self.__getstate__()

    def __getstate__(self):
        assert self._hmm != NULL

        cdef dict   state
        cdef object map_

        state = {
            "M": self._hmm.M,
            "t": self.transition_probabilities,
            "mat": self.match_emissions,
            "ins": self.insert_emissions,
            "name": self.name,
            "acc": self.accession,
            "desc": self.description,
            "ctime": self.creation_time,
            "comlog": self.command_line,
            "nseq": self._hmm.nseq,
            "eff_nseq": self._hmm.eff_nseq,
            "max_length": self._hmm.max_length,
            "checksum": self._hmm.checksum,
            "evparam": self.evalue_parameters.as_vector(),
            "cutoff": self.cutoffs.as_vector(),
            "compo": self.composition,
            "offset": self._hmm.offset,
            "alphabet": self.alphabet,
            "flags": self._hmm.flags,
            "consensus": self.consensus,
            "consensus_structure": self.consensus_structure,
            "consensus_accessibility": self.consensus_accessibility,
            "model_mask": self.model_mask,
            "reference": self.reference,
        }

        # copy alignment map if available
        if self._hmm.flags & libhmmer.p7_hmm.p7H_MAP:
            assert self._hmm.map != NULL
            state["map"] = map_ = array.array('i')
            for i in range(self._hmm.M + 1):
                map_.append(self._hmm.map[i])

        return state

    def __setstate__(self, dict state):
        assert self._hmm != NULL

        cdef int                 i
        cdef int                 M         = self._hmm.M
        cdef int                 K         = self._hmm.abc.K
        cdef const float[:, ::1] t         = state["t"]
        cdef const float[:, ::1] mat       = state["mat"]
        cdef const float[:, ::1] ins       = state["ins"]
        cdef const float[::1]    evparam   = state["evparam"]
        cdef const float[::1]    cutoff    = state["cutoff"]
        cdef const float[::1]    compo
        cdef const int[::1]      map_

        # check the HMM has the right dimensions
        if self._hmm.M != state["M"]:
            raise ValueError(f"HMM has a different node count ({self._hmm.M}, state has {state['M']})")
        if self.alphabet != state["alphabet"]:
            raise AlphabetMismatch(self.alphabet, state["alphabet"])
        self._hmm.flags = state["flags"]

        # attributes settable via a property
        self.name = state["name"]
        self.accession = state["acc"]
        self.description = state["desc"]
        self.command_line = state["comlog"]
        self.creation_time = state["ctime"]
        self.consensus = state["consensus"]
        self.consensus_structure = state["consensus_structure"]
        self.consensus_accessibility = state["consensus_accessibility"]
        self.model_mask = state["model_mask"]
        self.reference = state["reference"]

        # numerical values that can be set directly on the C struct
        self._hmm.nseq = state["nseq"]
        self._hmm.eff_nseq = state["eff_nseq"]
        self._hmm.max_length = state["max_length"]
        self._hmm.checksum = state["checksum"]
        self._hmm.offset = state["offset"]

        # vector data that must be copied
        assert self._hmm.t != NULL
        assert self._hmm.t[0] != NULL
        assert self._hmm.ins != NULL
        assert self._hmm.ins[0] != NULL
        assert self._hmm.mat != NULL
        assert self._hmm.mat[0] != NULL
        assert t.ndim == 2
        assert t.shape[0] == M+1
        assert t.shape[1] == p7H_NTRANSITIONS
        assert ins.ndim == 2
        assert ins.shape[0] == M+1
        assert ins.shape[1] == K
        assert mat.ndim == 2
        assert mat.shape[0] == M+1
        assert mat.shape[1] == K
        for i in range(M+1):
            memcpy(self._hmm.ins[i],  &ins[i][0], K                * sizeof(float))
            memcpy(self._hmm.mat[i],  &mat[i][0], K                * sizeof(float))
            memcpy(self._hmm.t[i],    &t[i][0],   p7H_NTRANSITIONS * sizeof(float))

        # arrays that must be copied (compo may be None in the stat dictionary)
        assert evparam.ndim == 1
        assert evparam.shape[0] == p7_NEVPARAM
        assert cutoff.ndim == 1
        assert cutoff.shape[0] == p7_NCUTOFFS
        memcpy(&self._hmm.evparam[0], &evparam[0], p7_NEVPARAM * sizeof(float))
        memcpy(&self._hmm.cutoff[0], &cutoff[0], p7_NCUTOFFS * sizeof(float))
        if self._hmm.flags & libhmmer.p7_hmm.p7H_COMPO:
            compo = state["compo"]
            assert compo.ndim == 1
            assert compo.shape[0] == K
            memcpy(&self._hmm.compo[0], &compo[0], K * sizeof(float))

        # copy alignment map only if it is available
        if self._hmm.flags & p7H_MAP == 0:
            if self._hmm.map != NULL:
                free(self._hmm.map)
                self._hmm.map = NULL
        else:
            map_ = state["map"]
            self._hmm.map = <int*> realloc(self._hmm.map, (M + 1) * sizeof(int))
            if self._hmm.map == NULL:
                raise AllocationError("int", sizeof(int), M+1)
            assert map_.ndim == 1
            assert map_.shape[0] == M + 1
            memcpy(&self._hmm.map[0], &map_[0], (M + 1) * sizeof(int))

    # --- Properties ---------------------------------------------------------

    @property
    def M(self):
        """`int`: The length of the model (i.e. the number of nodes).
        """
        assert self._hmm != NULL
        return self._hmm.M

    @property
    def name(self):
        """`bytes`: The name of the HMM.
        """
        assert self._hmm != NULL
        return None if self._hmm.name == NULL else <bytes> self._hmm.name

    @name.setter
    def name(self, bytes name not None):
        assert self._hmm != NULL

        cdef int   length = len(name)
        cdef char* name_  = <char*> name
        cdef int   err    = libhmmer.p7_hmm.p7_hmm_SetName(self._hmm, name_)

        if err == libeasel.eslEMEM:
            raise AllocationError("char", sizeof(char), length)
        elif err != libeasel.eslOK:
            raise UnexpectedError(err, "p7_hmm_SetName")

    @property
    def accession(self):
        """`bytes` or `None`: The accession of the HMM, if any.
        """
        assert self._hmm != NULL
        return None if self._hmm.acc == NULL else <bytes> self._hmm.acc

    @accession.setter
    def accession(self, bytes accession):
        assert self._hmm != NULL

        cdef char* acc    = NULL if accession is None else <char*> accession
        cdef int   err    = libhmmer.p7_hmm.p7_hmm_SetAccession(self._hmm, acc)
        cdef int   length = 0 if accession is None else len(accession)

        if err == libeasel.eslEMEM:
            raise AllocationError("char", sizeof(char), length)
        elif err != libeasel.eslOK:
            raise UnexpectedError(err, "p7_hmm_SetAccession")

    @property
    def checksum(self):
        """`int` or `None`: The 32-bit checksum of the HMM, if any.

        The checksum if calculated from the alignment the HMM was created
        from, and was introduced in more recent HMM formats. This means
        some `HMM` objects may have a non-`None` checksum.

        .. versionadded:: 0.2.1

        .. versionchanged:: 0.3.1
           Returns `None` if the HMM flag for the checksum is not set.

        """
        assert self._hmm != NULL
        if not (self._hmm.flags & libhmmer.p7_hmm.p7H_CHKSUM):
            return None
        return self._hmm.checksum

    @property
    def composition(self):
        """`~pyhmmer.easel.VectorF` or `None`: The model composition.

        May not be available for freshly-created HMMs. To get the mean
        residue composition emitted by the model, the `~HMM.set_composition`
        method must be called to compute the composition from occupancy.

        Note:
            Although the allocated buffer in the ``P7_HMM`` object is always
            the same dimension so that it can store the largest possible
            alphabet, we only expose the relevant residues. This means that
            the vector will be of size ``alphabet.K``:

                >>> dna = easel.Alphabet.dna()
                >>> dna.K
                4
                >>> hmm = plan7.HMM(dna, 100, b"test")
                >>> hmm.set_composition()
                >>> len(hmm.composition)
                4

        .. versionadded:: 0.4.0

        """
        assert self._hmm != NULL

        cdef VectorF comp

        if not (self._hmm.flags & libhmmer.p7_hmm.p7H_COMPO):
            return None
        comp = VectorF.__new__(VectorF)
        comp._data = &(self._hmm.compo[0])
        comp._owner = self
        comp._n = comp._shape[0] = self.alphabet.K
        return comp

    @property
    def consensus(self):
        """`str` or `None`: The consensus residue line of the HMM, if set.

        .. versionadded:: 0.3.0

        """
        assert self._hmm != NULL
        if not (self._hmm.flags & libhmmer.p7_hmm.p7H_CONS):
            return None
        assert self._hmm.consensus != NULL
        return PyUnicode_FromStringAndSize(&self._hmm.consensus[1], self._hmm.M)

    @consensus.setter
    def consensus(self, str consensus):
        self._set_annotation_line(consensus, &self._hmm.consensus, libhmmer.p7_hmm.p7H_CONS)

    @property
    def consensus_structure(self):
        """`str` or `None`: The consensus structure of the HMM, if any.

        .. versionadded:: 0.3.1

        """
        assert self._hmm != NULL
        if not (self._hmm.flags & libhmmer.p7_hmm.p7H_CS):
            return None
        assert self._hmm.cs != NULL
        return PyUnicode_FromStringAndSize(&self._hmm.cs[1], self._hmm.M)

    @consensus_structure.setter
    def consensus_structure(self, str cs):
        self._set_annotation_line(cs, &self._hmm.cs, libhmmer.p7_hmm.p7H_CS)

    @property
    def consensus_accessibility(self):
        """`str` or `None`: The consensus accessibility of the HMM, if any.

        .. versionadded:: 0.3.1

        """
        assert self._hmm != NULL
        if not (self._hmm.flags & libhmmer.p7_hmm.p7H_CA):
            return None
        assert self._hmm.ca != NULL
        return PyUnicode_FromStringAndSize(&self._hmm.ca[1], self._hmm.M)

    @consensus_accessibility.setter
    def consensus_accessibility(self, str ca):
        self._set_annotation_line(ca, &self._hmm.ca, libhmmer.p7_hmm.p7H_CA)

    @property
    def reference(self):
        """`str` or `None`: The reference line from the alignment, if any.

        This is relevant if the HMM was built from a multiple sequence
        alignment (e.g. by `Builder.build_msa`, or by an external
        ``hmmbuild`` pipeline run).

        .. versionadded:: 0.3.1

        """
        assert self._hmm != NULL
        if not (self._hmm.flags & libhmmer.p7_hmm.p7H_RF):
            return None
        assert self._hmm.rf != NULL
        return PyUnicode_FromStringAndSize(&self._hmm.rf[1], self._hmm.M)

    @reference.setter
    def reference(self, str rf):
        self._set_annotation_line(rf, &self._hmm.rf, libhmmer.p7_hmm.p7H_RF)

    @property
    def model_mask(self):
        """`str` or `None`: The model mask line from the alignment, if any.

        .. versionadded:: 0.3.1

        """
        assert self._hmm != NULL
        if not (self._hmm.flags & libhmmer.p7_hmm.p7H_MM):
            return None
        assert self._hmm.mm != NULL
        return PyUnicode_FromStringAndSize(&self._hmm.mm[1], self._hmm.M)

    @model_mask.setter
    def model_mask(self, str mm):
        self._set_annotation_line(mm, &self._hmm.mm, libhmmer.p7_hmm.p7H_MM)

    @property
    def description(self):
        """`bytes` or `None`: The description of the HMM, if any.
        """
        assert self._hmm != NULL
        return None if self._hmm.desc == NULL else <bytes> self._hmm.desc

    @description.setter
    def description(self, bytes description):
        assert self._hmm != NULL

        cdef char* desc   = NULL if description is None else <char*> description
        cdef int   err    = libhmmer.p7_hmm.p7_hmm_SetDescription(self._hmm, desc)
        cdef int   length = 0 if description is None else len(description)

        if err == libeasel.eslEMEM:
            raise AllocationError("char", sizeof(char), length)
        elif err != libeasel.eslOK:
            raise UnexpectedError(err, "p7_hmm_SetDescription")

    @property
    def transition_probabilities(self):
        r"""`~pyhmmer.easel.MatrixF`: The transition probabilities of the model.

        The property exposes a matrix of shape :math:`(M+1, 7)`, with one row
        per node (plus one extra row for the entry probabilities), and one
        column per transition.

        Columns correspond to the following transitions, in order:
        :math:`M_n \to M_{n+1}`, :math:`M_n \to I_n`, :math:`M_n \to D_{n+1}`,
        :math:`I_n \to I_n`, :math:`I_n \to M_{n+1}`, :math:`D_n \to D_{n+1}`,
        :math:`D_n \to M_{n+1}`. Use the `pyhmmer.plan7.Transitions` enum
        instead of hardcoded indices to make your code more legible.

        Example:
            >>> t = thioesterase.transition_probabilities
            >>> t[1, Transitions.MM]
            0.999...
            >>> t[0, Transitions.DM]  # 1 by convention for the first node
            1.0
            >>> t[-1, Transitions.MD] # 0 by convention for the last node
            0.0

        Caution:
            If editing this matrix manually, note that some invariants need
            to hold for the HMM to be valid: :math:`I_n`, :math:`M_n` and
            :math:`D_n` transition probabilities should only contain
            probabilities between 0 and 1, and sum to 1::

                >>> t = thioesterase.transition_probabilities
                >>> t[50, Transitions.MM] + t[50, Transitions.MI] + t[50, Transitions.MD]
                1.000...
                >>> t[50, Transitions.IM] + t[50, Transitions.II]
                1.000...
                >>> t[50, Transitions.DM] + t[50, Transitions.DD]
                1.000...

            Consider calling `HMM.validate` after manual edition.

        .. versionadded:: 0.3.1

        .. versionchanged:: 0.4.0
            This property is now a `~pyhmmer.easel.MatrixF`.

        """
        assert self._hmm != NULL
        cdef MatrixF mat = MatrixF.__new__(MatrixF)
        mat._m = mat._shape[0] = self._hmm.M + 1
        mat._n = mat._shape[1] = libhmmer.p7_hmm.p7H_NTRANSITIONS
        mat._owner = self
        mat._data = <void**> self._hmm.t
        return mat

    @property
    def match_emissions(self):
        """`~pyhmmer.easel.MatrixF`: The match emissions of the model.

        The property exposes a matrix of shape :math:`(M+1, K)`, with one
        row per node and one column per alphabet symbol.

        Note:
            Since the first row corresponds to the entry probabilities, the
            emissions are unused. By convention, it should still contain
            valid probabilities, so it will always be set as follow with
            1 probability for the first symbol, and 0 for the rest::

                >>> hmm = HMM(easel.Alphabet.dna(), 100, b"test")
                >>> hmm.match_emissions
                MatrixF([[1.0, 0.0, 0.0, 0.0], ...])

        Caution:
            If editing this matrix manually, note that rows must contain
            valid probabilities for the HMM to be valid: each row must
            contains values between 0 and 1, and sum to 1. Consider
            calling `HMM.validate` after manual edition.

        .. versionadded:: 0.3.1

        .. versionchanged:: 0.4.0
            This property is now a `~pyhmmer.easel.MatrixF`, and stores row 0.

        """
        assert self._hmm != NULL
        cdef MatrixF mat = MatrixF.__new__(MatrixF)
        mat._m = mat._shape[0] = self._hmm.M + 1
        mat._n = mat._shape[1] = self.alphabet.K
        mat._owner = self
        mat._data = <void**> self._hmm.mat
        return mat

    @property
    def insert_emissions(self):
        """`~pyhmmer.easel.MatrixF`: The insert emissions of the model.

        The property exposes a matrix of shape :math:`(M+1, K)`, with one
        row per node and one column per alphabet symbol.

        Caution:
            If editing this matrix manually, note that rows must contain
            valid probabilities for the HMM to be valid: each row must
            contains values between 0 and 1, and sum to 1. Consider
            calling `HMM.validate` after manual edition.

        .. versionadded:: 0.3.1

        .. versionchanged:: 0.4.0
            This property is now a `~pyhmmer.easel.MatrixF`, and stores row 0.

        """
        assert self._hmm != NULL
        cdef MatrixF mat = MatrixF.__new__(MatrixF)
        mat._m = mat._shape[0] = self._hmm.M + 1
        mat._n = mat._shape[1] = self.alphabet.K
        mat._owner = self
        mat._data = <void**> self._hmm.ins
        return mat

    @property
    def command_line(self):
        """`str` or `None`: The command line that built the model.

        For HMMs created with `~pyhmmer.plan7.Builder`, this defaults to
        `sys.argv`. It can however be set to any string, including multiline
        to show successive commands.

        Example:
            >>> print(thioesterase.command_line)
            hmmbuild Thioesterase.hmm Thioesterase.fa
            hmmcalibrate Thioesterase.hmm

        .. versionadded:: 0.3.1

        """
        assert self._hmm != NULL
        if self._hmm.comlog == NULL:
            return None
        return self._hmm.comlog.decode("ascii")

    @command_line.setter
    def command_line(self, object cli):
        assert self._hmm != NULL

        if cli is None:
            free(self._hmm.comlog)
            self._hmm.comlog = NULL
            return

        cdef bytes  cli_ = cli.encode("ascii")
        cdef size_t n    = strlen(cli_)

        if self._hmm.comlog == NULL:
            self._hmm.comlog = strndup(<const char*> cli_, n + 1)
        else:
            self._hmm.comlog = <char*> realloc(<void*> self._hmm.comlog, sizeof(char) * (n + 1))
            if self._hmm.comlog != NULL:
                strncpy(self._hmm.comlog, <char*> cli_, n+1)
        if self._hmm.comlog == NULL:
            raise AllocationError("char", sizeof(char), n+1)

    @property
    def nseq(self):
        """`int` or `None`: The number of training sequences used, if any.

        If the HMM was created from a multiple sequence alignment, this
        corresponds to the number of sequences in the MSA.

        Example:
            >>> thioesterase.nseq
            278

        .. versionadded:: 0.3.1

        """
        assert self._hmm != NULL
        return None if self._hmm.nseq == -1 else self._hmm.nseq

    @property
    def nseq_effective(self):
        """`float` or `None`: The number of effective sequences used, if any.

        .. versionadded:: 0.3.1

        """
        assert self._hmm != NULL
        return None if self._hmm.eff_nseq == -1 else self._hmm.eff_nseq

    @property
    def creation_time(self):
        """`datetime.datetime` or `None`: The creation time of the HMM, if any.

        Example:
            Get the creation time for any HMM::

                >>> thioesterase.creation_time
                datetime.datetime(2008, 11, 25, 17, 28, 32)

            Set the creation time manually to a different date and time::

                >>> ctime = datetime.datetime(2021, 8, 23, 23, 59, 19)
                >>> thioesterase.creation_time = ctime
                >>> thioesterase.creation_time
                datetime.datetime(2021, 8, 23, 23, 59, 19)

        Danger:
            Internally, ``libhmmer`` always uses ``asctime`` to generate a
            timestamp for the HMMs, so this property assumes that every
            creation time field can be parsed into a `datetime.datetime`
            object using the  ``"%a %b %d %H:%M:%S %Y"`` format.

        .. versionadded:: 0.4.6

        """
        assert self._hmm != NULL

        cdef size_t l
        cdef str    ctime

        if self._hmm.ctime == NULL:
            return None

        l = strlen(self._hmm.ctime)
        ctime = PyUnicode_DecodeASCII(self._hmm.ctime, l, NULL)
        return datetime.datetime.strptime(ctime,'%a %b %d %H:%M:%S %Y')

    @creation_time.setter
    def creation_time(self, object ctime):
        assert self._hmm != NULL

        cdef str    ty
        cdef bytes  formatted
        cdef size_t n

        if ctime is None:
            free(self._hmm.ctime)
            self._hmm.ctime = NULL
            return
        elif not isinstance(ctime, datetime.datetime):
            ty = type(ctime).__name__
            raise TypeError(f"Expected datetime.datetime or None, found {ty}")

        formatted = ctime.strftime('%a %b %e %H:%M:%S %Y').encode('ascii')
        n = len(formatted)

        if self._hmm.ctime == NULL:
            self._hmm.ctime = <char*> malloc(sizeof(char) * (n + 1))
        else:
            self._hmm.ctime = <char*> realloc(<void*> self._hmm.ctime, sizeof(char) * (n + 1))
        if self._hmm.ctime == NULL:
            raise AllocationError("char", sizeof(char), n+1)
        if self._hmm.ctime != NULL:
            strncpy(self._hmm.ctime, <const char*> formatted, n + 1)

    @property
    def evalue_parameters(self):
        """`EvalueParameters`: The e-value parameters for this HMM.
        """
        assert self._hmm != NULL
        cdef EvalueParameters ep = EvalueParameters.__new__(EvalueParameters)
        ep._evparams = &self._hmm.evparam
        ep._owner = self
        return ep

    @property
    def cutoffs(self):
        """`~plan7.Cutoffs`: The bitscore cutoffs for this HMM.
        """
        assert self._hmm != NULL
        cdef Cutoffs cutoffs = Cutoffs.__new__(Cutoffs)
        cutoffs._owner = self
        cutoffs._cutoffs = &self._hmm.cutoff
        cutoffs._flags = &self._hmm.flags
        cutoffs._is_profile = False
        return cutoffs

    # --- Methods ------------------------------------------------------------

    cpdef HMM copy(self):
        """Return a copy of the HMM with the exact same configuration.

        .. versionadded:: 0.3.0

        """
        assert self._hmm != NULL

        cdef HMM new = HMM.__new__(HMM)
        new.alphabet = self.alphabet

        with nogil:
            new._hmm = libhmmer.p7_hmm.p7_hmm_Clone(self._hmm)
        if new._hmm == NULL:
            raise AllocationError("P7_HMM", sizeof(P7_HMM))
        return new

    cpdef VectorF match_occupancy(self):
        """Calculate the match occupancy for each match state.

        Returns:
            `~pyhmmer.easel.VectorF`: A vector of size :math:`M+1`
            containing the probability that each match state is used
            in a sample glocal path through the model.

        .. versionadded:: 0.4.10

        """
        assert self._hmm != NULL
        cdef VectorF mocc = VectorF.zeros(self._hmm.M)
        with nogil:
            status = libhmmer.p7_hmm.p7_hmm_CalculateOccupancy(
                self._hmm,
                <float*> mocc._data,
                NULL
            )
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_hmm_CalculateOccupancy")
        return mocc

    cpdef double mean_match_entropy(self) except *:
        """Compute the mean entropy per HMM match state.

        The mean entropy per match state emission distribution is defined
        as:

        .. math::

             - \\frac{1}{M} \\sum_{k=1}^{M} \\sum_x p_k(x) \\log_2 p_k(x)

        where :math:`p_k(x)` is the emission probability for symbol
        :math:`x` from match state :math:`k`.

        Returns:
            `float`: The mean entropy, in bits.

        Example:
            >>> thioesterase.mean_match_entropy()
            3.0425...

        .. versionadded:: 0.4.10

        """
        assert self._hmm != NULL
        with nogil:
            return libhmmer.modelstats.p7_MeanMatchEntropy(self._hmm)

    cpdef double mean_match_information(self, Background background) except *:
        """Compute the mean information content of the HMM match states.

        The mean information content per match state emission distribution
        is defined as:

        .. math::

            \\frac{1}{M} \sum_{k=1}^{M}
                \\left[
                    \\sum_x f(x)   \\log_2 f(x)
                  - \\sum_x p_k(x) \log_2 p_k(x)
                \\right]

        in bits, where :math:`p_k(x)` is the emission probability for symbol
        :math:`x` from match state :math:`k`, and :math:`f(x)` is the
        background emission probability for :math:`x` from the null model.

        Arguments:
            background (`~plan7.Background`): The null model from which to
                get the background emission probabilities.

        Returns:
            `float`: The mean information content, in bits.

        Example:
            >>> background = plan7.Background(easel.Alphabet.amino())
            >>> thioesterase.mean_match_information(background)
            1.1330...

        .. versionadded:: 0.4.10

        """
        assert self._hmm != NULL
        assert background._bg != NULL
        with nogil:
            return libhmmer.modelstats.p7_MeanMatchInfo(self._hmm, background._bg)

    cpdef double mean_match_relative_entropy(self, Background background) except *:
        """Compute the mean relative entropy per HMM match state.

        The mean relative entropy per match state emission distribution
        is defined as:

        .. math::

            \\frac{1}{M} \\sum_{k=1}^{M} \\sum_x {
                p_k(x) \\log_2 \\frac{p_k(x)}{f(x)}
            }

        in bits, where :math:`p_k(x)` is the emission probability for symbol
        :math:`x` from match state :math:`k`, and :math:`f(x)` is the
        background emission probability for :math:`x` from the null model.

        Arguments:
            background (`~plan7.Background`): The null model from which to
                get the background emission probabilities.

        Returns:
            `float`: The mean relative entropy, in bits.

        Example:
            >>> background = plan7.Background(easel.Alphabet.amino())
            >>> thioesterase.mean_match_relative_entropy(background)
            1.1201...

        .. versionadded:: 0.4.10

        """
        assert self._hmm != NULL
        assert background._bg != NULL
        with nogil:
            return libhmmer.modelstats.p7_MeanMatchRelativeEntropy(self._hmm, background._bg)

    cpdef void renormalize(self):
        """Renormalize all parameter vectors (emissions and transitions).

        .. versionadded:: 0.4.0

        """
        assert self._hmm != NULL

        cdef int status
        with nogil:
            status = libhmmer.p7_hmm.p7_hmm_Renormalize(self._hmm)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_hmm_Renormalize")

    cpdef void scale(self, double scale, bint exponential=False):
        """Rescale counts by a factor in a model containing counts.

        This method only affects core probability model emissions and
        transitions.

        Arguments:
            scale (`float`): The scaling factor to use (:math:`1.0` for no
                scaling). Often computed using the ratio of effective
                sequences (:math:`\\frac{n_{eff}}{n_{seq}}`)
            exponential (`bool`): When set to `True`, use ``scale`` as an
                exponential factor (:math:`C_i = C_i^s`) instead of
                a multiplicative factor (:math:`C_i = C_i \\times s`),
                resulting in a non-uniform scaling across columns. This
                can be useful when some heavily fragmented sequences are
                used to reconstruct a family MSA.

        .. versionadded:: 0.4.0

        """
        assert self._hmm != NULL

        cdef int status
        with nogil:
            if exponential:
                status = libhmmer.p7_hmm.p7_hmm_ScaleExponential(self._hmm, scale)
            else:
                status = libhmmer.p7_hmm.p7_hmm_Scale(self._hmm, scale)
        if status != libeasel.eslOK:
            func = "p7_hmm_ScaleExponential" if exponential else "p7_hmm_Scale"
            raise UnexpectedError(status, func)

    cpdef void set_composition(self) except *:
        """Calculate and set the model composition.

        .. versionadded:: 0.4.0

        """
        assert self._hmm != NULL

        cdef int status
        with nogil:
            status = libhmmer.p7_hmm.p7_hmm_SetComposition(self._hmm)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_hmm_SetComposition")

    cpdef void set_consensus(self, DigitalSequence sequence = None) except *:
        """Set the model consensus sequence.

        If given a sequence, the HMM is assumed to originate from a
        single-sequence model, and the sequence is used as the consensus.
        Otherwise, the consensus is computed by extracting the consensus
        at each position.

        Residues in the consensus sequence are uppercased when their emission
        probabilities are above an arbitrary, alphabet-specific threshold
        (:math:`0.9` for nucleotide alphabets, :math:`0.5` for protein).

        Arguments:
            sequence (`~pyhmmer.easel.DigitalSequence`): A sequence in
                digital mode with :math:`M` residues and the same
                alphabet as the HMM.

        Raises:
            `ValueError`: When given a sequence with the wrong length.
            `~pyhmmer.errors.AlphabetMismatch`: When given a sequence in
                the wrong alphabet.

        .. versionadded:: 0.10.1

        """
        assert self._hmm != NULL

        cdef int     status
        cdef ESL_SQ* sq     = NULL

        if sequence is not None:
            sq = sequence._sq
            if sequence.alphabet != self.alphabet:
                raise AlphabetMismatch(self.alphabet, sequence.alphabet)
            if sq.n < self.M:
                raise ValueError(f"Expected `DigitalSequence` of length {self.M!r}, found {sq.n!r}")

        with nogil:
            status = libhmmer.p7_hmm.p7_hmm_SetConsensus(self._hmm, sq)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_hmm_SetConsensus")

    cpdef Profile to_profile(
        self,
        Background background = None,
        int L = 400,
        bint multihit = True,
        bint local = True
    ):
        """Create a new profile configured for this HMM.

        This method is a shortcut for creating a new `~pyhmmer.plan7.Profile`
        and calling `~pyhmmer.plan7.Profile.configure` for a given HMM.
        Prefer manually calling `~pyhmmer.plan7.Profile.configure` to recycle
        the profile buffer when running inside a loop.

        Arguments:
            background (`~pyhmmer.plan7.Background`, optional): The null
                background model. In `None` given, create a default one
                for the HMM alphabet.
            L (`int`): The expected target sequence length.
            multihit (`bool`): Whether or not to use multihit mode.
            local (`bool`): Whether to use local or global mode.

        .. versionadded:: 0.8.0

        """
        cdef Profile profile = Profile(alphabet=self.alphabet, M=self.M)
        background = Background(self.alphabet) if background is None else background
        profile.configure(self, background, L=L, multihit=multihit, local=local)
        return profile

    cpdef void validate(self, float tolerance=1e-4) except *:
        """Validate the HMM agains the Plan7 structural constraints.

        HMMs created with PyHMMER are always valid, whether they come from
        a `~pyhmmer.plan7.Builder` or from the `HMM` constructor. However,
        it is possible to manually edit the emission scores and transition
        probabilities, and the structural constrains may not hold. Call
        this method to make sure the HMM is valid after you are done editing
        it.

        Arguments:
            tolerance (`float`): The tolerance with which to check that the
                transition probabilities sum to 1.

        Raises:
            `ValueError`: When the HMM fails validation.

        .. versionadded:: 0.8.1

        """
        cdef char[eslERRBUFSIZE] errbuf

        if libhmmer.p7_hmm.p7_hmm_Validate(self._hmm, errbuf, tolerance) != libeasel.eslOK:
            err_msg = errbuf.decode("utf-8", "replace")
            raise InvalidHMM(self, err_msg)

    cpdef void write(self, object fh, bint binary=False) except *:
        """Write the HMM to a file handle.

        Arguments:
            fh (`io.IOBase`): A Python file handle, opened in binary mode
                (this must be the case even with ``binary=False``, since
                the C code will emit bytes in either case).
            binary (`bool`): Pass ``False`` to emit the file in ASCII mode
                using the latest supported HMMER format, or ``True`` to use
                the binary HMMER3 format.

        """
        cdef PyObject*  type
        cdef PyObject*  value
        cdef PyObject*  traceback
        cdef int        status
        cdef FILE*      file
        cdef P7_HMM*    hm     = self._hmm

        file = fopen_obj(fh, "w")

        if binary:
            status = libhmmer.p7_hmmfile.p7_hmmfile_WriteBinary(file, -1, hm)
        else:
            status = libhmmer.p7_hmmfile.p7_hmmfile_WriteASCII(file, -1, hm)

        if status == libeasel.eslOK:
            fclose(file)
        else:
            _reraise_error()
            raise UnexpectedError(status, "p7_hmmfile_WriteASCII")

    cpdef void zero(self) noexcept:
        """Set all parameters to zero, including model composition.
        """
        assert self._hmm != NULL
        with nogil:
            libhmmer.p7_hmm.p7_hmm_Zero(self._hmm)


cdef class HMMFile:
    """A wrapper around a file (or database), storing serialized HMMs.

    You can use this class to iterate on the HMMs inside a file,
    loading them into `~pyhmmer.plan7.HMM` objects.

    Example:
        Load the first HMM from an HMM file located on the
        local filesystem::

            >>> with HMMFile("tests/data/hmms/txt/PF02826.hmm") as hmm_file:
            ...     hmm = hmm_file.read()
            >>> hmm.name
            b'2-Hacid_dh_C'
            >>> hmm.accession
            b'PF02826.20'

        Load all the HMMs from an HMM file into a `list`::

            >>> with HMMFile("tests/data/hmms/txt/RREFam.hmm") as hmm_file:
            ...     hmms = list(hmm_file)
            >>> len(hmms)
            28
            >>> hmms[0].accession
            b'RREFam008.1'

    """

    _FORMATS = dict(HMM_FILE_FORMATS)
    _MAGIC = dict(HMM_FILE_MAGIC)

    # --- Constructor --------------------------------------------------------

    @staticmethod
    cdef P7_HMMFILE* _open_fileobj(object fh) except *:
        cdef int         status
        cdef char*       token
        cdef int         token_len
        cdef bytes       filename
        cdef object      fh_       = fh
        cdef P7_HMMFILE* hfp       = NULL

        # use buffered IO to be able to peek efficiently
        if not hasattr(fh, "peek"):
            fh_ = io.BufferedReader(fh)

        # attempt to allocate space for the P7_HMMFILE
        hfp = <P7_HMMFILE*> malloc(sizeof(P7_HMMFILE));
        if hfp == NULL:
            raise AllocationError("P7_HMMFILE", sizeof(P7_HMMFILE))

        # store options
        hfp.f            = fopen_obj(fh_, "r")
        hfp.do_gzip      = False
        hfp.do_stdin     = False
        hfp.newly_opened = True
        hfp.is_pressed   = False

        # set pointers as NULL for now
        hfp.parser    = NULL
        hfp.efp       = NULL
        hfp.ffp       = NULL
        hfp.pfp       = NULL
        hfp.ssi       = NULL
        hfp.fname     = NULL
        hfp.errbuf[0] = b"\0"

        # extract the filename if the file handle has a `name` attribute
        if getattr(fh, "name", None) is not None:
            filename = fh.name.encode()
            hfp.fname = strdup(filename)
            if hfp.fname == NULL:
                raise AllocationError("char", sizeof(char), strlen(filename))

        # check if the parser is in binary format,
        magic = int.from_bytes(fh_.peek(4)[:4], sys.byteorder)
        if magic in HMM_FILE_MAGIC:
            hfp.format = HMM_FILE_MAGIC[magic]
            hfp.parser = read_bin30hmm
            # NB: the file must be advanced, since read_bin30hmm assumes
            #     the binary tag has been skipped already, buf we only peeked
            #     so far; note that we advance without seeking or rewinding.
            fh_.read(4)
            return hfp

        # create and configure the file parser
        hfp.efp = libeasel.fileparser.esl_fileparser_Create(hfp.f)
        if hfp.efp == NULL:
            libhmmer.p7_hmmfile.p7_hmmfile_Close(hfp)
            raise AllocationError("ESL_FILEPARSER", sizeof(ESL_FILEPARSER))
        status = libeasel.fileparser.esl_fileparser_SetCommentChar(hfp.efp, b"#")
        if status != libeasel.eslOK:
            libhmmer.p7_hmmfile.p7_hmmfile_Close(hfp)
            raise UnexpectedError(status, "esl_fileparser_SetCommentChar")

        # get the magic string at the beginning
        status = libeasel.fileparser.esl_fileparser_NextLine(hfp.efp)
        if status == libeasel.eslEOF:
            raise EOFError("HMM file is empty")
        elif status != libeasel.eslOK:
            libhmmer.p7_hmmfile.p7_hmmfile_Close(hfp)
            raise UnexpectedError(status, "esl_fileparser_NextLine");
        status = libeasel.fileparser.esl_fileparser_GetToken(hfp.efp, &token, &token_len)
        if status != libeasel.eslOK:
            libhmmer.p7_hmmfile.p7_hmmfile_Close(hfp)
            raise UnexpectedError(status, "esl_fileparser_GetToken");

        # detect the format
        if token.startswith(b"HMMER3/"):
            hfp.parser = read_asc30hmm
            format = token[5:].decode("utf-8", "replace")
            if format in HMM_FILE_FORMATS:
                hfp.format = HMM_FILE_FORMATS[format]
            else:
                hfp.parser = NULL
        elif token.startswith(b"HMMER2.0"):
            hfp.parser = read_asc20hmm
            hfp.format = p7_hmmfile_formats_e.p7_HMMFILE_20

        # check the format tag was recognized
        if hfp.parser == NULL:
            text = token.decode("utf-8", "replace")
            libhmmer.p7_hmmfile.p7_hmmfile_Close(hfp)
            raise ValueError("Unrecognized format tag in HMM file: {!r}".format(text))

        # return the finalized P7_HMMFILE*
        return hfp

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._alphabet = None
        self._hfp = NULL
        self._name = None

    def __init__(self, object file, bint db = True):
        """__init__(self, file, db=True)\n--\n

        Create a new HMM reader from the given file.

        Arguments:
            file (`str`, `bytes`, `os.PathLike` or file-like object): Either
                the path to a file containing the HMMs to read, or a
                file-like object in **binary mode**.
            db (`bool`): Set to `False` to force the parser to ignore the
                pressed HMM database if it finds one. Defaults to `True`.

        """
        cdef int                 status
        cdef bytes               fspath
        cdef char[eslERRBUFSIZE] errbuf

        try:
            fspath = os.fsencode(file)
            self._name = os.fsdecode(fspath)
            if db:
                function = "p7_hmmfile_OpenE"
                status = libhmmer.p7_hmmfile.p7_hmmfile_Open(fspath, NULL, &self._hfp, errbuf)
            else:
                function = "p7_hmmfile_OpenENoDB"
                status = libhmmer.p7_hmmfile.p7_hmmfile_OpenNoDB(fspath, NULL, &self._hfp, errbuf)
        except TypeError:
            self._hfp = HMMFile._open_fileobj(file)
            status    = libeasel.eslOK

        if status == libeasel.eslENOTFOUND:
            raise FileNotFoundError(errno.ENOENT, f"No such file or directory: {file!r}")
        elif status == libeasel.eslEFORMAT:
            if fspath is not None:
                if os.path.isdir(fspath):
                    raise IsADirectoryError(errno.EISDIR, f"Is a directory: {file!r}")
                elif os.stat(file).st_size == 0:
                    raise EOFError("HMM file is empty")
            raise ValueError("format not recognized by HMMER")
        elif status != libeasel.eslOK:
            raise UnexpectedError(status, function)

        self._alphabet = Alphabet.__new__(Alphabet)
        self._alphabet._abc = NULL
        self._file = file

    def __dealloc__(self):
        if self._hfp:
            warnings.warn("unclosed HMM file", ResourceWarning)
            self.close()

    def __repr__(self):
        cdef str ty = type(self).__name__
        if self._name is not None:
            return f"{ty}({self._name!r})"
        else:
            return f"<{ty} file={self._file!r}>"

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def __iter__(self):
        return self

    def __next__(self):
        cdef HMM hmm = self.read()
        if hmm is None:
            raise StopIteration()
        return hmm

    # --- Properties ---------------------------------------------------------

    @property
    def closed(self):
        """`bool`: Whether the `HMMFile` is closed or not.
        """
        return self._hfp == NULL

    @property
    def name(self):
        """`str` or `None`: The path to the HMM file, if known.

        .. versionadded:: 0.7.0

        """
        return self._name

    # --- Python Methods -----------------------------------------------------

    cpdef void rewind(self) except *:
        """Rewind the file back to the beginning.
        """
        cdef int status
        if self._hfp == NULL:
            raise ValueError("I/O operation on closed file.")
        status = libhmmer.p7_hmmfile.p7_hmmfile_Position(self._hfp, 0)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_hmmfile_Position")
        # Manually rewind the MSV readers as well, if the file is pressed
        if self._hfp.is_pressed:
            rewind(self._hfp.ffp)
            rewind(self._hfp.pfp)

    cpdef HMM read(self):
        """Read the next HMM from the file.

        Returns:
            `HMM`: The next HMM in the file, or `None` if all HMMs were read
            from the file already.

        Raises:
            `ValueError`: When attempting to read a HMM from a closed
                file, or when the file could not be parsed.
            `~pyhmmer.errors.AllocationError`: When memory for the HMM could
                not be allocated successfully.
            `~pyhmmer.errors.AlphabetMismatch`: When the file contains HMMs
                in different alphabets, or in an alphabet that is different
                from the alphabet used to initialize the `HMMFile`.

        .. versionadded:: 0.4.11

        """
        cdef int status
        cdef HMM py_hmm
        cdef P7_HMM* hmm = NULL

        if self._hfp == NULL:
            raise ValueError("I/O operation on closed file.")

        # don't run in *nogil* because the file may call a file-like handle
        status = libhmmer.p7_hmmfile.p7_hmmfile_Read(self._hfp, &self._alphabet._abc, &hmm)

        if status == libeasel.eslOK:
            py_hmm = HMM.__new__(HMM)
            py_hmm.alphabet = self._alphabet # keep a reference to the alphabet
            py_hmm._hmm = hmm
            return py_hmm
        elif status == libeasel.eslEOF:
            return None
        elif status == libeasel.eslEMEM:
            raise AllocationError("P7_HMM", sizeof(P7_HMM))
        elif status == libeasel.eslESYS:
            raise OSError(self._hfp.errbuf.decode("utf-8", "replace"))
        elif status == libeasel.eslEFORMAT:
            raise ValueError("Invalid format in file: {}".format(self._hfp.errbuf.decode("utf-8", "replace")))
        elif status == libeasel.eslEINCOMPAT:
            raise AlphabetMismatch(self._alphabet)
        else:
            _reraise_error()
            raise UnexpectedError(status, "p7_hmmfile_Read")

    cpdef void close(self) except *:
        """Close the HMM file and free resources.

        This method has no effect if the file is already closed. It is called
        automatically if the `HMMFile` was used in a context::

            >>> with HMMFile("tests/data/hmms/bin/PKSI-AT.h3m") as hmm_file:
            ...     hmm = hmm_file.read()

        """
        if self._hfp:
            libhmmer.p7_hmmfile.p7_hmmfile_Close(self._hfp)
            self._hfp = NULL

    cpdef bint is_pressed(self) except *:
        """Check whether the HMM file is a pressed HMM database.

        A pressed database is an HMMER format to store optimized profiles
        in several files on the disk. It can be used to reduce the time
        needed to process sequences by cutting down the time needed to
        convert from an `HMM` to an `OptimizedProfile`.

        Example:
            >>> HMMFile("tests/data/hmms/txt/PKSI-AT.hmm").is_pressed()
            False
            >>> HMMFile("tests/data/hmms/bin/PKSI-AT.h3m").is_pressed()
            False
            >>> HMMFile("tests/data/hmms/db/PKSI-AT.hmm").is_pressed()
            True

        .. versionadded:: 0.4.11

        """
        if self._hfp == NULL:
            raise ValueError("I/O operation on closed file.")
        return self._hfp.is_pressed

    cpdef HMMPressedFile optimized_profiles(self):
        """Get an iterator over the `OptimizedProfile` in the HMM database.

        Returns:
            `~pyhmmer.plan7.HMMPressedFile`: An iterator over the optimized
            profiles in a pressed HMM database.

        .. versionadded:: 0.4.11

        """
        if self._hfp == NULL:
            raise ValueError("I/O operation on closed file.")
        if not self._hfp.is_pressed:
            raise ValueError("HMM file does not contain optimized profiles.")
        cdef HMMPressedFile optimized = HMMPressedFile.__new__(HMMPressedFile)
        optimized._alphabet = self._alphabet
        optimized._hmmfile = self
        optimized._hfp = self._hfp
        return optimized


cdef class HMMPressedFile:
    """An iterator over each `OptimizedProfile` in a pressed HMM database.

    This class cannot be instantiated: use the `HMMFile.optimized_profiles`
    to obtain an instance from an `HMMFile` that wraps a pressed HMM
    database.

    Example:
        Use a pressed HMM database to run a search pipeline using optimized
        profiles directly, instead of converting them from the text HMMs::

            >>> with HMMFile("tests/data/hmms/db/Thioesterase.hmm") as hfile:
            ...     models = hfile.optimized_profiles()
            ...     hits = next(pyhmmer.hmmsearch(models, proteins))

        In this example, `~pyhmmer.hmmer.hmmsearch` will receive an iterator
        of `OptimizedProfile` instead of an iterator of `HMM` is ``hfile``
        was passed directly; this lets `~pyhmmer.hmmer.hmmsearch` skip the
        conversion step before running the search pipeline.

    .. versionadded:: 0.4.11

    .. versionchanged:: 0.7.0
        Allow instantiating an `HMMPressedFile` from a filename.

    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._alphabet = None
        self._hmmfile = None
        self._hfp = NULL
        self._position = 0

    def __init__(self, object file):
        """__init__(self, file)\n--\n

        Create a new pressed file from the given filename.

        Arguments:
            file (`str`, `bytes` or `os.PathLike`): The path to the pressed
                HMM file containing the optimized profiles to read.

        """
        self._hmmfile = HMMFile(file, db=True)
        self._hfp = self._hmmfile._hfp
        if not self._hfp.is_pressed:
            raise ValueError("HMM file does not contain optimized profiles.")
        self._alphabet = self._hmmfile._alphabet

    def __iter__(self):
        return self

    def __next__(self):
        cdef OptimizedProfile om = self.read()
        if om is None:
            raise StopIteration()
        return om

    def __repr__(self):
        cdef str ty = type(self).__name__
        cdef str mod = type(self).__module__
        if self._hmmfile._name is not None:
            return f"{ty}({self._hmmfile._name!r})"
        else:
            return f"<{ty} file={self._hmmfile._file!r}>"

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def __len__(self):
        assert self._hfp.ssi != NULL
        return self._hfp.ssi.nprimary - self._position

    # --- Properties ---------------------------------------------------------

    @property
    def closed(self):
        """`bool`: Whether the `HMMPressedFile` is closed or not.

        .. versionadded:: 0.7.0

        """
        return self._hmmfile.closed

    @property
    def name(self):
        """`str`: The path to the HMM file.

        Note:
            Unlike `HMMFile.name` this attribute can never be `None`, since a
            `HMMPressedFile` object can only be created from a filename, and
            not from a file-like object.

        .. versionadded:: 0.7.0

        """
        return self._hmmfile.name


    # --- Methods ------------------------------------------------------------

    cpdef void rewind(self) except *:
        """Rewind the file back to the beginning.
        """
        self._hmmfile.rewind()
        self._position = 0

    cpdef OptimizedProfile read(self):
        """Read the next optimized profile from the file.

        Returns:
            `OptimizedProfile`: The next optimized profile in the file, or
            `None` if all profiles were read from the file already.

        Raises:
            `ValueError`: When attempting to read an optimized profile from a
                closed file, or when the file could not be parsed.
            `~pyhmmer.errors.AllocationError`: When memory for the
                `OptimizedProfile` could not be allocated successfully.
            `~pyhmmer.errors.AlphabetMismatch`: When the file contains
                optimized profiles in different alphabets, or in an alphabet
                that is different from the alphabet used to initialize the
                `HMMFile`.

        .. versionadded:: 0.4.11

        """
        cdef int              status
        cdef OptimizedProfile om     = OptimizedProfile.__new__(OptimizedProfile)

        if self._hfp == NULL:
            raise ValueError("I/O operation on closed file.")

        with nogil:
            status = p7_oprofile_ReadMSV(self._hfp, &self._alphabet._abc, &om._om)
            if status == libeasel.eslOK:
                status = p7_oprofile_ReadRest(self._hfp, om._om)

        if status == libeasel.eslOK:
            om.alphabet = self._alphabet
            self._position += 1
            return om
        elif status == libeasel.eslEOF:
            return None
        elif status == libeasel.eslEMEM:
            raise AllocationError("P7_OPROFILE", sizeof(P7_OPROFILE))
        elif status == libeasel.eslESYS:
            raise OSError(self._hfp.errbuf.decode("utf-8", "replace"))
        elif status == libeasel.eslEFORMAT:
            raise ValueError("Invalid format in file: {}".format(self._hfp.errbuf.decode("utf-8", "replace")))
        elif status == libeasel.eslEINCOMPAT:
            raise AlphabetMismatch(self._alphabet)
        else:
            _reraise_error()
            raise UnexpectedError(status, "p7_oprofile_ReadMSV")

    cpdef void close(self) except *:
        """Close the pressed file and free resources.

        This method has no effect if the file is already closed. It is called
        automatically if the `HMMFile` was used in a context::

            >>> with HMMPressedFile("tests/data/hmms/db/PKSI-AT.hmm") as hmm_db:
            ...     optimized_profile = hmm_db.read()

        """
        if self._hfp:
            libhmmer.p7_hmmfile.p7_hmmfile_Close(self._hfp)
            self._hfp = self._hmmfile._hfp = NULL

cdef class IterationResult:
    """The results of a single iteration from an `IterativeSearch`.

    Attributes:
        hmm (`~pyhmmer.plan7.HMM`): The HMM used to search for hits
            during this iteration.
        hits (`~pyhmmer.plan7.TopHits`): The hits found during this
            iteration.
        msa (`~pyhmmer.easel.DigitalMSA`): A multiple sequence alignment
            containing the hits from this iteration.
        converged (`bool`): A flag marking whether this iteration converged
            (no new hit found in the target sequences with respect to the
            pipeline inclusion thresholds).
        iteration (`int`): The number of iterations done so far. First
            iteration starts at *1*.

    """

    def __cinit__(
        self,
        HMM hmm,
        TopHits hits,
        DigitalMSA msa,
        bint converged,
        size_t iteration
    ):
        self.hmm = hmm
        self.hits = hits
        self.msa = msa
        self.converged = converged
        self.iteration = iteration

    def __iter__(self):
        yield from (self.hmm, self.hits, self.msa, self.converged, self.iteration)


cdef class IterativeSearch:
    """A helper class for running iterative queries like JackHMMER.

    See `Pipeline.iterate_seq` and `Pipeline.iterate_hmm` for more
    information.

    Attributes:
        pipeline (`~pyhmmer.plan7.Pipeline`): The pipeline object to use
            to get hits on each iteration.
        builder (`~pyhmmer.plan7.Builder`): The builder object for
            converting multiple sequence alignments obtained after each
            run to a `~pyhmmer.plan7.HMM`.
        query (`~pyhmmer.easel.DigitalSequence` or `~pyhmmer.plan7.HMM`):
            The query object to use for the first iteration.
        converged (`bool`): Whether the iterative search already converged
            or not.
        targets (`~pyhmmer.easel.DigitalSequenceBlock`): The targets sequences
            to search for homologs.
        ranking (`~pyhmmer.easel.KeyHash`): A mapping storing the rank of
            hits from previous iterations.
        iteration (`int`): The index of the last iteration done so far.

    Yields:
        `~pyhmmer.plan7.IterationResult`: A named tuple containing the hits,
        multiple sequence alignment and HMM for each iteration, as well as
        the iteration index and a flag marking whether the search converged.

    References:
        - Johnson, Steven L., Eddy, Sean R. & Portugaly, Elon.
          *Hidden Markov model speed heuristic and iterative HMM search
          procedure*. BMC Bioinformatics 11, 431 (18 August 2010).
          :doi:`10.1186/1471-2105-11-431`.

    """

    def __init__(
        self,
        Pipeline pipeline,
        Builder builder,
        object query,
        DigitalSequenceBlock targets,
        object select_hits = None,
    ):
        """__init__(self, pipeline, builder, query, targets, select_hits=None)\n--\n
        """
        self.pipeline = pipeline
        self.background = pipeline.background
        self.builder = builder
        self.query = query
        self.converged = False
        self.ranking = KeyHash()
        self.msa = None
        self.iteration = 0
        self.select_hits = select_hits
        self.targets = targets

    def __iter__(self):
        return self

    def __next__(self):
        cdef list    extra_sequences
        cdef list    extra_traces
        cdef HMM     hmm
        cdef TopHits hits
        cdef int     n_new
        cdef int     n_prev

        if self.converged:
            raise StopIteration
        elif self.iteration == 0:
            if isinstance(self.query, HMM):
                hmm = self.query
                n_prev = 1
                extra_sequences = None
                extra_traces = None
            else:
                hmm, _, _ = self.builder.build(self.query, self.background)
                n_prev = 1
                extra_sequences = [self.query]
                extra_traces = [Trace.from_sequence(self.query)]
        else:
            hmm, _, _ = self.builder.build_msa(self.msa, self.background)
            n_prev = len(self.msa.sequences)
            if isinstance(self.query, HMM):
                extra_sequences = None
                extra_traces = None
            else:
                extra_sequences = [self.query]
                extra_traces = [Trace.from_sequence(self.query)]

        hits = self._search_hmm(hmm)
        hits.sort(by="key")
        if self.select_hits is not None:
            self.select_hits(hits)
        n_new = hits.compare_ranking(self.ranking)

        self.msa = hits.to_msa(
            self.pipeline.alphabet,
            sequences=extra_sequences,
            traces=extra_traces,
            all_consensus_cols=True,
            digitize=True,
        )
        self.msa.name = self.query.name + f"-i{self.iteration+1}".encode("utf-8")
        self.msa.description = self.query.description or None
        self.msa.accession = self.query.accession or None
        self.msa.author = b"jackhmmer (pyHMMER)"

        if n_new == 0 and len(self.msa.sequences) <= n_prev:
            self.converged = True

        self.pipeline.clear()
        self.iteration += 1
        return IterationResult(hmm, hits, self.msa, self.converged, self.iteration)

    cpdef TopHits _search_hmm(self, HMM hmm):
        return self.pipeline.search_hmm(hmm, self.targets)


cdef class OptimizedProfile:
    """An optimized profile that uses platform-specific instructions.

    Optimized profiles store the match emissions and transition
    probabilities for a profile HMM so that they can be loaded in the
    SIMD code. Typically, matrices use aligned storage so that they can
    loaded efficiently, and are striped *à la* Farrar to compute
    pairwise scores for each sequence residue and profile node.

    Attributes:
        alphabet (`~pyhmmer.easel.Alphabet`): The alphabet for which this
            optimized profile is configured.

    References:
        - Farrar, Michael.
          *Striped Smith–Waterman Speeds Database Searches Six Times over
          Other SIMD Implementations*. Bioinformatics 23, no. 2
          (15 January 2007): 156–61. :doi:`10.1093/bioinformatics/btl582`.

    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._om = NULL
        self.alphabet = None

    def __init__(self, int M, Alphabet alphabet):
        """__init__(self, M, alphabet)\n--\n

        Create a new optimized profile from scratch.

        Once allocated, you must call the `~OptimizedProfile.convert`
        method with a `~plan7.Profile` object. It's actually easier to
        use `Profile.to_optimized` method to obtained a configured
        `OptimizedProfile` directly, unless you're explicitly trying
        to recycle memory.

        Arguments:
            M (`int`): The length of the model (i.e. the number of nodes).
            alphabet (`~pyhmmer.easel.Alphabet`): The alphabet of the model.

        """
        # store the alphabet to make sure it's not deallocated
        self.alphabet = alphabet
        # create a new optimized profile large enough to store M nodes
        with nogil:
            self._om = p7_oprofile.p7_oprofile_Create(M, alphabet._abc)
        if self._om == NULL:
            raise AllocationError("P7_OPROFILE", sizeof(P7_OPROFILE))

    def __dealloc__(self):
        p7_oprofile.p7_oprofile_Destroy(self._om)

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"<{ty} alphabet={self.alphabet!r} M={self.M!r} name={self.name!r}>"

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        cdef OptimizedProfile new
        if id(self) not in memo:
            new = memo[id(self)] = self.copy()
            new.alphabet = copy.deepcopy(self.alphabet, memo=memo)
            new._om.abc = new.alphabet._abc
        return memo[id(self)]

    def __eq__(self, object other):
        assert self._om != NULL

        if not isinstance(other, OptimizedProfile):
            return NotImplemented

        cdef char[eslERRBUFSIZE] errbuf
        cdef OptimizedProfile    op     = <Profile> other
        cdef int                 status = p7_oprofile_Compare(self._om, op._om, 0.0, errbuf)

        if status == libeasel.eslOK:
            return True
        elif status == libeasel.eslFAIL:
            return False
        else:
            raise UnexpectedError(status, "p7_oprofile_Compare")

    def __sizeof__(self):
        assert self._om != NULL
        return p7_oprofile_Sizeof(self._om) + sizeof(self)

    # --- Properties ---------------------------------------------------------

    @property
    def M(self):
        """`int`: The number of nodes in the model.

        .. versionadded:: 0.4.0

        """
        assert self._om != NULL
        return self._om.M

    @property
    def L(self):
        """`int`: The currently configured target sequence length.

        .. versionadded:: 0.4.0

        """
        assert self._om != NULL
        return self._om.L

    @L.setter
    def L(self, int L):
        assert self._om != NULL
        cdef int status
        with nogil:
            status = p7_oprofile.p7_oprofile_ReconfigLength(self._om, L)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_oprofile_ReconfigLength")

    @property
    def name(self):
        """`bytes` or `None`: The name of the profile, if any.

        .. versionadded:: 0.4.11

        """
        assert self._om != NULL
        return None if self._om.name == NULL else <bytes> self._om.name

    @property
    def accession(self):
        """`bytes` or `None`: The accession of the profile, if any.

        .. versionadded:: 0.4.11

        """
        assert self._om != NULL
        return None if self._om.acc == NULL else <bytes> self._om.acc

    @property
    def description(self):
        """`bytes` or `None`: The description of the profile, if any.

        .. versionadded:: 0.4.11

        """
        assert self._om != NULL
        return None if self._om.desc == NULL else <bytes> self._om.desc

    @property
    def consensus(self):
        """`str` or `None`: The consensus residue line of the profile, if any.

        .. versionadded:: 0.4.11

        """
        assert self._om != NULL
        if self._om.consensus[0] == b'\0':
            return None
        return PyUnicode_FromStringAndSize(&self._om.consensus[1], self._om.M)

    @property
    def consensus_structure(self):
        """`str` or `None`: The consensus structure of the profile, if any.

        .. versionadded:: 0.4.11

        """
        assert self._om != NULL
        if self._om.cs[0] == b'\0':
            return None
        return PyUnicode_FromStringAndSize(&self._om.cs[1], self._om.M)

    @property
    def reference(self):
        """`str` or `None`: The reference line from the alignment, if any.

        This is relevant if the profile was built from a multiple sequence
        alignment (e.g. by `Builder.build_msa`, or by an external
        ``hmmbuild`` pipeline run).

        .. versionadded:: 0.3.1

        """
        assert self._om != NULL
        assert self._om.rf != NULL
        if self._om.rf[0] == b'\0':
            return None
        return PyUnicode_FromStringAndSize(&self._om.rf[1], self._om.M)

    @property
    def model_mask(self):
        """`str` or `None`: The model mask line from the alignment, if any.

        .. versionadded:: 0.3.1

        """
        assert self._om != NULL
        assert self._om.mm != NULL
        if self._om.mm[0] == b'\0':
            return None
        return PyUnicode_FromStringAndSize(&self._om.mm[1], self._om.M)

    # --- MSV Filter ---

    @property
    def rbv(self):
        """`~pyhmmer.easel.MatrixU8`: The match scores for the MSV filter.
        """
        assert self._om != NULL

        cdef MatrixU8 mat = MatrixU8.__new__(MatrixU8)
        mat._m = mat._shape[0] = self.alphabet.Kp
        mat._n = mat._shape[1] = 16 * p7O_NQB(self._om.M)
        mat._owner = self
        mat._data = <void**> self._om.rbv
        return mat

    @property
    def sbv(self):
        """`~pyhmmer.easel.MatrixU8`: The match scores for the SSV filter.

        .. versionadded:: 0.4.0

        """
        assert self._om != NULL

        cdef int nqb = p7O_NQB(self._om.M)
        cdef int nqs = nqb + p7O_EXTRA_SB

        cdef MatrixU8 mat = MatrixU8.__new__(MatrixU8)
        mat._m = mat._shape[0] = self.alphabet.Kp
        mat._n = mat._shape[1] = 16 * nqs
        mat._owner = self
        mat._data = <void**> self._om.sbv
        return mat

    @property
    def tbm(self):
        r"""`int`: The constant cost for a :math:`B \to M_k` transition.

        .. versionadded:: 0.4.0

        """
        assert self._om != NULL
        return self._om.tbm_b

    @property
    def tec(self):
        r"""`int`: The constant cost for a :math:`E \to C` transition.

        .. versionadded:: 0.4.0

        """
        assert self._om != NULL
        return self._om.tec_b

    @property
    def tjb(self):
        """`int`: The constant cost for a :math:`NCJ` move.

        .. versionadded:: 0.4.0

        """
        assert self._om != NULL
        return self._om.tjb_b

    @property
    def scale_b(self):
        """`float`: The scale for MSV filter scores.

        .. versionadded:: 0.10.3

        """
        assert self._om != NULL
        return self._om.scale_b

    @property
    def base_b(self):
        """`int`: The offset for MSV filter scores.

        .. versionadded:: 0.10.3

        """
        assert self._om != NULL
        return self._om.base_b

    @property
    def bias_b(self):
        """`int`: The positive bias to emission scores.

        .. versionadded:: 0.4.0

        .. versionchanged:: 0.10.3
            Renamed from `bias`.

        """
        assert self._om != NULL
        return self._om.bias_b

    # --- ViterbiFilter ---

    # rwv
    # twv
    # xv

    @property
    def scale_w(self):
        """`float`: The scale for Viterbi filter scores.

        .. versionadded:: 0.10.3

        """
        assert self._om != NULL
        return self._om.scale_w

    @property
    def base_w(self):
        """`int`: The offset for Viterbi filter scores.

        .. versionadded:: 0.10.3

        """
        assert self._om != NULL
        return self._om.base_w

    @property
    def ddbound_w(self):
        """`int`: The threshold precalculated for lazy :math:`DD` evaluation.

        .. versionadded:: 0.10.3

        """
        assert self._om != NULL
        return self._om.ddbound_w

    @property
    def ncj_roundoff(self):
        """`float`: The missing precision on :math:`NN,CC,JJ` after rounding.

        .. versionadded:: 0.10.3

        """
        assert self._om != NULL
        return self._om.ncj_roundoff

    # --- Forward, Backard ---

    @property
    def rfv(self):
        """`~pyhmmer.easel.MatrixF`: The match scores for the Forward/Backward.

        .. versionadded:: 0.10.3

        """
        assert self._om != NULL

        cdef MatrixF mat = MatrixF.__new__(MatrixF)
        mat._m = mat._shape[0] = self.alphabet.Kp
        mat._n = mat._shape[1] = 4 * p7O_NQF(self._om.M)
        mat._owner = self
        mat._data = <void**> self._om.rfv
        return mat

    @property
    def tfv(self):
        """`~pyhmmer.easel.VectorF`: The transition scores for the Forward/Backard.

        .. versionadded:: 0.10.3

        """
        assert self._om != NULL

        cdef VectorF vec = VectorF.__new__(VectorF)
        vec._n = vec._shape[0] = 8 * 4 * p7O_NQF(self._om.M)
        vec._owner = self
        vec._data = <void*> self._om.tfv
        return vec

    @property
    def xf(self):
        """`~pyhmmer.easel.MatrixF`: The :math:`NECJ` transition costs.

        .. versionadded:: 0.10.3

        """
        assert self._om != NULL

        cdef MatrixF mat = MatrixF.__new__(MatrixF)
        mat._m = mat._shape[0] = p7O_NXSTATES
        mat._n = mat._shape[1] = p7O_NXTRANS
        mat._owner = self
        mat._data = <void**> self._om.xf
        return mat

    # --- Miscellaneous ---

    @property
    def offsets(self):
        """`~plan7.Offsets`: The disk offsets for this optimized profile.
        """
        assert self._om != NULL
        # expose the disk offsets as the `offsets` attribute
        cdef Offsets offsets = Offsets.__new__(Offsets)
        offsets._offs = &self._om.offs
        offsets._owner = self
        return offsets

    @property
    def evalue_parameters(self):
        """`~plan7.EvalueParameters`: The e-value parameters for this profile.
        """
        assert self._om != NULL
        cdef EvalueParameters ep = EvalueParameters.__new__(EvalueParameters)
        ep._evparams = &self._om.evparam
        ep._owner = self
        return ep

    @property
    def cutoffs(self):
        """`~plan7.Cutoffs`: The bitscore cutoffs for this profile.
        """
        assert self._om != NULL
        cdef Cutoffs cutoffs = Cutoffs.__new__(Cutoffs)
        cutoffs._owner = self
        cutoffs._cutoffs = &self._om.cutoff
        cutoffs._flags = NULL
        cutoffs._is_profile = True
        return cutoffs

    @property
    def compositions(self):
        """`~pyhmmer.easel.VectorF`: The per-model HMM filter composition.

        .. versionadded:: 0.11.0

        """
        assert self._om != NULL
        cdef VectorF vec = VectorF.__new__(VectorF)
        vec._n = vec._shape[0] = p7_MAXABET
        vec._owner = self
        vec._data = <void*> self._om.compo
        return vec

    @property
    def local(self):
        """`bool`: Whether the optimized profile is in local mode.

        .. versionadded:: 0.7.0

        """
        assert self._om != NULL
        return p7_oprofile.p7_oprofile_IsLocal(self._om)

    @property
    def multihit(self):
        """`bool`: Whether the optimized is in multihit mode.

        .. versionadded:: 0.7.0

        """
        assert self._om != NULL
        return self._om.nj == 1.0

    @multihit.setter
    def multihit(self, multihit):
        if multihit:
            if not self.multihit:
                p7_oprofile.p7_oprofile_ReconfigMultihit(self._om, self._om.L)
        else:
            if self.multihit:
                p7_oprofile.p7_oprofile_ReconfigUnihit(self._om, self._om.L)


    # --- Methods ------------------------------------------------------------

    cpdef OptimizedProfile copy(self):
        """Create an exact copy of the optimized profile.

        Note:
            The `Alphabet` referenced to by this object is not copied, use
            `copy.deepcopy` if this is the intended behaviour.

        """
        assert self._om != NULL
        cdef OptimizedProfile new = OptimizedProfile.__new__(OptimizedProfile)
        new.alphabet = self.alphabet
        with nogil:
            new._om = p7_oprofile.p7_oprofile_Copy(self._om)
        if new._om == NULL:
            raise AllocationError("P7_OPROFILE", sizeof(P7_OPROFILE))
        return new

    cpdef void write(self, object fh_filter, object fh_profile) except *:
        """Write an optimized profile to two separate files.

        HMMER implements an acceleration pipeline using several scoring
        algorithms. Parameters for MSV (the *Multi ungapped Segment Viterbi*)
        are saved independently to the ``fh_filter`` handle, while the rest of
        the profile is saved to ``fh_profile``.

        """
        cdef P7_OPROFILE* om     = self._om
        cdef int          status
        cdef FILE*        pfp
        cdef FILE*        ffp

        assert self._om != NULL

        pfp = fopen_obj(fh_profile, "w")
        ffp = fopen_obj(fh_filter, "w")
        status = p7_oprofile_Write(ffp, pfp, self._om)
        if status == libeasel.eslOK:
            fclose(ffp)
            fclose(pfp)
        else:
            raise UnexpectedError(status, "p7_oprofile_Write")

    cpdef void convert(self, Profile profile) except *:
        """Store the given profile into ``self`` as a platform-specific profile.

        Use this method to configure an optimized profile from a `Profile`
        while recycling the internal vector buffers.

        Raises:
            `ValueError`: When the optimized profile is too small to hold the
                profile, or when the standard and the optimized profiles are
                not compatible.

        See Also:
            The `Profile.to_optimized` method, which allows getting an
            `OptimizedProfile` directly from a profile without having to
            allocate first.

        """
        assert self._om != NULL
        assert profile._gm != NULL

        cdef int status

        if self._om.allocM < profile._gm.M:
            raise ValueError("Optimized profile is too small to hold profile")
        with nogil:
            status = p7_oprofile.p7_oprofile_Convert(profile._gm, self._om)
        if status == libeasel.eslEINVAL:
            raise ValueError("Standard and optimized profiles are not compatible.")
        elif status == libeasel.eslEMEM:
            raise AllocationError("P7_OPROFILE", sizeof(P7_OPROFILE))
        elif status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_oprofile_Convert")

    cpdef object ssv_filter(self, DigitalSequence seq):
        """Compute the SSV filter score for the given sequence.

        Arguments:
            seq (`~pyhmmer.easel.DigitalSequence`): The sequence in digital
                format for which to compute the SSV filter score.

        Returns:
            `float` or `None`: The SSV filter score for the sequence.

        Note:
            * `math.inf` may be returned if an overflow occurs that will also
              occur in the MSV filter. This is the case whenever
              :math:`\\text{base} - \\text{tjb} - \\text{tbm} \\ge 128`
            * `None` may be returned if the MSV filter score needs to be
              recomputed (because it may not overflow even though the SSV
              filter did).

        Raises:
            `~pyhmmer.errors.AlphabetMismatch`: When the alphabet of the
                sequence does not correspond to the profile alphabet.

        Caution:
            This method is not available on the PowerPC platform (calling
            it will raise a `NotImplementedError`).

        .. versionadded:: 0.4.0

        """
        assert self._om != NULL

        cdef float score
        cdef int status

        if HMMER_IMPL == "SSE":
            if self.alphabet != seq.alphabet:
                raise AlphabetMismatch(self.alphabet, seq.alphabet)
            with nogil:
                status = p7_SSVFilter(seq._sq.dsq, seq._sq.L, self._om, &score)
            if status == libeasel.eslOK:
                return score
            elif status == libeasel.eslERANGE:
                return math.inf
            elif status == libeasel.eslENORESULT:
                return None
            else:
                raise UnexpectedError(status, "p7_SSVFilter")
        else:
            raise NotImplementedError(f"p7_SSVFilter is not available on {HMMER_IMPL} platforms")

cdef class OptimizedProfileBlock:
    """A container for storing `OptimizedProfile` objects.

    This collections allows the scan loop of `Pipeline` objects to scan
    several target profiles without having to acquire the GIL for each
    new model. It does so by synchronizing a Python `list` storing
    `OptimizedProfile` objects with a C-contiguous array of pointers to the
    underlying struct of each object.

    .. versionadded:: 0.7.0

    """

    def __cinit__(self, Alphabet alphabet, *args, **kwargs):
        self._locks = NULL
        self._block = NULL
        self._storage = list()
        self.alphabet = alphabet

    def __init__(self, Alphabet alphabet, object iterable = ()):
        """__init__(self, alphabet, iterable=())\n--\n

        Create a new block from an iterable of `OptimizedProfile` objects.

        """
        if self._block == NULL:
            self._block = p7_oprofile_CreateBlock(8)
            if self._block == NULL:
                raise AllocationError("P7_OM_BLOCK", sizeof(P7_OM_BLOCK))
        if self._locks == NULL:
            self._locks = <PyThread_type_lock*> calloc(8, sizeof(PyThread_type_lock))
            if self._locks == NULL:
                raise AllocationError("PyThread_type_lock", sizeof(PyThread_type_lock), 8)
            for i in range(self._block.listSize):
                self._locks[i] = PyThread_allocate_lock()
        self.clear()
        self.extend(iterable)

    def __dealloc__(self):
        if self._locks != NULL:
            for i in range(self._block.listSize):
                PyThread_free_lock(self._locks[i])
            free(self._locks)
        if self._block != NULL:
            # avoid a double free of the sequence contents
            for i in range(self._block.listSize):
                self._block.list[i] = NULL
            p7_oprofile_DestroyBlock(self._block)

    def __len__(self):
        assert self._block != NULL
        return self._block.count

    def __contains__(self, object item):
        if not isinstance(item, OptimizedProfile):
            return False
        return item in self._storage

    def __getitem__(self, object index):
        if isinstance(index, slice):
            return type(self)(self.alphabet, self._storage[index])
        else:
            return self._storage[index]

    def __setitem__(self, object index, object optimized_profiles):
        cdef size_t           i
        cdef OptimizedProfile optimized_profile

        if isinstance(index, slice):
            optimized_profiles = list(optimized_profiles)
            for optimized_profile in optimized_profiles:
                if optimized_profile.alphabet != self.alphabet:
                    raise AlphabetMismatch(self.alphabet, optimized_profile.alphabet)
            self._storage[index] = optimized_profiles
            self._block.count = len(self._storage)
            self._allocate(self._block.count)
            for i, optimized_profile in enumerate(self._storage):
                self._block.list[i] = optimized_profile._om
        else:
            optimized_profile = optimized_profiles
            if optimized_profile.alphabet != self.alphabet:
                raise AlphabetMismatch(self.alphabet, optimized_profile.alphabet)
            self._storage[index] = optimized_profile
            self._block.list[index] = optimized_profile._om

    def __delitem__(self, object index):
        cdef size_t           i
        cdef OptimizedProfile optimized_profile

        if isinstance(index, slice):
            del self._storage[index]
            self._block.count = len(self._storage)
            self._allocate(self._block.count)
            for i, optimized_profile in enumerate(self._storage):
                self._block.list[i] = optimized_profile._om
        else:
            self.pop(index)

    def __reduce__(self):
        return type(self), (self.alphabet,), None, iter(self)

    def __copy__(self):
        return self.copy()

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"{ty}({self.alphabet!r}, {self._storage!r})"

    def __eq__(self, object other):
        cdef OptimizedProfileBlock other_
        if not isinstance(other, OptimizedProfileBlock):
            return NotImplemented
        other_ = other
        return self._storage == other_._storage

    def __sizeof__(self):
        return (
                sizeof(self)
            +   sizeof(P7_OM_BLOCK)
            +   self._block.listSize * sizeof(PyThread_type_lock*)
            +   self._block.listSize * sizeof(PyThread_type_lock)
            +   self._block.listSize * sizeof(P7_OPROFILE*)
        )

    # --- C methods ----------------------------------------------------------

    cdef void _allocate(self, size_t n) except *:
        """Allocate enough storage for at least ``n`` items.
        """
        assert self._block != NULL

        cdef ssize_t i
        cdef size_t  capacity = new_capacity(n, self._block.count)

        with nogil:
            self._block.list = <P7_OPROFILE**> realloc(self._block.list, capacity * sizeof(P7_OPROFILE*))
            self._locks = <PyThread_type_lock*> realloc(self._locks, capacity * sizeof(PyThread_type_lock))
        if self._locks == NULL:
            raise AllocationError("PyThread_type_lock", sizeof(PyThread_type_lock), capacity)
        else:
            for i in range(self._block.listSize, capacity):
                self._locks[i] = PyThread_allocate_lock()
        if self._block.list == NULL:
            self._block.listSize = 0
            raise AllocationError("P7_OPROFILE*", sizeof(P7_OPROFILE*), capacity)
        else:
            self._block.listSize = capacity

    # --- Python methods -----------------------------------------------------

    cpdef void append(self, OptimizedProfile optimized_profile) except *:
        """Append an optimized profile at the end of the block.
        """
        assert self._block != NULL
        if self.alphabet != optimized_profile.alphabet:
            raise AlphabetMismatch(self.alphabet, optimized_profile.alphabet)
        if self._block.count == self._block.listSize:
            self._allocate(self._block.count)
        self._storage.append(optimized_profile)
        self._block.list[self._block.count] = optimized_profile._om
        self._block.count += 1

    cpdef void clear(self) except *:
        """Remove all optimized profiles from the block.
        """
        assert self._block != NULL
        cdef ssize_t i
        self._storage.clear()
        self._block.count = 0

    cpdef void extend(self, object iterable) except *:
        """Extend block by appending optimized profiles from the iterable.
        """
        assert self._block != NULL
        cdef ssize_t hint = operator.length_hint(iterable)
        if self._block.count + hint > self._block.listSize:
            self._allocate(self._block.count + hint)
        for optimized_profile in iterable:
            self.append(optimized_profile)

    cpdef size_t index(self, OptimizedProfile optimized_profile, ssize_t start=0, ssize_t stop=sys.maxsize) except *:
        """Return the index of the first occurence of ``optimized_profile``.

        Raises:
            `ValueError`: When the block does not contain ``optimized_profile``.

        """
        return self._storage.index(optimized_profile, start, stop)

    cpdef void insert(self, ssize_t index, OptimizedProfile optimized_profile) except *:
        """Insert a new optimized profile in the block before ``index``.
        """
        assert self._block != NULL

        if self.alphabet != optimized_profile.alphabet:
            raise AlphabetMismatch(self.alphabet, optimized_profile.alphabet)

        if index < 0:
            index = 0
        elif index > self._block.count:
            index = self._block.count

        if self._block.count == self._block.listSize - 1:
            self._allocate(self._block.listSize + 1)

        if index != self._block.count:
            memmove(&self._block.list[index + 1], &self._block.list[index], (self._block.count - index)*sizeof(P7_OPROFILE*))

        self._storage.insert(index, optimized_profile)
        self._block.list[index] = optimized_profile._om
        self._block.count += 1

    cpdef OptimizedProfile pop(self, ssize_t index=-1):
        """Remove and return an optimized profile from the block.
        """
        assert self._block != NULL
        cdef ssize_t index_ = index

        if self._block.count == 0:
            raise IndexError("pop from empty block")
        if index_ < 0:
            index_ += self._block.count
        if index_ < 0 or index_ >= self._block.count:
            raise IndexError(index)

        # remove item from storage
        item = self._storage.pop(index_)

        # update pointers in the reference array
        self._block.count -= 1
        if index_ < self._block.count:
            memmove(&self._block.list[index_], &self._block.list[index_ + 1], (self._block.count - index_)*sizeof(P7_OPROFILE*))
        return item

    cpdef void remove(self, OptimizedProfile optimized_profile) except *:
        """Remove the first occurence of the given optimized profile.
        """
        self.pop(self.index(optimized_profile))

    cpdef OptimizedProfileBlock copy(self):
        """Return a copy of the optimized profile block.

        Note:
            The optimized profiles internally refered to by this collection
            are not copied. Use `copy.deepcopy` if you also want to duplicate
            the internal storage of each optimized profile.

        """
        assert self._block != NULL
        cdef OptimizedProfileBlock new = OptimizedProfileBlock.__new__(OptimizedProfileBlock, self.alphabet)
        new._storage = self._storage.copy()
        new._block = p7_oprofile_CreateBlock(self._block.count)
        if new._block == NULL:
            raise AllocationError("P7_OM_BLOCK", sizeof(P7_OM_BLOCK))
        memcpy(new._block.list, self._block.list, self._block.count * sizeof(ESL_SQ*))
        new._block.count = self._block.count
        new._locks = <PyThread_type_lock*> calloc(self._block.count, sizeof(PyThread_type_lock))
        if new._locks == NULL:
            raise AllocationError("PyThread_type_lock", sizeof(PyThread_type_lock), new._block.count)
        for i in range(new._block.listSize):
            new._locks[i] = PyThread_allocate_lock()
        return new


@cython.freelist(8)
@cython.no_gc_clear
cdef class Offsets:
    """A mutable view over the disk offsets of a profile.
    """

    def __cinit__(self):
        self._owner = None
        self._offs  = NULL

    def __init__(self, object owner):
        """__init__(self, owner)\n--\n
        """
        cdef str ty

        if isinstance(owner, Profile):
            self._offs = &(<Profile> owner)._gm.offs
            self._owner = owner
        elif isinstance(owner, OptimizedProfile):
            self._offs = &(<OptimizedProfile> owner)._om.offs
            self._owner = owner
        else:
            ty = type(owner).__name__
            raise TypeError(f"expected Profile or OptimizedProfile, found {ty}")

    def __copy__(self):
        assert self._offs != NULL
        cdef Offsets copy = Offsets.__new__(Offsets)
        copy._offs = self._offs
        copy._owner = self._owner
        return copy

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"<{ty} model={self.model!r} filter={self.filter!r} profile={self.profile!r}>"

    # --- Properties ---------------------------------------------------------

    @property
    def model(self):
        assert self._offs != NULL
        cdef off_t model = self._offs[0][<int> p7_offsets_e.p7_MOFFSET]
        return None if model == -1 else model

    @model.setter
    def model(self, object model):
        assert self._offs != NULL
        self._offs[0][<int> p7_offsets_e.p7_MOFFSET] = -1 if model is None else model

    @property
    def filter(self):
        assert self._offs != NULL
        cdef off_t filter = self._offs[0][<int> p7_offsets_e.p7_FOFFSET]
        return None if filter == -1 else filter

    @filter.setter
    def filter(self, object filter):
        assert self._offs != NULL
        self._offs[0][<int> p7_offsets_e.p7_FOFFSET] = -1 if filter is None else filter

    @property
    def profile(self):
        assert self._offs != NULL
        cdef off_t profile = self._offs[0][<int> p7_offsets_e.p7_POFFSET]
        return None if profile == -1 else profile

    @profile.setter
    def profile(self, object profile):
        assert self._offs != NULL
        self._offs[0][<int> p7_offsets_e.p7_POFFSET] = -1 if profile is None else profile


cdef double   DEFAULT_F1      = 0.02
cdef double   DEFAULT_F2      = 1e-3
cdef double   DEFAULT_F3      = 1e-5
cdef uint32_t DEFAULT_SEED    = 42
cdef double   DEFAULT_E       = 10.0
cdef double   DEFAULT_DOME    = 10.0
cdef double   DEFAULT_INCE    = 0.01
cdef double   DEFAULT_INCDOME = 0.01
cdef size_t   HMMER_TARGET_LIMIT = 100000

cdef class Pipeline:
    """An HMMER3 accelerated sequence/profile comparison pipeline.

    The Plan7 pipeline handles the platform-accelerated comparison of
    sequences to profile HMMs. It performs either a *search* (comparing
    a single query profile to a target sequence database) or a *scan*
    (comparing a single query sequence to a target profile database). The
    two methods are yielding equivalent results: if you have a collection
    of :math:`M` sequences and :math:`N` HMMs to compare, doing a search
    or a scan should give the same raw scores. The E-values will however
    be different if ``Z`` and  ``domZ`` where not set manually: :math:`Z`
    will be set to :math:`M` for a *search*, and to :math:`N` for a scan.

    The main reason for which you should choose *search* or *scan* is the
    relative size of the sequences and HMMs databases. In the original
    HMMER3 code, the memory was managed in a way that you never had to
    load the entirety of the target sequences in memory. In PyHMMER, the
    methods accept both reading the target database from a file, or loading
    it entirely into memory. A *scan* is always slower than a *search*
    because of the overhead introduced when reconfiguring a profile for a
    new sequence.

    Attributes:
        alphabet (`~pyhmmer.easel.Alphabet`): The alphabet for which the
            pipeline is configured.
        background (`~pyhmmer.plan7.Background`): The null background model
            to use to compute scores.
        randomness (`~pyhmmer.easel.Randomness`): The random number generator
            being used by the pipeline.

    .. versionchanged:: 0.4.2
       Added the ``randomness`` attribute.

    """

    M_HINT = 100         # default model size
    L_HINT = 100         # default sequence size
    _BIT_CUTOFFS = dict(PIPELINE_BIT_CUTOFFS)

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._pli = NULL
        self.alphabet = None
        self.profile = None
        self.background = None
        self.randomness = None

    def __init__(
        self,
        Alphabet alphabet,
        Background background = None,
        *,
        bint bias_filter=True,
        bint null2=True,
        uint32_t seed=DEFAULT_SEED,
        object Z=None,
        object domZ=None,
        double F1=DEFAULT_F1,
        double F2=DEFAULT_F2,
        double F3=DEFAULT_F3,
        double E=DEFAULT_E,
        object T=None,
        double domE=DEFAULT_DOME,
        object domT=None,
        double incE=DEFAULT_INCE,
        object incT=None,
        double incdomE=DEFAULT_INCDOME,
        object incdomT=None,
        str bit_cutoffs=None,
    ):
        """__init__(self, alphabet, background=None, *, bias_filter=True, null2=True, seed=42, Z=None, domZ=None, F1=0.02, F2=1e-3, F3=1e-5, E=10.0, T=None, domE=10.0, domT=None, incE=0.01, incT=None, incdomE=0.01, incdomT=None, bit_cutoffs=None)\n--\n

        Instantiate and configure a new accelerated comparison pipeline.

        Arguments:
            alphabet (`~pyhmmer.easel.Alphabet`): The biological alphabet the
                of the HMMs and sequences that are going to be compared. Used
                to build the background model.
            background (`~pyhmmer.plan7.Background`, optional): The background
                model to use with the pipeline, or ``None`` to create and use
                a default one. *The pipeline needs ownership of the background
                model, so any background model passed there will be copied.*

        Keyword Arguments:
            bias_filter (`bool`): Whether or not to enable composition bias
                filter. Defaults to `True`.
            null2 (`bool`): Whether or not to compute biased composition score
                corrections. Defaults to `True`.
            seed (`int`, optional): The seed to use with the random number
                generator. Pass *0* to use a one-time arbitrary seed, or
                `None` to keep the default seed from HMMER.
            Z (`int`, optional): The effective number of comparisons done,
                for E-value calculation. Leave as `None` to auto-detect by
                counting the number of sequences queried.
            domZ (`int`, optional): The number of significant sequences found,
                for domain E-value calculation. Leave as `None` to auto-detect
                by counting the number of sequences reported.
            F1 (`float`): The MSV filter threshold.
            F2 (`float`): The Viterbi filter threshold.
            F3 (`float`): The uncorrected Forward filter threshold.
            E (`float`): The per-target E-value threshold for reporting
                a hit.
            T (`float`, optional): The per-target bit score threshold for
                reporting a hit. *If given, takes precedence over* ``E``.
            domE (`float`): The per-domain E-value threshold for reporting
                a domain hit.
            domT (`float`, optional): The per-domain bit score threshold for
                reporting a domain hit. *If given, takes precedence over*
                ``domE``.
            incE (`float`): The per-target E-value threshold for including
                a hit in the resulting `TopHits`.
            incT (`float`, optional): The per-target bit score threshold for
                including a hit in the resulting `TopHits`. *If given, takes
                precedence over* ``incE``.
            incdomE (`float`): The per-domain E-value threshold for including
                a domain in the resulting `TopHits`.
            incdomT (`float`, optional): The per-domain bit score thresholds
                for including a domain in the resulting `TopHits`. *If given,
                takes precedence over* ``incdomE``.
            bit_cutoffs (`str`, optional): The model-specific thresholding
                option to use for reporting hits. With `None` (the default),
                use global pipeline options; otherwise pass one of
                ``"noise"``, ``"gathering"`` or ``"trusted"`` to use the
                appropriate cutoffs.

        Hint:
            In order to run the pipeline in slow/max mode, disable the bias
            filter, and set the three filtering parameters to :math:`1.0`::

                >>> dna = easel.Alphabet.dna()
                >>> max_pli = Pipeline(dna, bias_filter=False, F1=1.0, F2=1.0, F3=1.0)

        .. versionchanged:: 0.4.6
           Added keywords arguments to set the E-value thresholds.

        """
        # extract default parameters from the class attributes
        cdef int m_hint = self.M_HINT
        cdef int l_hint = self.L_HINT

        # store a reference to the alphabet to avoid deallocation
        self.alphabet = alphabet

        # use the background model or create a default one
        if background is None:
            self.background = Background(alphabet)
        else:
            self.background = background.copy()

        # allocate the pipeline
        with nogil:
            self._pli = libhmmer.p7_pipeline.p7_pipeline_Create(
                NULL,
                m_hint,
                l_hint,
                False, # long targets
                p7_pipemodes_e.p7_SEARCH_SEQS
            )
        if self._pli == NULL:
            raise AllocationError("P7_PIPELINE", sizeof(P7_PIPELINE))

        # create a Randomness object exposing the internal pipeline RNG
        self.randomness = Randomness.__new__(Randomness)
        self.randomness._owner = self
        self.randomness._rng = self._pli.r

        # create an empty profile and optimized profile
        self.profile = Profile(m_hint, self.alphabet)
        self.opt = OptimizedProfile(m_hint, self.alphabet)

        # configure the pipeline with the additional keyword arguments
        self.null2 = null2
        self.bias_filter = bias_filter
        self.Z = Z
        self.domZ = domZ
        self.seed = seed
        self.F1 = F1
        self.F2 = F2
        self.F3 = F3
        self.E = E
        self.T = T
        self.domE = domE
        self.domT = domT
        self.incE = incE
        self.incT = incT
        self.incdomE = incdomE
        self.incdomT = incdomT

        # setup the model-specific reporting cutoffs
        self._save_cutoff_parameters()
        self.bit_cutoffs = bit_cutoffs

    def __dealloc__(self):
        libhmmer.p7_pipeline.p7_pipeline_Destroy(self._pli)

    # --- Properties ---------------------------------------------------------

    @property
    def Z(self):
        """`float` or `None`: The number of effective targets searched.

        It is used to compute the independent e-value for each domain, and
        for an entire hit. If `None`, the parameter number will be set
        automatically after all the comparisons have been done. Otherwise,
        it can be set to an arbitrary number.

        """
        return self._Z

    @Z.setter
    def Z(self, object Z):
        assert self._pli != NULL
        if Z is None:
            self._pli.Z       = 0.0
            self._pli.Z_setby = p7_zsetby_e.p7_ZSETBY_NTARGETS
            self._Z           = None
        else:
            self._pli.Z_setby = p7_zsetby_e.p7_ZSETBY_OPTION
            self._pli.Z = self._Z = Z

    @property
    def domZ(self):
        """`float` or `None`: The number of significant targets.

        It is used to compute the conditional e-value for each domain. If
        `None`, the parameter number will be set automatically after all
        the comparisons have been done, and all hits have been thresholded.
        Otherwise, it can be set to an arbitrary number.

        """
        return self._domZ

    @domZ.setter
    def domZ(self, object domZ):
        assert self._pli != NULL
        if domZ is None:
            self._pli.domZ       = 0.0
            self._pli.domZ_setby = p7_zsetby_e.p7_ZSETBY_NTARGETS
            self._domZ           = None
        else:
            self._pli.domZ_setby = p7_zsetby_e.p7_ZSETBY_OPTION
            self._pli.domZ = self._domZ = domZ

    @property
    def seed(self):
        """`int`: The seed given at pipeline initialization.

        Setting this attribute to a different value will cause the random
        number generator to be reseeded immediately.

        .. versionadded:: 0.2.0

        .. versionchanged:: 0.4.2
           Avoid shadowing initial null seeds given on pipeline initialization.

        """
        return self._seed

    @seed.setter
    def seed(self, uint32_t seed):
        self._seed = seed
        self._pli.do_reseeding = self._pli.ddef.do_reseeding = seed != 0
        self.randomness.seed(seed)

    @property
    def null2(self):
        """`bool`: Whether or not to enable the *null2* score correction.

        .. versionadded:: 0.4.1

        """
        assert self._pli != NULL
        return self._pli.do_null2

    @null2.setter
    def null2(self, bint null2):
        assert self._pli != NULL
        self._pli.do_null2 = null2

    @property
    def bias_filter(self):
        """`bool`: Whether or not to enable the biased comp HMM filter.

        .. versionadded:: 0.4.1

        """
        assert self._pli != NULL
        return self._pli.do_biasfilter

    @bias_filter.setter
    def bias_filter(self, bint bias_filter):
        assert self._pli != NULL
        self._pli.do_biasfilter = bias_filter

    @property
    def F1(self):
        """`float`: The MSV filter threshold.

        .. versionadded:: 0.4.1

        """
        assert self._pli != NULL
        return self._pli.F1

    @F1.setter
    def F1(self, double F1):
        assert self._pli != NULL
        self._pli.F1 = F1

    @property
    def F2(self):
        """`float`: The Viterbi filter threshold.

        .. versionadded:: 0.4.1

        """
        assert self._pli != NULL
        return self._pli.F2

    @F2.setter
    def F2(self, double F2):
        assert self._pli != NULL
        self._pli.F2 = F2

    @property
    def F3(self):
        """`float`: The uncorrected Forward filter threshold.

        .. versionadded:: 0.4.1

        """
        assert self._pli != NULL
        return self._pli.F3

    @F3.setter
    def F3(self, double F3):
        assert self._pli != NULL
        self._pli.F3 = F3

    @property
    def E(self):
        """`float`: The per-target E-value threshold for reporting a hit.

        .. versionadded:: 0.4.6

        """
        assert self._pli != NULL
        return self._pli.E

    @E.setter
    def E(self, double E):
        assert self._pli != NULL
        self._pli.E = E

    @property
    def T(self):
        """`float` or `None`: The per-target score threshold for reporting a hit.

        If set to a non-`None` value, this threshold takes precedence over
        the per-target E-value threshold (`Pipeline.E`).

        .. versionadded:: 0.4.8

        """
        assert self._pli != NULL
        return None if self._pli.by_E else self._pli.T

    @T.setter
    def T(self, object T):
        assert self._pli != NULL
        if T is None:
            self._pli.T = 0.0
            self._pli.by_E = True
        else:
            self._pli.T = T
            self._pli.by_E = False

    @property
    def domE(self):
        """`float`: The per-domain E-value threshold for reporting a hit.

        .. versionadded:: 0.4.6

        """
        assert self._pli != NULL
        return self._pli.domE

    @domE.setter
    def domE(self, double domE):
        assert self._pli != NULL
        self._pli.domE = domE

    @property
    def domT(self):
        """`float` or `None`: The per-domain score threshold for reporting a hit.

        If set to a non-`None` value, this threshold takes precedence over
        the per-domain E-value threshold (`Pipeline.domE`).

        .. versionadded:: 0.4.8

        """
        assert self._pli != NULL
        return None if self._pli.dom_by_E else self._pli.domT

    @domT.setter
    def domT(self, object domT):
        assert self._pli != NULL
        if domT is None:
            self._pli.domT = 0.0
            self._pli.dom_by_E = True
        else:
            self._pli.domT = domT
            self._pli.dom_by_E = False

    @property
    def incE(self):
        """`float`: The per-target E-value threshold for including a hit.

        .. versionadded:: 0.4.6

        """
        assert self._pli != NULL
        return self._pli.incE

    @incE.setter
    def incE(self, double incE):
        assert self._pli != NULL
        self._pli.incE = incE

    @property
    def incT(self):
        """`float` or `None`: The per-target score threshold for including a hit.

        If set to a non-`None` value, this threshold takes precedence over
        the per-target E-value inclusion threshold (`Pipeline.incE`).

        .. versionadded:: 0.4.8

        """
        assert self._pli != NULL
        return None if self._pli.inc_by_E else self._pli.incT

    @incT.setter
    def incT(self, object incT):
        assert self._pli != NULL
        if incT is None:
            self._pli.incT = 0.0
            self._pli.inc_by_E = True
        else:
            self._pli.incT = incT
            self._pli.inc_by_E = False

    @property
    def incdomE(self):
        """`float`: The per-domain E-value threshold for including a hit.

        .. versionadded:: 0.4.6

        """
        assert self._pli != NULL
        return self._pli.incdomE

    @incdomE.setter
    def incdomE(self, double incdomE):
        assert self._pli != NULL
        self._pli.incdomE = incdomE

    @property
    def incdomT(self):
        """`float` or `None`: The per-domain score threshold for including a hit.

        If set to a non-`None` value, this threshold takes precedence over
        the per-domain E-value inclusion threshold (`Pipeline.incdomE`).

        .. versionadded:: 0.4.8

        """
        assert self._pli != NULL
        return None if self._pli.incdom_by_E else self._pli.incdomT

    @incdomT.setter
    def incdomT(self, object incdomT):
        assert self._pli != NULL
        if incdomT is None:
            self._pli.incdomT = 0.0
            self._pli.incdom_by_E = True
        else:
            self._pli.incdomT = incdomT
            self._pli.incdom_by_E = False

    @property
    def bit_cutoffs(self):
        """`str` or `None`: The model-specific thresholding option, if any.

        .. versionadded:: 0.4.6

        """
        assert self._pli != NULL
        return next(
            (k for k,v in PIPELINE_BIT_CUTOFFS.items() if v == self._pli.use_bit_cutoffs),
            None
        )

    @bit_cutoffs.setter
    def bit_cutoffs(self, str bit_cutoffs):
        assert self._pli != NULL
        if bit_cutoffs is not None:
            #
            flag = PIPELINE_BIT_CUTOFFS.get(bit_cutoffs)
            if flag is None:
                raise InvalidParameter("bit_cutoffs", bit_cutoffs, choices=list(PIPELINE_BIT_CUTOFFS) + [None])
            self._pli.use_bit_cutoffs = flag
            # save previous values before overwriting
            self._save_cutoff_parameters()
            # disable thresholding by E and T values
            self._pli.T = self._pli.domT = 0.0
            self._pli.by_E = self._pli.dom_by_E = False
            self._pli.incT = self._pli.incdomT = 0.0
            self._pli.inc_by_E = self._pli.incdom_by_E = False
        else:
            # disable bit cutoffs
            self._pli.use_bit_cutoffs = False
            # restore previous configuration
            self._restore_cutoff_parameters()

    # --- Utils --------------------------------------------------------------

    cdef int _save_cutoff_parameters(self) except 1:
        assert self._pli != NULL
        self._cutoff_save = {
            "T": self._pli.T,
            "domT": self._pli.domT,
            "by_E": self._pli.by_E,
            "dom_by_E": self._pli.dom_by_E,
            "incT": self._pli.incT,
            "incdomT": self._pli.incdomT,
            "inc_by_E": self._pli.inc_by_E,
            "incdom_by_E": self._pli.incdom_by_E,
        }
        return 0

    cdef int _restore_cutoff_parameters(self) except 1:
        assert self._pli != NULL
        self._pli.T = self._cutoff_save['T']
        self._pli.domT = self._cutoff_save['domT']
        self._pli.by_E = self._cutoff_save['by_E']
        self._pli.dom_by_E = self._cutoff_save['dom_by_E']
        self._pli.incT = self._cutoff_save['incT']
        self._pli.incdomT = self._cutoff_save['incdomT']
        self._pli.inc_by_E = self._cutoff_save['inc_by_E']
        self._pli.incdom_by_E = self._cutoff_save['incdom_by_E']

    cdef P7_OPROFILE* _get_om_from_query(self, object query, int L = L_HINT) except NULL:
        assert self._pli != NULL

        if isinstance(query, OptimizedProfile):
            return (<OptimizedProfile> query)._om

        if isinstance(query, HMM):
            # reallocate the profile if it is too small, otherwise just clear it
            if self.profile._gm.allocM < query.M:
                libhmmer.p7_profile.p7_profile_Destroy(self.profile._gm)
                self.profile._gm = libhmmer.p7_profile.p7_profile_Create(query.M, self.alphabet._abc)
                if self.profile._gm == NULL:
                    raise AllocationError("P7_PROFILE", sizeof(P7_OPROFILE))
            else:
                self.profile.clear()
            # configure the profile from the query HMM
            self.profile.configure(<HMM> query, self.background, L)
            # use the local profile as a query
            query = self.profile

        if isinstance(query, Profile):
            # reallocate the optimized profile if it is too small
            if self.opt._om.allocM < query.M:
                p7_oprofile.p7_oprofile_Destroy(self.opt._om)
                self.opt._om = p7_oprofile.p7_oprofile_Create(query.M, self.alphabet._abc)
                if self.opt._om == NULL:
                    raise AllocationError("P7_OPROFILE", sizeof(P7_OPROFILE))
            # convert the profile to an optimized one
            self.opt.convert(self.profile)
            # use the temporary optimized profile
            return self.opt._om

        else:
            ty = type(query).__name__
            raise TypeError(f"Expected HMM, Profile or OptimizedProfile, found {ty}")

    @staticmethod
    cdef int _missing_cutoffs(const P7_PIPELINE* pli, const P7_OPROFILE* om) except 1 nogil:
        with gil:
            bit_cutoffs = next(
                (k for k,v in PIPELINE_BIT_CUTOFFS.items() if v == pli.use_bit_cutoffs),
                None
            )
            model_name = PyUnicode_DecodeASCII(om.name, strlen(om.name), "replace")
            raise MissingCutoffs(model_name, bit_cutoffs)

    # --- Methods ------------------------------------------------------------

    cpdef list arguments(self):
        """Generate an ``argv``-like list with pipeline options as CLI flags.

        Example:
            >>> alphabet = easel.Alphabet.amino()
            >>> plan7.Pipeline(alphabet).arguments()
            []
            >>> plan7.Pipeline(alphabet, F1=0.01).arguments()
            ['--F1', '0.01']

        .. versionadded:: 0.6.0

        """
        cdef list argv = []
        cdef str  cutoff

        if self._pli.use_bit_cutoffs:
            cutoff = self.bit_cutoffs
            if cutoff == "gathering":
                argv.append("--cut_ga")
            elif cutoff == "noise":
                argv.append("--cut_nc")
            elif cutoff == "trusted":
                argv.append("--cut_tc")
            else:
                raise ValueError(f"Unknown bit cutoffs: {cutoff!r}")
        else:
            # options controlling reporting thresholds
            if self._pli.by_E:
                if self._pli.E != DEFAULT_E:
                    argv.append("-E")
                    argv.append(str(self._pli.E))
            else:
                argv.append("-T")
                argv.append(str(self._pli.T))
            if self._pli.dom_by_E:
                if self._pli.domE != DEFAULT_DOME:
                    argv.append("--domE")
                    argv.append(str(self._pli.domE))
            else:
                argv.append("--domT")
                argv.append(str(self._pli.domT))
            # options controlling inclusion thresholds
            if self._pli.inc_by_E:
                if self._pli.incE != DEFAULT_INCE:
                    argv.append("--incE")
                    argv.append(str(self._pli.incE))
            else:
                argv.append("--incT")
                argv.append(str(self._pli.incT))
            if self._pli.incdom_by_E:
                if self._pli.incdomE != DEFAULT_INCDOME:
                    argv.append("--incdomE")
                    argv.append(str(self._pli.incdomE))
            else:
                argv.append("--incdomT")
                argv.append(str(self._pli.incdomT))

        if self._pli.F1 != DEFAULT_F1:
            argv.append("--F1")
            argv.append(str(self.F1))
        if self._pli.F2 != DEFAULT_F2:
            argv.append("--F2")
            argv.append(str(self.F2))
        if self._pli.F3 != DEFAULT_F3:
            argv.append("--F3")
            argv.append(str(self.F3))

        if not self.bias_filter:
            argv.append("--nobias")
        if not self.null2:
            argv.append("--nonull2")

        if self._pli.Z_setby == p7_zsetby_e.p7_ZSETBY_OPTION:
            argv.append("-Z")
            argv.append(str(self._pli.Z))
        if self._pli.domZ_setby == p7_zsetby_e.p7_ZSETBY_OPTION:
            argv.append("--domZ")
            argv.append(str(self._pli.domZ))

        if self._seed != DEFAULT_SEED:
            argv.append("--seed")
            argv.append(str(self._seed))

        return argv

    cpdef void clear(self):
        """Reset the pipeline configuration to its default state.
        """
        assert self._pli != NULL

        cdef int      status
        cdef uint32_t seed

        # reset the Z and domZ values from the CLI if needed
        if self._pli.Z_setby == p7_zsetby_e.p7_ZSETBY_NTARGETS:
            self.Z = None
        if self._pli.domZ_setby == p7_zsetby_e.p7_ZSETBY_NTARGETS:
            self.domZ = None

        # reinitialize the random number generator, even if
        # `self._pli.do_reseeding` is False, because a true
        # deallocation/reallocation of a P7_PIPELINE would reinitialize
        # it unconditionally.
        self.randomness.seed(self._seed)
        # reinitialize the domaindef
        libhmmer.p7_domaindef.p7_domaindef_Reuse(self._pli.ddef)

        # Reset accounting values
        self._pli.nmodels         = 0
        self._pli.nseqs           = 0
        self._pli.nres            = 0
        self._pli.nnodes          = 0
        self._pli.n_past_msv      = 0
        self._pli.n_past_bias     = 0
        self._pli.n_past_vit      = 0
        self._pli.n_past_fwd      = 0
        self._pli.n_output        = 0
        self._pli.pos_past_msv    = 0
        self._pli.pos_past_bias   = 0
        self._pli.pos_past_vit    = 0
        self._pli.pos_past_fwd    = 0
        self._pli.pos_output      = 0
        self._pli.W               = 0
        self._pli.hfp             = NULL
        self._pli.errbuf[0]       = b'\0'

    cpdef TopHits search_hmm(
        self,
        object query,
        SearchTargets sequences
    ):
        """Run the pipeline using a query HMM against a sequence database.

        Arguments:
            query (`HMM`, `Profile` or `OptimizedProfile`): The object to use
                to query the sequence database.
            sequences (`DigitalSequenceBlock` or `SequenceFile`): The target
                sequences to query with the HMM, either pre-loaded in memory
                inside a `~pyhmmer.easel.DigitalSequenceBlock`, or to be
                read iteratively from a `SequenceFile` opened in digital mode.

        Returns:
            `~pyhmmer.plan7.TopHits`: the hits found in the sequence database.

        Raises:
            `ValueError`: When the pipeline is configured to use model-specific
                reporting thresholds but the `HMM` query doesn't have the right
                cutoffs available.
            `~pyhmmer.errors.AlphabetMismatch`: When the alphabet of the
                current pipeline does not match the alphabet of the given
                HMM.

        Hint:
            This method corresponds to running ``hmmsearch`` with the
            ``query`` HMM against the ``sequences`` database.

        .. versionadded:: 0.2.0

        .. versionchanged:: 0.4.9
            Query can now be a `Profile` or an `OptimizedProfile`.

        .. versionchanged:: 0.7.0
            Targets can be inside a `DigitalSequenceBlock` or a `SequenceFile`.

        """
        assert self._pli != NULL

        cdef size_t       L
        cdef str          ty
        cdef P7_OPROFILE* om
        cdef int          status
        cdef int          allocM
        cdef TopHits      hits   = TopHits(query)

        # check that the sequence file is in digital mode
        if SearchTargets is SequenceFile:
            if sequences.name is None:
                raise ValueError("can only use a `SequenceFile` backed by a file for reading targets")
            if not sequences.digital:
                raise ValueError("target sequences file is not in digital mode")
        # check that all alphabets are consistent
        if not self.alphabet._eq(query.alphabet):
            raise AlphabetMismatch(self.alphabet, query.alphabet)
        if not self.alphabet._eq(sequences.alphabet):
            raise AlphabetMismatch(self.alphabet, sequences.alphabet)
        # get a length estimate for profile configuration
        if SearchTargets is DigitalSequenceBlock:
            L = self.L_HINT if not sequences else sequences._refs[0].L
            if sequences and len(sequences.largest()) > HMMER_TARGET_LIMIT:
               raise ValueError(f"sequence length over comparison pipeline limit ({HMMER_TARGET_LIMIT})")
        else:
            L = self.L_HINT
        # get the optimized profile from the query
        om = self._get_om_from_query(query, L=L)

        with nogil:
            # make sure the pipeline is set to search mode and ready for a new HMM
            self._pli.mode = p7_pipemodes_e.p7_SEARCH_SEQS
            self._pli.nseqs = 0
            # run the search loop on all database sequences while recycling memory
            if SearchTargets is DigitalSequenceBlock:
                Pipeline._search_loop(
                    self._pli,
                    om,
                    self.background._bg,
                    <const ESL_SQ**> sequences._refs,
                    sequences._length,
                    hits._th,
                )
            elif SearchTargets is SequenceFile:
                Pipeline._search_loop_file(
                    self._pli,
                    om,
                    self.background._bg,
                    sequences._sqfp,
                    hits._th,
                )
            else:
                raise NotImplementedError("Pipeline.search_hmm")
            # sort hits and set bookkeeping attributes
            hits._sort_by_key()
            hits._threshold(self)

        # record the query metadata
        hits._query = query
        # return the hits
        return hits

    cpdef TopHits search_msa(
        self,
        DigitalMSA query,
        object sequences,
        Builder builder = None,
    ):
        """Run the pipeline using a query alignment against a sequence database.

        Arguments:
            query (`~pyhmmer.easel.DigitalMSA`): The multiple sequence
                alignment to use to query the sequence database.
            sequences (`DigitalSequenceBlock` or `SequenceFile`): The target
                sequences to query with the alignment, either pre-loaded in
                memory inside a `pyhmmer.easel.DigitalSequenceBlock`, or to be
                read iteratively from a `SequenceFile` opened in digital mode.
            builder (`~pyhmmer.plan7.Builder`, optional): A HMM builder to
                use to convert the query to a `~pyhmmer.plan7.HMM`. If
                `None` is given, it will use a default one.

        Returns:
            `~pyhmmer.plan7.TopHits`: the hits found in the sequence database.

        Raises:
            `~pyhmmer.errors.AlphabetMismatch`: When the alphabet of the
                current pipeline does not match the alphabet of the given
                query.
            `ValueError`: When the `Builder` fails to create an HMM from
                the given `DigitalMSA` query.

        Hint:
            This method corresponds to running ``phmmer`` with the
            ``query`` multiple sequence alignment against the ``sequences``
            database.

        Caution:
            Internally, this method will create a new HMM from the query MSA
            using the `Builder.build_msa` method. HMMER requires that every
            HMM has a name, so the `Builder` will attempt to use the name
            of the query MSA to name the HMM. Passing an MSA without a name
            will result in an error.

        .. versionadded:: 0.3.0

        .. versionchanged:: 0.7.0
            Targets can be inside a `DigitalSequenceBlock` or a `SequenceFile`.

        """
        assert self._pli != NULL

        cdef HMM              hmm
        cdef OptimizedProfile opt
        cdef Profile          profile
        cdef TopHits          hits

        # check the pipeline was configured with the same alphabet
        if not self.alphabet._eq(query.alphabet):
            raise AlphabetMismatch(self.alphabet, query.alphabet)
        # use a default HMM builder if none was given
        builder = Builder(self.alphabet, seed=self.seed) if builder is None else builder
        # build the HMM and the profile from the query MSA
        hmm, profile, opt = builder.build_msa(query, self.background)
        if isinstance(sequences, DigitalSequenceBlock):
            hits = self.search_hmm[DigitalSequenceBlock](opt, sequences)
        elif isinstance(sequences, SequenceFile):
            hits = self.search_hmm[SequenceFile](opt, sequences)
        else:
            ty = type(sequences).__name__
            raise TypeError(f"Expected DigitalSequenceBlock or SequenceFile, found {ty}")
        # record query metadata
        hits._query = query
        return hits

    cpdef TopHits search_seq(
        self,
        DigitalSequence query,
        object sequences,
        Builder builder = None,
    ):
        """Run the pipeline using a query sequence against a sequence database.

        Arguments:
            query (`~pyhmmer.easel.DigitalSequence`): The sequence object to
                use to query the sequence database.
            sequences (`DigitalSequenceBlock` or `SequenceFile`): The target
                sequences to query with the query sequence, either pre-loaded
                in memory inside a `pyhmmer.easel.DigitalSequenceBlock`, or to
                be read iteratively from a `SequenceFile` opened in digital
                mode.
            builder (`~pyhmmer.plan7.Builder`, optional): A HMM builder to
                use to convert the query to a `~pyhmmer.plan7.HMM`. If
                `None` is given, it will use a default one.

        Returns:
            `~pyhmmer.plan7.TopHits`: the hits found in the sequence database.

        Raises:
            `~pyhmmer.errors.AlphabetMismatch`: When the alphabet of the
                current pipeline does not match the alphabet of the given
                query.
            `ValueError`: When the `Builder` fails to create an HMM from
                the given `DigitalSequence` query.

        Hint:
            This method corresponds to running ``phmmer`` with the
            ``query`` sequence against the ``sequences`` database.

        .. versionadded:: 0.2.0

        .. versionchanged:: 0.7.0
            Targets can be inside a `DigitalSequenceBlock` or a `SequenceFile`.

        """
        assert self._pli != NULL

        cdef HMM              hmm
        cdef OptimizedProfile opt
        cdef Profile          profile
        cdef TopHits          hits

        # check the pipeline was configure with the same alphabet
        if not self.alphabet._eq(query.alphabet):
            raise AlphabetMismatch(self.alphabet, query.alphabet)
        # make sure the pipeline is set to search mode and ready for a new HMM
        self._pli.mode = p7_pipemodes_e.p7_SEARCH_SEQS
        # use a default HMM builder if none was given
        builder = Builder(self.alphabet, seed=self.seed) if builder is None else builder
        # build the HMM and the profile from the query sequence
        hmm, profile, opt = builder.build(query, self.background)
        if isinstance(sequences, DigitalSequenceBlock):
            hits = self.search_hmm[DigitalSequenceBlock](opt, sequences)
        elif isinstance(sequences, SequenceFile):
            hits = self.search_hmm[SequenceFile](opt, sequences)
        else:
            ty = type(sequences).__name__
            raise TypeError(f"Expected DigitalSequenceBlock or SequenceFile, found {ty}")
        # record query metadata
        hits._query = query
        return hits

    @staticmethod
    cdef int _search_loop(
              P7_PIPELINE* pli,
              P7_OPROFILE* om,
              P7_BG*       bg,
        const ESL_SQ**     sq,
        const size_t       n_targets,
              P7_TOPHITS*  th,
    ) except 1 nogil:
        """Run the low-level search loop while the GIL is released.

        Arguments:
            pli (``P7_PIPELINE*``): A raw pointer to the pipeline underlying
                the `Pipeline` object running the search.
            om (``P7_OPROFILE*``): A raw pointer to the optimized profile
                being used to query the pipeline.
            bg (``P7_BG*``): A raw pointer to the background model used to
                compute the p-value and E-value for each alignment.
            sq (``ESL_SQ**``): An array of pointers to the target sequences
                being queried.
            n_targets (``size_t``): The number of target sequences inside the
                ``sq`` array.
            th (``P7_TOPHITS*``): A raw pointer to the hits underlying the
                `TopHits` object storing the results.

        """
        cdef int    status
        cdef size_t t

        # configure the pipeline for the current HMM
        status = libhmmer.p7_pipeline.p7_pli_NewModel(pli, om, bg)
        if status == libeasel.eslEINVAL:
            Pipeline._missing_cutoffs(pli, om)
        elif status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_pli_NewModel")

        # run the inner loop on all sequences
        for t in range(n_targets):
            # configure the profile, background and pipeline for the new sequence
            status = libhmmer.p7_pipeline.p7_pli_NewSeq(pli, sq[t])
            if status != libeasel.eslOK:
                raise UnexpectedError(status, "p7_pli_NewSeq")
            status = libhmmer.p7_bg.p7_bg_SetLength(bg, sq[t].n)
            if status != libeasel.eslOK:
                raise UnexpectedError(status, "p7_bg_SetLength")
            status = p7_oprofile.p7_oprofile_ReconfigLength(om, sq[t].n)
            if status != libeasel.eslOK:
                raise UnexpectedError(status, "p7_oprofile_ReconfigLength")
            # run the pipeline on the target sequence
            status = libhmmer.p7_pipeline.p7_Pipeline(pli, om, bg, sq[t], NULL, th)
            if status == libeasel.eslEINVAL:
                Pipeline._missing_cutoffs(pli, om)
            elif status == libeasel.eslERANGE:
                raise OverflowError("numerical overflow in the optimized vector implementation")
            elif status != libeasel.eslOK:
                raise UnexpectedError(status, "p7_Pipeline")
            # clear pipeline for reuse for next target
            libhmmer.p7_pipeline.p7_pipeline_Reuse(pli)

        # Return 0 to indicate success
        return 0

    @staticmethod
    cdef int _search_loop_file(
              P7_PIPELINE* pli,
              P7_OPROFILE* om,
              P7_BG*       bg,
              ESL_SQFILE*  sqfp,
              P7_TOPHITS*  th,
    ) except 1 nogil:
        """Run the low-level search loop while the GIL is released.

        Arguments:
            pli (``P7_PIPELINE*``): A raw pointer to the pipeline underlying
                the `Pipeline` object running the search.
            om (``P7_OPROFILE*``): A raw pointer to the optimized profile
                being used to query the pipeline.
            bg (``P7_BG*``): A raw pointer to the background model used to
                compute the p-value and E-value for each alignment.
            sq (``ESL_SQFILE``): A sequence file from which to read the
                target sequences being queried.
            th (``P7_TOPHITS*``): A raw pointer to the hits underlying the
                `TopHits` object storing the results.

        """
        cdef int     status
        cdef ESL_SQ* dbsq

        # allocate a temporary sequence to read targets into
        dbsq = libeasel.sq.esl_sq_CreateDigital(om.abc)
        if dbsq == NULL:
            raise AllocationError("ESL_SQ", sizeof(ESL_SQ))

        try:
            # configure the pipeline for the current HMM
            status = libhmmer.p7_pipeline.p7_pli_NewModel(pli, om, bg)
            if status == libeasel.eslEINVAL:
                Pipeline._missing_cutoffs(pli, om)
            elif status != libeasel.eslOK:
                raise UnexpectedError(status, "p7_pli_NewModel")
            # run the inner loop on all sequences
            while True:
                # read the next sequence
                status = libeasel.sqio.esl_sqio_Read(sqfp, dbsq)
                if status == libeasel.eslEOF:
                    break
                elif status == libeasel.eslEFORMAT:
                    raise ValueError("Could not parse file")
                elif status != libeasel.eslOK:
                    raise UnexpectedError(status, "p7_oprofile_ReadMSV")
                # check sequence length
                if dbsq.L > HMMER_TARGET_LIMIT:
                    raise ValueError(f"sequence length over comparison pipeline limit ({HMMER_TARGET_LIMIT})")
                # configure the profile, background and pipeline for the new sequence
                status = libhmmer.p7_pipeline.p7_pli_NewSeq(pli, dbsq)
                if status != libeasel.eslOK:
                    raise UnexpectedError(status, "p7_pli_NewSeq")
                status = libhmmer.p7_bg.p7_bg_SetLength(bg, dbsq.n)
                if status != libeasel.eslOK:
                    raise UnexpectedError(status, "p7_bg_SetLength")
                status = p7_oprofile.p7_oprofile_ReconfigLength(om, dbsq.n)
                if status != libeasel.eslOK:
                    raise UnexpectedError(status, "p7_oprofile_ReconfigLength")
                # run the pipeline on the target sequence
                status = libhmmer.p7_pipeline.p7_Pipeline(pli, om, bg, dbsq, NULL, th)
                if status == libeasel.eslEINVAL:
                    Pipeline._missing_cutoffs(pli, om)
                elif status == libeasel.eslERANGE:
                    raise OverflowError("numerical overflow in the optimized vector implementation")
                elif status != libeasel.eslOK:
                    raise UnexpectedError(status, "p7_Pipeline")
                # clear pipeline and sequence for reuse for next target
                libhmmer.p7_pipeline.p7_pipeline_Reuse(pli)
                libeasel.sq.esl_sq_Reuse(dbsq)

        finally:
            libeasel.sq.esl_sq_Destroy(dbsq)

        # Return 0 to indicate success
        return 0

    cpdef TopHits scan_seq(
        self,
        DigitalSequence query,
        ScanTargets targets,
    ):
        """Run the pipeline using a query sequence against a profile database.

        Arguments:
            query (`~pyhmmer.easel.DigitalSequence`): The sequence object to
                use to query the profile database.
            targets (`OptimizedProfileBlock` or `HMMPressedFile`): The
                optimized profiles to query, either pre-loaded in memory as
                an `OptimizedProfileBlock`, or to be read iteratively from
                a pressed HMM file.

        Returns:
            `~pyhmmer.plan7.TopHits`: the hits found in the profile database.

        Raises:
            `~pyhmmer.errors.AlphabetMismatch`: When the alphabet of the
                current pipeline does not match the alphabet of the given
                query or profiles.

        Caution:
            In the current version, this method is not optimized to use
            the *pressed* database, even if it exists. This will cause the
            MSV and SSV filters to be rebuilt at each iteration, which could
            be slow. Consider at least pre-fetching the HMM database if
            calling this method several times in a row.

        Hint:
            This method corresponds to running ``hmmscan`` with the
            ``query`` sequence against the ``hmms`` database.

        .. versionadded:: 0.4.0

        .. versionchanged:: 0.7.0
            Require optimized profiles to be inside an `OptimizedProfileBlock`.

        """
        cdef int                  allocM
        cdef Profile              profile
        cdef TopHits              hits    = TopHits(query)

        assert self._pli != NULL

        # check the pipeline was configure with the same alphabet
        if not self.alphabet._eq(query.alphabet):
            raise AlphabetMismatch(self.alphabet, query.alphabet)
        if ScanTargets is OptimizedProfileBlock:
            if not self.alphabet._eq(targets.alphabet):
                raise AlphabetMismatch(self.alphabet, targets.alphabet)

        # verify the length
        if query._sq.n > HMMER_TARGET_LIMIT:
            raise ValueError("sequence length over comparison pipeline limit (100,000)")

        with nogil:
            # make sure the pipeline is set to scan mode and ready for a new sequence
            self._pli.mode = p7_pipemodes_e.p7_SCAN_MODELS
            self._pli.nmodels = 0
            # run the search loop on all database sequences while recycling memory
            if ScanTargets is OptimizedProfileBlock:
                Pipeline._scan_loop(
                    self._pli,
                    query._sq,
                    self.background._bg,
                    targets._block.list,
                    targets._block.count,
                    hits._th,
                    targets._locks,
                )
            else:
                Pipeline._scan_loop_file(
                    self._pli,
                    query._sq,
                    self.background._bg,
                    targets._hfp,
                    hits._th,
                )
            # threshold hits
            hits._sort_by_key()
            hits._threshold(self)

        # record the query metadata
        hits._query = query
        # return the hits
        return hits

    @staticmethod
    cdef int _scan_loop(
              P7_PIPELINE*        pli,
        const ESL_SQ*             sq,
              P7_BG*              bg,
              P7_OPROFILE**       om,
        const size_t              n_targets,
              P7_TOPHITS*         th,
              PyThread_type_lock* locks,
    ) except 1 nogil:
        cdef int    status
        cdef size_t t

        # configure the pipeline for the current sequence
        status = libhmmer.p7_pipeline.p7_pli_NewSeq(pli, sq)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_pli_NewSeq")

        # run the inner loop on all HMMs
        for t in range(n_targets):
            # configure the background and pipeline for the new optimized profile
            status = libhmmer.p7_pipeline.p7_pli_NewModel(pli, om[t], bg)
            if status == libeasel.eslEINVAL:
                Pipeline._missing_cutoffs(pli, om[t])
            elif status != libeasel.eslOK:
                raise UnexpectedError(status, "p7_pli_NewModel")
            status = libhmmer.p7_bg.p7_bg_SetLength(bg, sq.n)
            if status != libeasel.eslOK:
                raise UnexpectedError(status, "p7_bg_SetLength")
            # use a critical section here because `p7_oprofile_ReconfigLength`
            # will modify the optimized profile, so there is a risk of data race
            # if the scan loop is run in parallel on the same target profiles.
            # configure the profile
            if not PyThread_acquire_lock(locks[t], WAIT_LOCK):
                raise RuntimeError("Failed to acquire lock")
            try:
                status = p7_oprofile.p7_oprofile_ReconfigLength(om[t], sq.n)
                if status != libeasel.eslOK:
                    raise UnexpectedError(status, "p7_oprofile_ReconfigLength")
                # run the pipeline on the query sequence
                status = libhmmer.p7_pipeline.p7_Pipeline(pli, om[t], bg, sq, NULL, th)
                if status == libeasel.eslEINVAL:
                    Pipeline._missing_cutoffs(pli, om[t])
                elif status == libeasel.eslERANGE:
                    raise OverflowError("numerical overflow in the optimized vector implementation")
                elif status != libeasel.eslOK:
                    raise UnexpectedError(status, "p7_Pipeline")
            finally:
                PyThread_release_lock(locks[t])
            # clear pipeline for reuse for next target
            libhmmer.p7_pipeline.p7_pipeline_Reuse(pli)

        # Return 0 to indicate success
        return 0

    @staticmethod
    cdef int _scan_loop_file(
              P7_PIPELINE*     pli,
        const ESL_SQ*          sq,
              P7_BG*           bg,
              P7_HMMFILE*      hfp,
              P7_TOPHITS*      th,
    ) except 1 nogil:
        cdef int           status
        cdef P7_OPROFILE*  om
        cdef ESL_ALPHABET* abc    = <ESL_ALPHABET*> sq.abc

        # configure the pipeline for the current sequence
        status = libhmmer.p7_pipeline.p7_pli_NewSeq(pli, sq)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_pli_NewSeq")

        # run the inner loop on all HMMs
        while True:
            # read the next profile from the file
            status = p7_oprofile_ReadMSV(hfp, &abc, &om)
            if status == libeasel.eslOK:
                status = p7_oprofile_ReadRest(hfp, om)
            if status == libeasel.eslEOF:
                break
            elif status != libeasel.eslOK:
                raise UnexpectedError(status, "p7_oprofile_ReadMSV")
            # make sure that the optimized profile memory is not leaked
            try:
                # configure the background and pipeline for the new optimized profile
                status = libhmmer.p7_pipeline.p7_pli_NewModel(pli, om, bg)
                if status == libeasel.eslEINVAL:
                    Pipeline._missing_cutoffs(pli, om)
                elif status != libeasel.eslOK:
                    raise UnexpectedError(status, "p7_pli_NewModel")
                status = libhmmer.p7_bg.p7_bg_SetLength(bg, sq.n)
                if status != libeasel.eslOK:
                    raise UnexpectedError(status, "p7_bg_SetLength")
                # configure the profile
                status = p7_oprofile.p7_oprofile_ReconfigLength(om, sq.n)
                if status != libeasel.eslOK:
                    raise UnexpectedError(status, "p7_oprofile_ReconfigLength")
                # run the pipeline on the query sequence
                status = libhmmer.p7_pipeline.p7_Pipeline(pli, om, bg, sq, NULL, th)
                if status == libeasel.eslEINVAL:
                    Pipeline._missing_cutoffs(pli, om)
                elif status == libeasel.eslERANGE:
                    raise OverflowError("numerical overflow in the optimized vector implementation")
                elif status != libeasel.eslOK:
                    raise UnexpectedError(status, "p7_Pipeline")
            finally:
               p7_oprofile_Destroy(om)
               om = NULL

            # clear pipeline for reuse for next target
            libhmmer.p7_pipeline.p7_pipeline_Reuse(pli)

        # Return 0 to indicate success
        return 0

    cpdef IterativeSearch iterate_hmm(
        self,
        HMM query,
        DigitalSequenceBlock sequences,
        Builder builder = None,
        object select_hits = None,
    ):
        """Run the pipeline to find homologous sequences to a query HMM.

        Arguments:
            query (`~pyhmmer.plan7.HMM`): The sequence object to use to
                query the sequence database.
            sequences (`~pyhmmer.easel.DigitalSequenceBlock`): The target
                sequences to query with the HMM.
            builder (`~pyhmmer.plan7.Builder`, optional): A HMM builder to
                use to convert the query and subsequent alignments to a
                `~pyhmmer.plan7.HMM`. If `None` is given, this method will
                create one with the default parameters.
            select_hits (callable, optional): A function or callable object
                for manually selecting hits during each iteration. It should
                take a single `~pyhmmer.plan7.TopHits` argument and change the
                inclusion of every individual hits by setting the
                `~pyhmmer.plan7.Hit.included` and `~pyhmmer.plan7.Hit.dropped`
                flags of each `~pyhmmer.plan7.Hit` manually.

        Returns:
            `~pyhmmer.plan7.IterativeSearch`: An iterator object yielding
            the hits, sequence alignment, and HMM for each iteration.

        Raises:
            `~pyhmmer.errors.AlphabetMismatch`: When the alphabet of the
                current pipeline does not match the alphabet of the given
                query or database sequences.

        Hint:
            This method corresponds to running ``jackhmmer`` with the
            ``query`` sequence against the ``sequences`` database.

        Caution:
            Default values used for ``jackhmmer`` do not correspond to the
            default parameters used for creating a pipeline in the other
            cases. To have truly identical results to the ``jackhmmer``
            results in default mode, create the `Pipeline` object
            with ``incE=0.001`` and ``incdomE=0.001``.

        See Also:
            The `~Pipeline.iterate_seq`, which does the same operation with
            a query sequence instead of a query HMM, and contains more
            details and examples.

        .. versionchanged:: 0.7.0
            Targets must now be inside a `~pyhmmer.easel.DigitalSequenceBlock`.

        """
        # check that alphabets are consistent
        if not self.alphabet._eq(query.alphabet):
            raise AlphabetMismatch(self.alphabet, query.alphabet)
        if not self.alphabet._eq(sequences.alphabet):
            raise AlphabetMismatch(self.alphabet, sequences.alphabet)
        # check that builder is in hand architecture, not fast
        if builder is None:
            builder = Builder(self.alphabet, seed=self.seed, architecture="hand")
        elif builder.architecture != "hand":
            raise ValueError("`iterate_seq` only supports a builder with 'hand' architecture")
        # return the iterator
        return IterativeSearch(self, builder, query, sequences, select_hits)

    cpdef IterativeSearch iterate_seq(
        self,
        DigitalSequence query,
        DigitalSequenceBlock sequences,
        Builder builder = None,
        object select_hits = None,
    ):
        """Run the pipeline to find homologous sequences to a query sequence.

        This method implements an iterative search over a database to find
        all sequences homologous to a query sequence. It is very sensitive
        to the pipeline inclusion thresholds (``incE`` and ``incdomE``).

        Since this method returns an iterator, the local results of each
        iteration will be available for inspection before starting the next
        one. The ``select_hits`` callback in particular can be used for
        manually including/excluding hits in each iteration, which is not
        supported in the original ``jackhmmer``, but available on the HMMER
        `web client <https://www.ebi.ac.uk/Tools/hmmer/search/jackhmmer>`_.

        Arguments:
            query (`~pyhmmer.easel.DigitalSequence`): The sequence object to
                use to query the sequence database.
            sequences (`~pyhmmer.easel.DigitalSequenceBlock`): The target
                sequences to query with the query sequence.
            builder (`~pyhmmer.plan7.Builder`, optional): A HMM builder to
                use to convert the query and subsequent alignments to a
                `~pyhmmer.plan7.HMM`. If `None` is given, this method will
                create one with the default parameters.
            select_hits (callable, optional): A function or callable object
                for manually selecting hits during each iteration. It should
                take a single `TopHits` argument and change the inclusion of
                individual hits with the `Hit.include` and `Hit.drop` methods.

        Returns:
            `~pyhmmer.plan7.IterativeSearch`: An iterator object yielding
            the hits, sequence alignment, and HMM for each iteration.

        Raises:
            `~pyhmmer.errors.AlphabetMismatch`: When the alphabet of the
                current pipeline does not match the alphabet of the given
                query or database sequences.

        Hint:
            This method corresponds to running ``jackhmmer`` with the
            ``query`` sequence against the ``sequences`` database.

        Caution:
            Default values used for ``jackhmmer`` do not correspond to the
            default parameters used for creating a pipeline in the other
            cases. To have truly identical results to the ``jackhmmer``
            results in default mode, create the `Pipeline` object
            with ``incE=0.001`` and ``incdomE=0.001``.

        Example:
            Starting from a pipeline and a query sequence, let's first
            obtain the iterator over the successive results::

                >>> abc = easel.Alphabet.amino()
                >>> pli = plan7.Pipeline(abc, incE=1e-3, incdomE=1e-3)
                >>> iterator = pli.iterate_seq(reductase, proteins)

            Once this is ready, we can keep iterating until we converge::

                >>> converged = False
                >>> while not converged:
                ...     _, hits, _, converged, _ = next(iterator)
                ...     print(f"Hits: {len(hits)}  Converged: {converged}")
                Hits: 1  Converged: False
                Hits: 2  Converged: False
                Hits: 2  Converged: True

            To prevent diverging searches from running infinitely, you
            could wrap the search in a ``for`` loop instead, using a
            number of maximum iterations as the upper boundary::

                >>> iterator = pli.iterate_seq(reductase, proteins)
                >>> max_iterations = 10
                >>> for n in range(max_iterations):
                ...     iteration = next(iterator)
                ...     if iteration.converged:
                ...         break

        .. versionadded:: 0.6.0

        .. versionchanged:: 0.7.0
            Targets must now be inside a `~pyhmmer.easel.DigitalSequenceBlock`.

        """
        # check that alphabets are consistent
        if not self.alphabet._eq(query.alphabet):
            raise AlphabetMismatch(self.alphabet, query.alphabet)
        if not self.alphabet._eq(sequences.alphabet):
            raise AlphabetMismatch(self.alphabet, sequences.alphabet)
        # check that builder is in hand architecture, not fast
        if builder is None:
            builder = Builder(self.alphabet, seed=self.seed, architecture="hand")
        elif builder.architecture != "hand":
            raise ValueError("`iterate_seq` only supports a builder with 'hand' architecture")
        # return the iterator
        return IterativeSearch(self, builder, query, sequences, select_hits)


cdef double DEFAULT_LONG_F1           = 0.02
cdef double DEFAULT_LONG_F2           = 3e-3
cdef double DEFAULT_LONG_F3           = 3e-5
cdef int    DEFAULT_LONG_B1           = 100
cdef int    DEFAULT_LONG_B2           = 240
cdef int    DEFAULT_LONG_B3           = 1000
cdef int    DEFAULT_LONG_BLOCK_LENGTH = 0x40000

cdef class LongTargetsPipeline(Pipeline):
    """An HMMER3 pipeline tuned for long targets.

    The default HMMER3 pipeline is configured not to accept target sequences
    longer than 100,000 residues. Although there is no strong limitation for
    this threshold, comparing a sequence of :math:`L` residues to a profile
    with :math:`M` nodes requires the allocation of a :math:`L \times M`
    dynamic programming matrix.

    For sequences too long, it's actually more efficient memory-wise to use
    a sliding window to match the profile to the sequence. The usual
    comparison pipeline is then used to perform the comparison on each
    window, and results are merged once the entire sequence is done being
    processed. The context size :math:`C` is large enough to accommodate
    for the entire profile, so that there is no risk of missing a hit in
    the overlaps between windows. The window size :math:`W` can be
    changed with the ``block_length`` argument when instantiating a new
    `~pyhmmer.plan7.LongTargetsPipeline` object.

    .. versionadded:: 0.4.9

    .. versionadded:: 0.10.8
       The ``window_length`` and ``window_beta`` keyword arguments.

    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._window_length = -1
        self._idlens = NULL

    def __dealloc__(self):
        libhmmer.nhmmer.idlen_list_destroy(self._idlens)

    def __init__(
        self,
        Alphabet alphabet,
        Background background = None,
        *,
        double F1=DEFAULT_LONG_F1,
        double F2=DEFAULT_LONG_F2,
        double F3=DEFAULT_LONG_F3,
        str strand=None,
        int B1=DEFAULT_LONG_B1,
        int B2=DEFAULT_LONG_B2,
        int B3=DEFAULT_LONG_B3,
        int block_length=DEFAULT_LONG_BLOCK_LENGTH,
        object window_length=None,
        object window_beta=None,
        **kwargs,
    ):
        """__init__(self, alphabet, background=None, *, F1=0.02, F2=3e-3, F3=3e-5, strand=None, B1=100, B2=240, B3=1000, block_length=0x40000, window_length=None, window_beta=None, **kwargs)\n--\n

        Instantiate and configure a new long targets pipeline.

        Arguments:
            alphabet (`~pyhmmer.easel.Alphabet`): The biological alphabet the
                of the HMMs and sequences that are going to be compared. Used
                to build the background model. **A nucleotide alphabet is
                expected**.
            background (`~pyhmmer.plan7.Background`, optional): The background
                model to use with the pipeline, or ``None`` to create and use
                a default one. *The pipeline needs ownership of the background
                model, so any background model passed there will be copied.*

        Keyword Arguments:
            strand (`str`, optional): The strand to use when processing
                nucleotide sequences. Use ``"watson"`` to use only the
                coding strand, ``"crick"`` to use only the reverse strand,
                or leave as `None` to process both strands.
            B1 (`int`): The window length to use for the biased-composition
                modifier of the MSV filter.
            B2 (`int`): The window length to use for the biased-composition
                modifier of the Viterbi filter.
            B3 (`int`): The window length to use for the biased-composition
                modifier of the Forward filter.
            block_length (`int`): The number of residues to use as the
                window size :math:`W` when reading blocks from the long
                target sequences.
            window_length (`int`): The window length to use to compute
                E-values.
            window_beta (`float`): The tail mass at which window length
                is determined.
            **kwargs: Any additional parameter will be passed to the
                `~pyhmmer.plan7.Pipeline` constructor.

        """
        # check that a nucleotide alphabet is given
        if not alphabet.is_nucleotide():
            raise AlphabetMismatch(Alphabet.dna(), alphabet)

        # create the pipeline
        super().__init__(alphabet, background, F1=F1, F2=F2, F3=F3, **kwargs)

        # set the options for long targets
        self._pli.long_targets = True
        self._pli.block_length = block_length
        self.strand = strand
        self.B1 = B1
        self.B2 = B2
        self.B3 = B3
        self.window_length = window_length
        if window_beta is None:
            self.window_beta = libhmmer.p7_builder.p7_DEFAULT_WINDOW_BETA
        else:
            self.window_beta = window_beta

        # allocate storage for the ID/length list
        self._idlens = libhmmer.nhmmer.idlen_list_init(64)
        if self._idlens == NULL:
            raise AllocationError("ID_LENGTH_LIST", sizeof(ID_LENGTH_LIST))

    # --- Properties ---------------------------------------------------------

    @property
    def B1(self):
        """`int`: The window length for biased-composition modifier of the MSV.
        """
        assert self._pli != NULL
        return self._pli.B1

    @B1.setter
    def B1(self, int B1):
        assert self._pli != NULL
        self._pli.B1 = B1

    @property
    def B2(self):
        """`int`: The window length for biased-composition modifier of the Viterbi.
        """
        assert self._pli != NULL
        return self._pli.B2

    @B2.setter
    def B2(self, int B2):
        assert self._pli != NULL
        self._pli.B2 = B2

    @property
    def B3(self):
        """`int`: The window length for biased-composition modifier of the Forward.
        """
        assert self._pli != NULL
        return self._pli.B3

    @B3.setter
    def B3(self, int B3):
        assert self._pli != NULL
        self._pli.B3 = B3

    @property
    def strand(self):
        """`str` or `None`: The strand to process, or `None` for both.
        """
        assert self._pli != NULL
        return (None, "watson", "crick")[self._pli.strands]

    @strand.setter
    def strand(self, str strand):
        assert self._pli != NULL
        if strand is None:
            self._pli.strands = p7_strands_e.p7_STRAND_BOTH
        elif strand == "watson":
            self._pli.strands = p7_strands_e.p7_STRAND_TOPONLY
        elif strand == "crick":
            self._pli.strands = p7_strands_e.p7_STRAND_BOTTOMONLY
        else:
            raise InvalidParameter("strand", strand, choices=["watson", "crick", None])

    @property
    def window_length(self):
        """`int` or `None`: The window length for nucleotide sequences.
        """
        return None if self._window_length == -1 else self._window_length

    @window_length.setter
    def window_length(self, object window_length):
        if window_length is None:
            self._window_length = -1
        elif window_length > 0:
            self._window_length = window_length
        else:
            raise InvalidParameter("window_length", window_length, hint="integer greater than 0 or None")

    @property
    def window_beta(self):
        """`float`: The tail mass at which window length is determined.
        """
        return self._window_beta

    @window_beta.setter
    def window_beta(self, double window_beta):
        if window_beta > 1 or window_beta < 0:
            raise InvalidParameter("window_beta", window_beta, hint="real number between 0 and 1")
        self._window_beta = window_beta

    # --- Methods ------------------------------------------------------------

    cpdef list arguments(self):
        """Generate an ``argv``-like list with pipeline options as CLI flags.

        Example:
            >>> alphabet = easel.Alphabet.dna()
            >>> plan7.LongTargetsPipeline(alphabet).arguments()
            []
            >>> plan7.LongTargetsPipeline(alphabet, B1=200).arguments()
            ['--B1', '200']

        .. versionadded:: 0.6.0

        """
        cdef list argv = []
        cdef str  cutoff

        if self._pli.use_bit_cutoffs:
            cutoff = self.bit_cutoffs
            if cutoff == "gathering":
                argv.append("--cut_ga")
            elif cutoff == "noise":
                argv.append("--cut_nc")
            elif cutoff == "trusted":
                argv.append("--cut_tc")
            else:
                raise InvalidParameter("bit_cutoffs", cutoff, choices=["gathering", "trusted", "noise"])
        else:
            # options controlling reporting thresholds
            if self._pli.by_E:
                if self._pli.E != DEFAULT_E:
                    argv.append("-E")
                    argv.append(str(self._pli.E))
            else:
                argv.append("-T")
                argv.append(str(self._pli.T))
            if self._pli.dom_by_E:
                if self._pli.domE != DEFAULT_DOME:
                    argv.append("--domE")
                    argv.append(str(self._pli.domE))
            else:
                argv.append("--domT")
                argv.append(str(self._pli.domT))
            # options controlling inclusion thresholds
            if self._pli.inc_by_E:
                if self._pli.incE != DEFAULT_INCE:
                    argv.append("--incE")
                    argv.append(str(self._pli.incE))
            else:
                argv.append("--incT")
                argv.append(str(self._pli.incT))
            if self._pli.incdom_by_E:
                if self._pli.incdomE != DEFAULT_INCDOME:
                    argv.append("--incdomE")
                    argv.append(str(self._pli.incdomE))
            else:
                argv.append("--incdomT")
                argv.append(str(self._pli.incdomT))

        if self._pli.F1 != DEFAULT_LONG_F1:
            argv.append("--F1")
            argv.append(str(self.F1))
        if self._pli.F2 != DEFAULT_LONG_F2:
            argv.append("--F2")
            argv.append(str(self.F2))
        if self._pli.F3 != DEFAULT_LONG_F3:
            argv.append("--F3")
            argv.append(str(self.F3))

        if self._pli.B1 != DEFAULT_LONG_B1:
            argv.append("--B1")
            argv.append(str(self.B1))
        if self._pli.B2 != DEFAULT_LONG_B2:
            argv.append("--B2")
            argv.append(str(self.B2))
        if self._pli.B3 != DEFAULT_LONG_B3:
            argv.append("--B3")
            argv.append(str(self.B3))

        if not self.bias_filter:
            argv.append("--nobias")
        if not self.null2:
            argv.append("--nonull2")

        if self._pli.Z_setby == p7_zsetby_e.p7_ZSETBY_OPTION:
            argv.append("-Z")
            argv.append(str(self._pli.Z))
        if self._pli.domZ_setby == p7_zsetby_e.p7_ZSETBY_OPTION:
            argv.append("--domZ")
            argv.append(str(self._pli.domZ))

        if self._seed != DEFAULT_SEED:
            argv.append("--seed")
            argv.append(str(self._seed))

        if self._pli.block_length != DEFAULT_LONG_BLOCK_LENGTH:
            argv.append("--block_length")
            argv.append(str(self._pli.block_length))

        if self._pli.strands == p7_strands_e.p7_STRAND_TOPONLY:
            argv.append("--watson")
        elif self._pli.strands == p7_strands_e.p7_STRAND_BOTTOMONLY:
            argv.append("--crick")

        return argv

    cpdef TopHits scan_seq(
        self,
        DigitalSequence query,
        ScanTargets targets,
    ):
        """Run the pipeline using a query sequence against a profile database.

        This is currently unsupported for `LongTargetsPipeline`.

        """
        raise NotImplementedError("Cannot run a database scan with the long target pipeline")



    cpdef TopHits search_hmm(
        self,
        object query,
        SearchTargets sequences
    ):
        """Run the pipeline using a query HMM against a sequence database.

        Arguments:
            query (`HMM`, `Profile` or `OptimizedProfile`): The object to use
                to query the sequence database.
            sequences (`~pyhmmer.easel.DigitalSequenceBlock`): The target
                sequences to query with the HMM.

        Returns:
            `~pyhmmer.plan7.TopHits`: the hits found in the sequence database.

        Raises:
            `ValueError`: When the pipeline is configured to use model-specific
                reporting thresholds but the `HMM` query doesn't have the right
                cutoffs available.
            `~pyhmmer.errors.AlphabetMismatch`: When the alphabet of the
                current pipeline does not match the alphabet of the given
                HMM.

        Hint:
            This method corresponds to running ``nhmmer`` with the
            ``query`` HMM against the ``sequences`` database.

        .. versionchanged:: 0.7.0
            Targets must now be inside a `~pyhmmer.easel.DigitalSequenceBlock`.

        """
        assert self._pli != NULL

        cdef size_t               L
        cdef uint64_t             j
        cdef ssize_t              nseqs
        cdef int                  status
        cdef int64_t              res_count
        cdef int                  allocM
        cdef HMM                  hmm
        cdef int                  max_length
        cdef ScoreData            scoredata      = ScoreData.__new__(ScoreData)
        cdef TopHits              hits           = TopHits(query)
        cdef P7_HIT*              hit            = NULL
        cdef P7_OPROFILE*         om             = NULL

        # perform specific checks for targets stored in a file
        if SearchTargets is SequenceFile:
            if sequences.name is None:
                raise ValueError("can only use a `SequenceFile` backed by a file for reading targets")
            if not sequences.digital:
                raise ValueError("target sequences file is not in digital mode")

        # check the pipeline was configured with the same alphabet
        if not self.alphabet._eq(query.alphabet):
            raise AlphabetMismatch(self.alphabet, query.alphabet)
        if not self.alphabet._eq(sequences.alphabet):
            raise AlphabetMismatch(self.alphabet, sequences.alphabet)

        # get a length estimate for profile configuration
        if SearchTargets is DigitalSequenceBlock:
            L = self.L_HINT if not sequences else sequences._refs[0].L
        else:
            L = self.L_HINT
        # get the optimized profile from the query
        om = self._get_om_from_query(query, L=L)

        # compute max length based on window length to use for E-values
        max_length = om.max_length
        if self._window_length > 0:
            max_length = om.max_length = self._window_length
        elif isinstance(query, HMM):
            hmm = query
            libhmmer.p7_builder.p7_Builder_MaxLength(hmm._hmm, self._window_beta)
            max_length = hmm._hmm.max_length
        else:
            raise TypeError("Cannot use `Profile` or `OptimizedProfile` query without `max_length` set")

        with nogil:
            # make sure the pipeline is set to search mode and ready for a new HMM
            self._pli.mode = p7_pipemodes_e.p7_SEARCH_SEQS
            self._pli.nseqs = 0
            # reset the ID/length list
            libhmmer.nhmmer.idlen_list_clear(self._idlens)
            # create the score data struct
            scoredata.Kp = self.profile._gm.abc.Kp
            scoredata._sd = libhmmer.p7_scoredata.p7_hmm_ScoreDataCreate(om, NULL)
            if scoredata._sd == NULL:
                raise AllocationError("P7_SCOREDATA", sizeof(P7_SCOREDATA))
            # run the search loop on all database sequences while recycling memory
            if SearchTargets is DigitalSequenceBlock:
                LongTargetsPipeline._search_loop_longtargets(
                    self._pli,
                    om,
                    self.background._bg,
                    <const ESL_SQ**> sequences._refs,
                    sequences._length,
                    hits._th,
                    scoredata._sd,
                    self._idlens,
                )
            else:
                LongTargetsPipeline._search_loop_longtargets_file(
                    self._pli,
                    om,
                    self.background._bg,
                    sequences._sqfp,
                    hits._th,
                    scoredata._sd,
                    self._idlens,
                )
            # threshold with user-provided Z if any and compute E-values
            if self._Z is None:
                res_count = self._pli.nres
            else:
                res_count = <int64_t> (1000000 * self._pli.Z)
                if self._pli.strands == p7_strands_e.p7_STRAND_BOTH:
                    res_count *= 2
            libhmmer.p7_tophits.p7_tophits_ComputeNhmmerEvalues(hits._th, res_count, max_length)
            # assign target sequence lengths
            hits._sort_by_seqidx()
            status = libhmmer.nhmmer.idlen_list_assign(self._idlens, hits._th)
            if status != libeasel.eslOK:
                raise UnexpectedError(status, "idlen_list_assign")
            # remove target duplicates from the hits (using the ones with best E-value)
            libhmmer.p7_tophits.p7_tophits_RemoveDuplicates(hits._th, True)
            # sort hits, threshold and record pipeline configuration
            hits._sort_by_key()
            hits._threshold(self)
            # tally up total number of hits and target coverage
            hits._pli.n_output = hits._pli.pos_output = 0
            for j in range(hits._th.N):
                if (hits._th.hit[j].flags & p7_hitflags_e.p7_IS_REPORTED) or (hits._th.hit[j].flags & p7_hitflags_e.p7_IS_INCLUDED):
                    hits._pli.n_output += 1
                    hits._pli.pos_output += 1 + llabs(hits._th.hit[j].dcl[0].jali - hits._th.hit[j].dcl[0].iali)

        # record the query metadata
        hits._query = query
        # return the hits
        return hits

    cpdef TopHits search_seq(
        self,
        DigitalSequence query,
        object sequences,
        Builder builder = None,
    ):
        """Run the pipeline using a query sequence against a sequence database.

        Arguments:
            query (`~pyhmmer.easel.DigitalSequence`): The sequence object to
                use to query the sequence database.
            sequences (`~pyhmmer.easel.DigitalSequenceBlock`): The target
                sequences to query with the query sequence.
            builder (`~pyhmmer.plan7.Builder`, optional): A HMM builder to
                use to convert the query to a `~pyhmmer.plan7.HMM`. If
                `None` is given, it will use a default one.

        Returns:
            `~pyhmmer.plan7.TopHits`: the hits found in the sequence database.

        Raises:
            `~pyhmmer.errors.AlphabetMismatch`: When the alphabet of the
                current pipeline does not match the alphabet of the given
                query.

        Hint:
            This method corresponds to running ``nhmmer`` with the
            ``query`` sequence against the ``sequences`` database.

        .. versionchanged:: 0.7.0
            Targets must now be inside a `~pyhmmer.easel.DigitalSequenceBlock`.

        """
        assert self._pli != NULL

        cdef HMM              hmm
        cdef TopHits          hits

        if not self.alphabet._eq(query.alphabet):
            raise AlphabetMismatch(self.alphabet, query.alphabet)

        if builder is None:
            builder = Builder(
                self.alphabet,
                seed=self.seed,
                window_length=self.window_length,
                window_beta=self.window_beta,
            )
        elif builder.window_length != self.window_length:
            raise ValueError("builder and long targets pipeline have different window lengths")
        elif builder.window_beta != self.window_beta:
            raise ValueError("builder and long targets pipeline have different window beta")

        hmm = builder.build(query, self.background)[0]
        assert hmm._hmm.max_length != -1
        if isinstance(sequences, DigitalSequenceBlock):
            hits = self.search_hmm[DigitalSequenceBlock](hmm, sequences)
        elif isinstance(sequences, SequenceFile):
            hits = self.search_hmm[SequenceFile](hmm, sequences)
        else:
            ty = type(sequences).__name__
            raise TypeError(f"Expected DigitalSequenceBlock or SequenceFile, found {ty}")

        hits._query = query
        return hits

    cpdef TopHits search_msa(
        self,
        DigitalMSA query,
        object sequences,
        Builder builder = None,
    ):
        """Run the pipeline using a query alignment against a sequence database.

        Arguments:
            query (`~pyhmmer.easel.DigitalMSA`): The multiple sequence
                alignment to use to query the sequence database.
            sequences (`~pyhmmer.easel.DigitalSequenceBlock`): The target
                sequences to query with the query alignment.
            builder (`~pyhmmer.plan7.Builder`, optional): A HMM builder to
                use to convert the query to a `~pyhmmer.plan7.HMM`. If
                `None` is given, it will use a default one.

        Returns:
            `~pyhmmer.plan7.TopHits`: the hits found in the sequence database.

        Raises:
            `~pyhmmer.errors.AlphabetMismatch`: When the alphabet of the
                current pipeline does not match the alphabet of the given
                query.

        Hint:
            This method corresponds to running ``nhmmer`` with the
            ``query`` multiple sequence alignment against the ``sequences``
            database.

        .. versionchanged:: 0.7.0
            Targets must now be inside a `~pyhmmer.easel.DigitalSequenceBlock`.

        """
        assert self._pli != NULL

        cdef HMM              hmm
        cdef TopHits          hits

        if not self.alphabet._eq(query.alphabet):
            raise AlphabetMismatch(self.alphabet, query.alphabet)

        if builder is None:
            builder = Builder(
                self.alphabet,
                seed=self.seed,
                window_length=self.window_length,
                window_beta=self.window_beta,
            )
        elif builder.window_length != self.window_length:
            raise ValueError("builder and long targets pipeline have different window lengths")
        elif builder.window_beta != self.window_beta:
            raise ValueError("builder and long targets pipeline have different window beta")

        hmm = builder.build_msa(query, self.background)[0]
        assert hmm._hmm.max_length != -1
        if isinstance(sequences, DigitalSequenceBlock):
            hits = self.search_hmm[DigitalSequenceBlock](hmm, sequences)
        elif isinstance(sequences, SequenceFile):
            hits = self.search_hmm[SequenceFile](hmm, sequences)
        else:
            ty = type(sequences).__name__
            raise TypeError(f"Expected DigitalSequenceBlock or SequenceFile, found {ty}")

        hits._query = query
        return hits

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
    ) except 1 nogil:
        cdef int     status
        cdef int     nres
        cdef size_t  t
        cdef int64_t j      = 0
        cdef int64_t rem    = 0
        cdef int64_t i      = 0
        cdef int64_t C      = om.max_length
        cdef int64_t W      = pli.block_length

        if C <= 0:
            raise InvalidParameter("C", C, hint="strictly positive integer")
        if W <= 0:
            raise InvalidParameter("W", W, hint="strictly positive integer greater than C")
        if W <= C:
            # TODO: see if this is handled in nhmmer, and how to emulate it
            raise InvalidParameter("W", W, hint="strictly positive integer greater than C")

        # configure the pipeline for the current HMM
        status = libhmmer.p7_pipeline.p7_pli_NewModel(pli, om, bg)
        if status == libeasel.eslEINVAL:
            Pipeline._missing_cutoffs(pli, om)
        elif status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_pli_NewModel")

        # create a temporary sequence in digital mode to store the current window
        tmpsq = libeasel.sq.esl_sq_CreateDigital(om.abc)
        if tmpsq == NULL:
            raise AllocationError("ESL_SQ", sizeof(ESL_SQ))

        # run the inner loop on all sequences
        for t in range(n_targets):
            # record length
            libhmmer.nhmmer.idlen_list_add(idlens, t, sq[t].L)

            # initialize the sequence window storage
            tmpsq.idx = t
            tmpsq.L = -1
            libeasel.sq.esl_sq_SetAccession(tmpsq, sq[t].acc)
            libeasel.sq.esl_sq_SetName(tmpsq, sq[t].name)
            libeasel.sq.esl_sq_SetDesc(tmpsq, sq[t].desc)
            libeasel.sq.esl_sq_SetSource(tmpsq, sq[t].name)
            libeasel.sq.esl_sq_GrowTo(tmpsq, min(W+C, sq[t].n))

            # iterate over successive windows of width W, keeping C residues
            # from the previous iteration as context
            # NB: this is basically the PyRex equivalent for this:
            #     ```
            #     for i in range(0, sq[n].n, W-C)
            #     ```
            #     but Cython refuses to optimize that because it cannot be
            #     sure that W-C is strictly positive, even though we do check
            #     that it is.
            for i from 0 <= i < sq[t].n by W - C:
                # update window coordinates in the temporary sequence
                tmpsq.C     = 0 if i == 0 else min(C, sq[t].n - i)
                tmpsq.W     = min(W, sq[t].n - i - tmpsq.C)
                tmpsq.n     = tmpsq.C + tmpsq.W
                tmpsq.start = i + 1
                tmpsq.end   = i + tmpsq.n
                # DEBUG: show loop state
                # printf("[i=%li] C=%li W=%li n=%li start=%li end=%li\n", i, tmpsq.C, tmpsq.W, tmpsq.n, tmpsq.start, tmpsq.end);

                # copy sequence digits from the target sequence to the temporary sequence
                memcpy(&tmpsq.dsq[1], &sq[t].dsq[i+1], tmpsq.n)
                tmpsq.dsq[0] = tmpsq.dsq[tmpsq.n+1] = libeasel.eslDSQ_SENTINEL

                # configure the profile, background and pipeline for the new sequence
                status = libhmmer.p7_pipeline.p7_pli_NewSeq(pli, tmpsq)
                if status != libeasel.eslOK:
                    raise UnexpectedError(status, "p7_pli_NewSeq")

                # process coding strand
                if pli.strands != p7_strands_e.p7_STRAND_BOTTOMONLY:
                    # account for overlapping region of windows
                    pli.nres -= tmpsq.C
                    # run the pipeline on the forward strand
                    status = libhmmer.p7_pipeline.p7_Pipeline_LongTarget(pli, om, scoredata, bg, th, pli.nseqs, tmpsq, p7_complementarity_e.p7_NOCOMPLEMENT, NULL, NULL, NULL)
                    if status == libeasel.eslEINVAL:
                        Pipeline._missing_cutoffs(pli, om)
                    elif status == libeasel.eslERANGE:
                        raise OverflowError("numerical overflow in the optimized vector implementation")
                    elif status != libeasel.eslOK:
                        raise UnexpectedError(status, "p7_Pipeline_LongTarget")
                    # clear pipeline for reuse for next target
                    libhmmer.p7_pipeline.p7_pipeline_Reuse(pli)
                else:
                    pli.nres -= tmpsq.n

                # process reverse strand
                if pli.strands != p7_strands_e.p7_STRAND_TOPONLY:
                    # reverse complement sequence
                    libeasel.sq.esl_sq_ReverseComplement(tmpsq)
                    # run the pipeline on the reverse strand
                    status = libhmmer.p7_pipeline.p7_Pipeline_LongTarget(pli, om, scoredata, bg, th, pli.nseqs, tmpsq, p7_complementarity_e.p7_COMPLEMENT, NULL, NULL, NULL)
                    if status == libeasel.eslEINVAL:
                        Pipeline._missing_cutoffs(pli, om)
                    elif status == libeasel.eslERANGE:
                        raise OverflowError("numerical overflow in the optimized vector implementation")
                    elif status != libeasel.eslOK:
                        raise UnexpectedError(status, "p7_Pipeline_LongTarget")
                    # clear pipeline for reuse for next target
                    libhmmer.p7_pipeline.p7_pipeline_Reuse(pli)
                    pli.nres += tmpsq.W

            # clear the allocated sequence and advance to next sequence
            libeasel.sq.esl_sq_Reuse(tmpsq)
            pli.nseqs += 1

        # Free temporary data
        libeasel.sq.esl_sq_Destroy(tmpsq)

        # Return 0 to indicate success
        return 0
    @staticmethod
    cdef int _search_loop_longtargets_file(
        P7_PIPELINE*    pli,
        P7_OPROFILE*    om,
        P7_BG*          bg,
        ESL_SQFILE*     sqfp,
        P7_TOPHITS*     th,
        P7_SCOREDATA*   scoredata,
        ID_LENGTH_LIST* idlens
    ) except 1 nogil:
        cdef int     status
        #cdef int     nres
        #cdef size_t  t
        cdef int     seq_id  = 0
        #cdef int64_t index   = 0
        #cdef int64_t j       = 0
        #cdef int64_t rem     = 0
        #cdef int64_t i       = 0
        cdef int64_t C       = om.max_length
        cdef int64_t W       = pli.block_length
        cdef ESL_SQ* dbsq    = NULL
        cdef ESL_SQ* dbsq_rc = NULL

        if C <= 0:
            raise InvalidParameter("C", C, hint="strictly positive integer")
        if W <= 0:
            raise InvalidParameter("W", W, hint="strictly positive integer greater than C")
        if W <= C:
            # TODO: see if this is handled in nhmmer, and how to emulate it?
            raise InvalidParameter("W", W, hint="strictly positive integer greater than C")

        # configure the pipeline for the current HMM
        status = libhmmer.p7_pipeline.p7_pli_NewModel(pli, om, bg)
        if status == libeasel.eslEINVAL:
            Pipeline._missing_cutoffs(pli, om)
        elif status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_pli_NewModel")

        try:
            # allocate temporary memory for the sequences
            dbsq = libeasel.sq.esl_sq_CreateDigital(om.abc)
            if dbsq == NULL:
                raise AllocationError("ESL_SQ", sizeof(ESL_SQ))
            if om.abc.complement != NULL:
                dbsq_rc = libeasel.sq.esl_sq_CreateDigital(om.abc)
                if dbsq_rc == NULL:
                    raise AllocationError("ESL_SQ", sizeof(ESL_SQ))
            # read the first sequence window
            status = libeasel.sqio.esl_sqio_ReadWindow(sqfp, 0, W, dbsq)
            if status == libeasel.eslEOF:
                return 0
            elif status != libeasel.eslOK:
                raise UnexpectedError(status, "esl_sqio_ReadWindow")
            # run the search loop
            while status == libeasel.eslOK:
                # reconfigure the pipeline
                dbsq.idx = seq_id
                libhmmer.p7_pipeline.p7_pli_NewSeq(pli, dbsq)
                # process forward strand
                if pli.strands != p7_strands_e.p7_STRAND_BOTTOMONLY:
                    pli.nres -= dbsq.C
                    status = libhmmer.p7_pipeline.p7_Pipeline_LongTarget(pli, om, scoredata, bg, th, pli.nseqs, dbsq, p7_complementarity_e.p7_NOCOMPLEMENT, NULL, NULL, NULL)
                    if status == libeasel.eslEINVAL:
                        Pipeline._missing_cutoffs(pli, om)
                    elif status == libeasel.eslERANGE:
                        raise OverflowError("numerical overflow in the optimized vector implementation")
                    elif status != libeasel.eslOK:
                        raise UnexpectedError(status, "p7_Pipeline_LongTarget")
                    libhmmer.p7_pipeline.p7_pipeline_Reuse(pli)
                else:
                    pli.nres -= dbsq.n
                # process reverse strand
                if pli.strands != p7_strands_e.p7_STRAND_BOTTOMONLY and dbsq.abc.complement != NULL:
                    libeasel.sq.esl_sq_Copy(dbsq, dbsq_rc)
                    libeasel.sq.esl_sq_ReverseComplement(dbsq_rc)
                    libhmmer.p7_pipeline.p7_Pipeline_LongTarget(pli, om, scoredata, bg, th, pli.nseqs, dbsq_rc, p7_complementarity_e.p7_COMPLEMENT, NULL, NULL, NULL)
                    if status == libeasel.eslEINVAL:
                        Pipeline._missing_cutoffs(pli, om)
                    elif status == libeasel.eslERANGE:
                        raise OverflowError("numerical overflow in the optimized vector implementation")
                    elif status != libeasel.eslOK:
                        raise UnexpectedError(status, "p7_Pipeline_LongTarget")
                    libhmmer.p7_pipeline.p7_pipeline_Reuse(pli)
                    pli.nres += dbsq_rc.W
                # read next window
                status = libeasel.sqio.esl_sqio_ReadWindow(sqfp, C, W, dbsq)
                if status == libeasel.eslEOD:
                    libhmmer.nhmmer.idlen_list_add(idlens, dbsq.idx, dbsq.L)
                    pli.nseqs += 1
                    libeasel.sq.esl_sq_Reuse(dbsq)
                    status = libeasel.sqio.esl_sqio_ReadWindow(sqfp, 0, W, dbsq)
                    seq_id += 1
                if status != libeasel.eslOK and status != libeasel.eslEOF:
                    raise UnexpectedError(status, "esl_sqio_ReadWindow")
        finally:
            libeasel.sq.esl_sq_Destroy(dbsq)
            libeasel.sq.esl_sq_Destroy(dbsq_rc)

        return 0


cdef class Profile:
    """A Plan7 search profile.

    Attributes:
        alphabet (`~pyhmmer.easel.Alphabet`): The alphabet for which this
            profile was configured.

    .. versionchanged:: 0.4.6
       Added the `~Profile.evalue_parameters` and `~Profile.cutoffs` attributes.

    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._gm = NULL
        self.alphabet = None

    def __init__(self, int M, Alphabet alphabet):
        """__init__(self, M, alphabet)\n--\n

        Create a new profile for the given ``alphabet``.

        Arguments:
            M (`int`): The length of the profile, i.e. the number of nodes.
            alphabet (`~pyhmmer.easel.Alphabet`): The alphabet to use with
                this profile.

        """
        # store the alphabet to make sure it's not deallocated
        self.alphabet = alphabet
        # create a new profile large enough to store M nodes
        with nogil:
            self._gm = libhmmer.p7_profile.p7_profile_Create(M, alphabet._abc)
        if not self._gm:
            raise AllocationError("P7_PROFILE", sizeof(P7_PROFILE))

    def __dealloc__(self):
        libhmmer.p7_profile.p7_profile_Destroy(self._gm)

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"<{ty} alphabet={self.alphabet!r} M={self.M!r} name={self.name!r}>"

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        cdef Profile new
        if id(self) not in memo:
            new = memo[id(self)] = self.copy()
            new.alphabet = copy.deepcopy(self.alphabet, memo=memo)
            new._gm.abc = new.alphabet._abc
        return memo[id(self)]

    def __eq__(self, object other):
        assert self._gm != NULL

        if not isinstance(other, Profile):
            return NotImplemented

        cdef Profile      p = <Profile> other
        cdef int     status = libhmmer.p7_profile.p7_profile_Compare(self._gm, p._gm, 0.0)

        if status == libeasel.eslOK:
            return True
        elif status == libeasel.eslFAIL:
            return False
        else:
            raise UnexpectedError(status, "p7_profile_Compare")

    def __sizeof__(self):
        assert self._gm != NULL
        return libhmmer.p7_profile.p7_profile_Sizeof(self._gm) + sizeof(self)

    # --- Properties ---------------------------------------------------------

    @property
    def M(self):
        """`int`: The length of the profile (i.e. the number of nodes).

        .. versionadded:: 0.3.0

        """
        assert self._gm != NULL
        return self._gm.M

    @property
    def L(self):
        """`int`: The current configured target sequence length.
        """
        assert self._gm != NULL
        return self._gm.L

    @L.setter
    def L(self, int L):
        assert self._gm != NULL
        cdef int status
        with nogil:
            status = libhmmer.modelconfig.p7_ReconfigLength(self._gm, L)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_ReconfigLength")

    @property
    def accession(self):
        """`bytes` or `None`: The accession of the profile, if any.

        .. versionadded:: 0.3.0

        """
        assert self._gm != NULL
        return None if self._gm.acc == NULL else <bytes> self._gm.acc

    @property
    def name(self):
        """`bytes` or `None`: The name of the profile, if any.

        .. versionadded:: 0.3.0

        """
        assert self._gm != NULL
        return None if self._gm.name == NULL else <bytes> self._gm.name

    @property
    def description(self):
        """`bytes` or `None`: The description of the profile, if any.

        .. versionadded:: 0.3.0

        """
        assert self._gm != NULL
        return None if self._gm.desc == NULL else <bytes> self._gm.desc

    @property
    def consensus(self):
        """`str` or `None`: The consensus residue line of the profile, if set.

        .. versionadded:: 0.4.1

        """
        assert self._gm != NULL
        if self._gm.consensus[0] == b'\0':
            return None
        return (&self._gm.consensus[1]).decode("ascii")

    @property
    def consensus_structure(self):
        """`str` or `None`: The consensus structure of the profile, if any.

        .. versionadded:: 0.4.1

        """
        assert self._gm != NULL
        if self._gm.cs[0] == b'\0':
            return None
        return (&self._gm.cs[1]).decode("ascii")

    @property
    def offsets(self):
        """`~pyhmmer.plan7.Offsets`: The disk offsets for this profile.
        """
        assert self._gm != NULL
        cdef Offsets offsets = Offsets.__new__(Offsets)
        offsets._offs = &self._gm.offs
        offsets._owner = self
        return offsets

    @property
    def evalue_parameters(self):
        """`~plan7.EvalueParameters`: The e-value parameters for this profile.
        """
        assert self._gm != NULL
        cdef EvalueParameters ep = EvalueParameters.__new__(EvalueParameters)
        ep._evparams = &self._gm.evparam
        ep._owner = self
        return ep

    @property
    def cutoffs(self):
        """`~plan7.Cutoffs`: The bitscore cutoffs for this profile.
        """
        assert self._gm != NULL
        cdef Cutoffs cutoffs = Cutoffs.__new__(Cutoffs)
        cutoffs._owner = self
        cutoffs._cutoffs = &self._gm.cutoff
        cutoffs._flags = NULL
        cutoffs._is_profile = True
        return cutoffs

    @property
    def local(self):
        """`bool`: Whether the profile is in local mode.

        .. versionadded:: 0.7.0

        """
        assert self._gm != NULL
        return libhmmer.p7_profile.p7_profile_IsLocal(self._gm)

    @property
    def multihit(self):
        """`bool`: Whether the profile is in multihit mode.

        .. versionadded:: 0.7.0

        """
        assert self._gm != NULL
        return libhmmer.p7_profile.p7_profile_IsMultihit(self._gm)

    @multihit.setter
    def multihit(self, multihit):
        if multihit:
            if not self.multihit:
                libhmmer.modelconfig.p7_ReconfigMultihit(self._gm, self._gm.L)
        else:
            if self.multihit:
                libhmmer.modelconfig.p7_ReconfigUnihit(self._gm, self._gm.L)

    # --- Methods ------------------------------------------------------------

    cpdef void clear(self) except *:
        """Clear internal buffers to reuse the profile without reallocation.
        """
        cdef int status
        with nogil:
            status = libhmmer.p7_profile.p7_profile_Reuse(self._gm)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_profile_Reuse")

    cpdef void configure(
        self,
        HMM hmm,
        Background background,
        int L=400,
        bint multihit=True,
        bint local=True
    ) except *:
        """Configure a search profile using the given models.

        Arguments:
            hmm (`~pyhmmer.plan7.HMM`): The model HMM with core probabilities.
            background (`~pyhmmer.plan7.Background`): The null background model.
            L (`int`): The expected target sequence length.
            multihit (`bool`): Whether or not to use multihit modes.
            local (`bool`): Whether or not to use non-local modes.

        """
        assert self._gm != NULL

        cdef int         mode
        cdef int         status
        cdef P7_HMM*     hm     = hmm._hmm
        cdef P7_BG*      bg     = background._bg
        cdef P7_PROFILE* gm     = self._gm

        # check alphabets are correct
        if not self.alphabet._eq(hmm.alphabet):
            raise AlphabetMismatch(self.alphabet, hmm.alphabet)
        elif not self.alphabet._eq(background.alphabet):
            raise AlphabetMismatch(self.alphabet, background.alphabet)
        elif hm.M > self._gm.allocM:
            raise ValueError("Profile is too small to hold HMM")
        # get the proper flag
        if multihit:
            mode = p7_LOCAL if local else p7_GLOCAL
        else:
            mode = p7_UNILOCAL if local else p7_UNIGLOCAL
        # configure the model
        with nogil:
            status = libhmmer.modelconfig.p7_ProfileConfig(hm, bg, gm, L, mode)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_ProfileConfig")

    cpdef Profile copy(self):
        """Return a copy of the profile with the exact same configuration.
        """
        assert self._gm != NULL

        cdef int status
        cdef Profile new = Profile.__new__(Profile)
        new.alphabet = self.alphabet

        # allocate a new profile
        with nogil:
            new._gm = libhmmer.p7_profile.p7_profile_Create(self._gm.allocM, self.alphabet._abc)
        if not new._gm:
            raise AllocationError("P7_PROFILE", sizeof(P7_PROFILE))

        # copy the current profile to the new profile
        with nogil:
            status = libhmmer.p7_profile.p7_profile_Copy(self._gm, new._gm)
        if status == libeasel.eslOK:
            return new
        else:
            raise UnexpectedError(status, "p7_profile_Copy")

    cpdef OptimizedProfile to_optimized(self):
        """Convert the profile to a platform-specific optimized profile.

        Returns:
            `OptimizedProfile`: The platform-specific optimized profile
            built using the configuration of this profile.

        """
        assert self._gm != NULL
        cdef OptimizedProfile opt = OptimizedProfile(self._gm.M, self.alphabet)
        opt.convert(self)
        return opt


cdef class ScoreData:
    """A compact representation of substitution scores and maximal extensions.

    .. versionadded:: 0.4.9

    """

    def __cinit__(self):
        self._sd = NULL

    def __init__(self, OptimizedProfile om, Profile gm = None):
        """__init__(self, om, gm=None)\n--\n
        """
        cdef P7_PROFILE*  _gm = NULL if gm is None else gm._gm
        cdef P7_OPROFILE* _om = om._om

        self.Kp = gm.alphabet.Kp
        with nogil:
            self._sd = libhmmer.p7_scoredata.p7_hmm_ScoreDataCreate(_om, _gm)
        if self._sd == NULL:
            raise AllocationError("P7_SCOREDATA", sizeof(P7_SCOREDATA))

    def __dealloc__(self):
        libhmmer.p7_scoredata.p7_hmm_ScoreDataDestroy(self._sd)
        self._sd = NULL

    def __copy__(self):
        return self.copy()

    cpdef ScoreData copy(self):
        """Create a copy of this score data object.
        """
        assert self._sd != NULL
        cdef ScoreData new = ScoreData.__new__(ScoreData)
        new.Kp = self.Kp
        with nogil:
            new._sd = libhmmer.p7_scoredata.p7_hmm_ScoreDataClone(self._sd, self.Kp)
        if new._sd == NULL:
            raise AllocationError("P7_SCOREDATA", sizeof(P7_SCOREDATA))
        return new


cdef class TopHits:
    """An immutable ranked list of top-scoring hits.

    `TopHits` are thresholded using the parameters from the pipeline, and are
    sorted by key when you obtain them from a `Pipeline` instance::

        >>> abc = thioesterase.alphabet
        >>> hits = Pipeline(abc).search_hmm(thioesterase, proteins)
        >>> hits.is_sorted(by="key")
        True

    Use `len` to query the number of top hits, and the usual indexing notation
    to extract a particular `Hit`::

        >>> len(hits)
        1
        >>> hits[0].name
        b'938293.PRJEB85.HG003687_113'

    .. versionadded:: 0.6.1
       `pickle` protocol support.

    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._th = NULL
        self._query = None
        memset(&self._pli, 0, sizeof(P7_PIPELINE))

    def __init__(self, object query not None):
        """__init__(self, query)\n--\n

        Create an empty `TopHits` instance.

        """
        self._query = query
        with nogil:
            # free allocated memory (in case __init__ is called more than once)
            libhmmer.p7_tophits.p7_tophits_Destroy(self._th)
            # allocate top hits
            self._th = libhmmer.p7_tophits.p7_tophits_Create()
            if self._th == NULL:
                raise AllocationError("P7_TOPHITS", sizeof(P7_TOPHITS))
            # clear pipeline configuration
            memset(&self._pli, 0, sizeof(P7_PIPELINE))

    def __dealloc__(self):
        libhmmer.p7_tophits.p7_tophits_Destroy(self._th)

    def __bool__(self):
        return len(self) > 0

    def __len__(self):
        assert self._th != NULL
        return self._th.N

    def __getitem__(self, index):
        if not (self._th.is_sorted_by_sortkey or self._th.is_sorted_by_seqidx):
            for i in range(self._th.N): self._th.hit[i] = &self._th.unsrt[i]
        if index < 0:
            index += self._th.N
        if index >= self._th.N or index < 0:
            raise IndexError("list index out of range")
        return Hit(self, index)

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        if id(self) not in memo:
            memo[id(self)] = self.copy()
        return memo[id(self)]

    def __add__(TopHits self, TopHits other):
        return self.merge(other)

    def __reduce__(self):
        return TopHits, (self.query,), self.__getstate__()

    def __getstate__(self):
        assert self._th != NULL

        cdef size_t    i
        cdef ptrdiff_t offset
        cdef Hit       hit    = None
        cdef list      unsrt  = []
        cdef list      hits   = []

        if self._th.N > 0:
            hit = Hit.__new__(Hit, self, 0)
            hit.hits = self

        for i in range(self._th.N):
            hit._hit = &self._th.unsrt[i]
            state = hit.__getstate__()
            unsrt.append(state)

        for i in range(self._th.N):
            offset = (<ptrdiff_t> self._th.hit[i] - <ptrdiff_t> &self._th.unsrt[0]) // sizeof(P7_HIT)
            hits.append(offset)

        return {
            "unsrt": unsrt,
            "hit": hits,
            "Nalloc": self._th.Nalloc,
            "N": self._th.N,
            "nreported": self._th.nreported,
            "nincluded": self._th.nincluded,
            "is_sorted_by_sortkey": self._th.is_sorted_by_sortkey,
            "is_sorted_by_seqidx": self._th.is_sorted_by_seqidx,
            "pipeline": {
                "by_E": self._pli.by_E,
                "E": self._pli.E,
                "T": self._pli.T,
                "dom_by_E": self._pli.dom_by_E,
                "domE": self._pli.domE,
                "domT": self._pli.domT,
                "use_bit_cutoffs": self._pli.use_bit_cutoffs,
                "inc_by_E": self._pli.inc_by_E,
                "incE": self._pli.incE,
                "incT": self._pli.incT,
                "incdom_by_E": self._pli.incdom_by_E,
                "incdomE": self._pli.incdomE,
                "incdomT": self._pli.incdomT,
                "Z": self._pli.Z,
                "domZ": self._pli.domZ,
                "Z_setby": self._pli.Z_setby,
                "domZ_setby": self._pli.domZ_setby,
                "do_max": self._pli.do_max,
                "F1": self._pli.F1,
                "F2": self._pli.F2,
                "F3": self._pli.F3,
                "B1": self._pli.B1,
                "B2": self._pli.B2,
                "B3": self._pli.B3,
                "do_biasfilter": self._pli.do_biasfilter,
                "do_null2": self._pli.do_null2,
                "nmodels": self._pli.nmodels,
                "nseqs": self._pli.nseqs,
                "nres": self._pli.nres,
                "nnodes": self._pli.nnodes,
                "n_past_msv": self._pli.n_past_msv,
                "n_past_bias": self._pli.n_past_bias,
                "n_past_vit": self._pli.n_past_vit,
                "n_past_fwd": self._pli.n_past_fwd,
                "n_output": self._pli.n_output,
                "pos_past_msv": self._pli.pos_past_msv,
                "pos_past_bias": self._pli.pos_past_bias,
                "pos_past_vit": self._pli.pos_past_vit,
                "pos_past_fwd": self._pli.pos_past_fwd,
                "pos_output": self._pli.pos_output,
                "mode": self._pli.mode,
                "long_targets": self._pli.long_targets,
                "strands": self._pli.strands,
                "W": self._pli.W,
                "block_length": self._pli.block_length
            }
        }

    def __setstate__(self, dict state):
        cdef int      status
        cdef size_t   i
        cdef uint32_t n
        cdef size_t   offset
        cdef VectorU8 hit_state

        # deallocate current data if needed
        if self._th != NULL:
            libhmmer.p7_tophits.p7_tophits_Destroy(self._th)
            self._th = NULL

        # allocate a new `P7_TOPHITS` but using exact-sized buffers from the state
        self._th = <P7_TOPHITS*> malloc(sizeof(P7_TOPHITS))
        if self._th == NULL:
            raise AllocationError("P7_TOPHITS", sizeof(P7_TOPHITS))

        # copy numbers
        self._th.N = self._th.Nalloc = state["N"]
        self._th.nreported = state["nreported"]
        self._th.nincluded = state["nincluded"]
        self._th.is_sorted_by_seqidx = state["is_sorted_by_seqidx"]
        self._th.is_sorted_by_sortkey = state["is_sorted_by_sortkey"]

        # allocate memory for hits
        self._th.unsrt = <P7_HIT*> calloc(state["N"], sizeof(P7_HIT))
        if self._th.unsrt == NULL:
            raise AllocationError("P7_HIT", sizeof(P7_HIT), state["N"])
        self._th.hit = <P7_HIT**> calloc(state["N"], sizeof(P7_HIT*))
        if self._th.hit == NULL:
            raise AllocationError("P7_HIT*", sizeof(P7_HIT*), state["N"])

        # setup sorted array
        assert len(state["hit"]) == <Py_ssize_t> self._th.N
        for i, offset in enumerate(state["hit"]):
            self._th.hit[i] = &self._th.unsrt[offset]

        # deserialize hits
        assert len(state["unsrt"]) == <Py_ssize_t> self._th.N
        for i, hit_state in enumerate(state["unsrt"]):
            n = 0
            status = libhmmer.p7_hit.p7_hit_Deserialize(
                <const uint8_t*> hit_state._data,
                &n,
                &self._th.unsrt[i],
            )
            if status != libeasel.eslOK:
                raise UnexpectedError(status, "p7_hit_Deserialize")

        # copy pipeline configuration
        self._pli.by_E = state["pipeline"]["by_E"]
        self._pli.E = state["pipeline"]["E"]
        self._pli.T = state["pipeline"]["T"]
        self._pli.dom_by_E = state["pipeline"]["dom_by_E"]
        self._pli.domE = state["pipeline"]["domE"]
        self._pli.domT = state["pipeline"]["domT"]
        self._pli.use_bit_cutoffs = state["pipeline"]["use_bit_cutoffs"]
        self._pli.inc_by_E = state["pipeline"]["inc_by_E"]
        self._pli.incE = state["pipeline"]["incE"]
        self._pli.incT = state["pipeline"]["incT"]
        self._pli.incdom_by_E = state["pipeline"]["incdom_by_E"]
        self._pli.incdomE = state["pipeline"]["incdomE"]
        self._pli.incdomT = state["pipeline"]["incdomT"]
        self._pli.Z = state["pipeline"]["Z"]
        self._pli.domZ = state["pipeline"]["domZ"]
        self._pli.Z_setby = state["pipeline"]["Z_setby"]
        self._pli.domZ_setby = state["pipeline"]["domZ_setby"]
        self._pli.do_max = state["pipeline"]["do_max"]
        self._pli.F1 = state["pipeline"]["F1"]
        self._pli.F2 = state["pipeline"]["F2"]
        self._pli.F3 = state["pipeline"]["F3"]
        self._pli.B1 = state["pipeline"]["B1"]
        self._pli.B2 = state["pipeline"]["B2"]
        self._pli.B3 = state["pipeline"]["B3"]
        self._pli.do_biasfilter = state["pipeline"]["do_biasfilter"]
        self._pli.do_null2 = state["pipeline"]["do_null2"]
        self._pli.nmodels = state["pipeline"]["nmodels"]
        self._pli.nseqs = state["pipeline"]["nseqs"]
        self._pli.nres = state["pipeline"]["nres"]
        self._pli.nnodes = state["pipeline"]["nnodes"]
        self._pli.n_past_msv = state["pipeline"]["n_past_msv"]
        self._pli.n_past_bias = state["pipeline"]["n_past_bias"]
        self._pli.n_past_vit = state["pipeline"]["n_past_vit"]
        self._pli.n_past_fwd = state["pipeline"]["n_past_fwd"]
        self._pli.n_output = state["pipeline"]["n_output"]
        self._pli.pos_past_msv = state["pipeline"]["pos_past_msv"]
        self._pli.pos_past_bias = state["pipeline"]["pos_past_bias"]
        self._pli.pos_past_vit = state["pipeline"]["pos_past_vit"]
        self._pli.pos_past_fwd = state["pipeline"]["pos_past_fwd"]
        self._pli.pos_output = state["pipeline"]["pos_output"]
        self._pli.mode = state["pipeline"]["mode"]
        self._pli.long_targets = state["pipeline"]["long_targets"]
        self._pli.strands = state["pipeline"]["strands"]
        self._pli.W = state["pipeline"]["W"]
        self._pli.block_length = state["pipeline"]["block_length"]

    # --- Properties ---------------------------------------------------------

    @property
    def mode(self):
        """`str`: Whether the hits were obtained in ``scan`` or ``search`` mode.

        .. versionadded:: 0.9.0

        """
        return "search" if self._pli.mode == p7_pipemodes_e.p7_SEARCH_SEQS else "scan"

    @property
    def query(self):
        """`object`: The query object these hits were obtained for.

        The actual type of `TopHits.query` depends on the query that was given
        to the `Pipeline`, or the `~pyhmmer.hmmer` function, that created the
        object::

            >>> hits = next(pyhmmer.hmmsearch(thioesterase, proteins))
            >>> hits.query is thioesterase
            True

        Hint:
            This property replaces the ``query_name``, ``query_accession``
            and ``query_length`` properties that were deprecated in 
            *v0.10.15* and removed in *v0.11.0*.

        .. versionadded 0.10.15

        """
        return self._query

    @property
    def Z(self):
        """`float`: The effective number of targets searched.
        """
        return self._pli.Z

    @property
    def domZ(self):
        """`float`: The effective number of significant targets searched.
        """
        return self._pli.domZ

    @property
    def E(self):
        """`float`: The per-target E-value threshold for reporting a hit.

        .. versionadded:: 0.5.0

        """
        return self._pli.E

    @property
    def T(self):
        """`float` or `None`: The per-target score threshold for reporting a hit.

        .. versionadded:: 0.5.0

        """
        return None if self._pli.by_E else self._pli.T

    @property
    def domE(self):
        """`float`: The per-domain E-value threshold for reporting a hit.

        .. versionadded:: 0.5.0

        """
        return self._pli.domE

    @property
    def domT(self):
        """`float` or `None`: The per-domain score threshold for reporting a hit.

        .. versionadded:: 0.5.0

        """
        return None if self._pli.dom_by_E else self._pli.domT

    @property
    def incE(self):
        """`float`: The per-target E-value threshold for including a hit.

        .. versionadded:: 0.5.0

        """
        return self._pli.incE

    @property
    def incT(self):
        """`float` or `None`: The per-target score threshold for including a hit.

        .. versionadded:: 0.4.8

        """
        return None if self._pli.inc_by_E else self._pli.incT

    @property
    def incdomE(self):
        """`float`: The per-domain E-value threshold for including a hit.

        .. versionadded:: 0.5.0

        """
        return self._pli.incdomE

    @property
    def incdomT(self):
        """`float` or `None`: The per-domain score threshold for including a hit.

        .. versionadded:: 0.5.0

        """
        return None if self._pli.incdom_by_E else self._pli.incdomT

    @property
    def bit_cutoffs(self):
        """`str` or `None`: The model-specific thresholding option, if any.

        .. versionadded:: 0.5.0

        """
        return next(
            (k for k,v in PIPELINE_BIT_CUTOFFS.items() if v == self._pli.use_bit_cutoffs),
            None
        )

    @property
    def searched_models(self):
        """`int`: The number of models searched.

        .. versionadded:: 0.5.0

        """
        return self._pli.nmodels

    @property
    def searched_nodes(self):
        """`int`: The number of model nodes searched.

        .. versionadded:: 0.5.0

        """
        return self._pli.nnodes

    @property
    def searched_sequences(self):
        """`int`: The number of sequences searched.

        .. versionadded:: 0.5.0

        """
        return self._pli.nseqs

    @property
    def searched_residues(self):
        """`int`: The number of residues searched.

        .. versionadded:: 0.5.0

        """
        return self._pli.nres

    @property
    def long_targets(self):
        """`bool`: Whether these hits were produced by a long targets pipeline.

        .. versionadded:: 0.5.0

        """
        return self._pli.long_targets

    @property
    def strand(self):
        """`str` or `None`: The strand these hits were obtained from.

        Is always `None` when the hits were not obtained from a long targets
        pipeline, or when the long targets pipeline was configured to search
        both strands.

        .. versionadded:: 0.5.0

        """
        if self._pli.long_targets:
            if self._pli.strands == p7_strands_e.p7_STRAND_TOPONLY:
                return "watson"
            elif self._pli.strands == p7_strands_e.p7_STRAND_BOTTOMONLY:
                return "crick"
        return None

    @property
    def block_length(self):
        """`int` or `None`: The block length these hits were obtained with.

        Is always `None` when the hits were not obtained from a long targets
        pipeline.

        .. versionadded:: 0.5.0

        """
        return self._pli.block_length if self._pli.long_targets else None

    @property
    def included(self):
        """iterator of `Hit`: An iterator over the hits marked as *included*.

        .. versionadded:: 0.7.0

        """
        return SizedIterator(
            self._th.nincluded,
            filter(operator.attrgetter("included"), self)
        )

    @property
    def reported(self):
        """iterator of `Hit`: An iterator over the hits marked as *reported*.

        .. versionadded:: 0.7.0

        """
        return SizedIterator(
            self._th.nreported,
            filter(operator.attrgetter("reported"), self)
        )

    # --- Utils --------------------------------------------------------------

    cdef int _threshold(self, Pipeline pipeline) except 1 nogil:
        # threshold the top hits with the given pipeline numbers
        cdef int status = libhmmer.p7_tophits.p7_tophits_Threshold(self._th, pipeline._pli)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_tophits_Threshold")
        # record pipeline configuration
        # NOTE: the pointers on the copy are set to NULL by precaution, but since
        # `pipeline._pli` is allocated by the Python object directly, it will not
        # be deallocated, so there should be no risk of double free nevertheless.
        memcpy(&self._pli, pipeline._pli, sizeof(P7_PIPELINE))
        self._pli.oxf = self._pli.oxb = self._pli.fwd = self._pli.bck = NULL
        self._pli.r = NULL
        self._pli.ddef = NULL
        self._pli.hfp = NULL
        return 0

    cdef int _sort_by_key(self) except 1 nogil:
        cdef int status = libhmmer.p7_tophits.p7_tophits_SortBySortkey(self._th)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_tophits_SortBySortkey")
        return 0

    cdef int _sort_by_seqidx(self) except 1 nogil:
        cdef int status = libhmmer.p7_tophits.p7_tophits_SortBySeqidxAndAlipos(self._th)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_tophits_SortBySeqidxAndAlipos")
        return 0

    cdef void _check_threshold_parameters(self, const P7_PIPELINE* other) except *:
        # check comparison counters are consistent
        if self._pli.long_targets and not other.long_targets:
            raise ValueError("Trying to merge a `TopHits` from a long targets pipeline to a `TopHits` from a regular pipeline.")
        if self._pli.Z_setby != other.Z_setby:
            raise ValueError("Trying to merge `TopHits` with `Z` values obtained with different methods.")
        elif self._pli.Z_setby != p7_zsetby_e.p7_ZSETBY_NTARGETS and self._pli.Z != other.Z:
            raise ValueError("Trying to merge `TopHits` obtained from pipelines manually configured to different `Z` values.")
        if self._pli.domZ_setby != other.domZ_setby:
            raise ValueError("Trying to merge `TopHits` with `domZ` values obtained with different methods.")
        elif self._pli.domZ_setby != p7_zsetby_e.p7_ZSETBY_NTARGETS and self._pli.domZ != other.domZ:
            raise ValueError("Trying to merge `TopHits` obtained from pipelines manually configured to different `domZ` values.")
        # check threshold modes are consistent
        if self._pli.by_E != other.by_E:
            raise ValueError(f"Trying to merge `TopHits` obtained from pipelines with different reporting threshold modes: {self._pli.by_E} != {other.by_E}")
        elif self._pli.dom_by_E != other.dom_by_E:
            raise ValueError("Trying to merge `TopHits` obtained from pipelines with different domain reporting threshold modes")
        elif self._pli.inc_by_E != other.inc_by_E:
            raise ValueError("Trying to merge `TopHits` obtained from pipelines with different inclusion threshold modes")
        elif self._pli.incdom_by_E != other.incdom_by_E:
            raise ValueError("Trying to merge `TopHits` obtained from pipelines with different domain inclusion threshold modes")
        # check inclusion and reporting threshold are the same
        if (self._pli.by_E and self._pli.E != other.E) or (not self._pli.by_E and self._pli.T != other.T):
            raise ValueError("Trying to merge `TopHits` obtained from pipelines with different reporting thresholds.")
        elif (self._pli.inc_by_E and self._pli.incE != other.incE) or (not self._pli.inc_by_E and self._pli.incT != other.incT):
            raise ValueError("Trying to merge `TopHits` obtained from pipelines with different inclusion thresholds.")
        elif (self._pli.dom_by_E and self._pli.domE != other.domE) or (not self._pli.dom_by_E and self._pli.domT != other.domT):
            raise ValueError("Trying to merge `TopHits` obtained from pipelines with different domain reporting thresholds.")
        elif (self._pli.incdom_by_E and self._pli.incdomE != other.incdomE) or (not self._pli.incdom_by_E and self._pli.incdomT != other.incdomT):
            raise ValueError("Trying to merge `TopHits` obtained from pipelines with different domain inclusion thresholds.")

    # --- Methods ------------------------------------------------------------

    cpdef TopHits copy(self):
        """Create a copy of this `TopHits` instance.

        .. versionadded:: 0.5.0

        """
        assert self._th != NULL
        assert self._th.N >= 0

        cdef TopHits copy = TopHits.__new__(TopHits)

        # record query metatada
        copy._query = self._query

        with nogil:
            # copy pipeline configuration
            memcpy(&copy._pli, &self._pli, sizeof(P7_PIPELINE))
            # copy top hits
            # WARN(@althonos): `p7_tophits_Clone` will fail when called
            #                  on an empty P7_TOPHITS as it attempts a
            #                  zero-allocation of the target arrays, it
            #                  should be reported as a bug at some point.
            if self._th.N > 0:
                copy._th = libhmmer.p7_tophits.p7_tophits_Clone(self._th)
            else:
                copy._th = libhmmer.p7_tophits.p7_tophits_Create()
        if copy._th == NULL:
            raise AllocationError("P7_TOPHITS", sizeof(P7_TOPHITS))

        return copy

    cpdef int compare_ranking(self, KeyHash ranking) except -1:
        """Compare current top hits to previous top hits ranking.

        This method is used by ``jackhmmer`` to record the hits obtained
        during each iteration, so that the inner loop can converge.

        Arguments:
            ranking (`~pyhmmer.easel.KeyHash`): A keyhash containing the
                ranks of the top hits from a previous run.

        Returns:
            `int`: The number of new hits found in this iteration.

        .. versionadded:: 0.6.0

        """
        assert self._th != NULL
        assert ranking._kh != NULL

        cdef int new_hits
        cdef int status

        with nogil:
            status = libhmmer.p7_tophits.p7_tophits_CompareRanking(
                self._th,
                ranking._kh,
                &new_hits
            )
        if status == libeasel.eslEMEM:
            raise AllocationError("ESL_KEYHASH", sizeof(ESL_KEYHASH))
        elif status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_tophits_CompareRanking")

        return new_hits

    cpdef void sort(self, str by="key") except *:
        """Sort hits in the current instance using the given method.

        Arguments:
            by (`str`): The comparison method to use to compare hits.
                Allowed values are: ``key`` (the default) to sort by key, or
                ``seqidx`` to sort by sequence index and alignment position.

        """
        assert self._th != NULL
        if by == "key":
            with nogil:
                self._sort_by_key()
        elif by == "seqidx":
            with nogil:
                self._sort_by_seqidx()
        else:
            raise InvalidParameter("by", by, choices=["key", "seqidx"])

    cpdef bint is_sorted(self, str by="key") except *:
        """Check whether or not the hits are sorted with the given method.

        See `~pyhmmer.plan7.TopHits.sort` for a list of allowed values for
        the ``by`` argument.

        """
        assert self._th != NULL
        if by == "key":
            return self._th.is_sorted_by_sortkey
        elif by == "seqidx":
            return self._th.is_sorted_by_seqidx
        else:
            raise InvalidParameter("by", by, choices=["key", "seqidx"])

    cpdef MSA to_msa(
        self,
        Alphabet alphabet,
        list sequences=None,
        list traces=None,
        bint trim=False,
        bint digitize=False,
        bint all_consensus_cols=False
    ):
        """Create multiple alignment of all included domains.

        Arguments:
            alphabet (`~pyhmmer.easel.Alphabet`): The alphabet of the
                HMM this `TopHits` was obtained from. It is required to
                convert back hits to single sequences.
            sequences (`list` of `~pyhmmer.easel.Sequence`, optional): A list
                of additional sequences to include in the alignment.
            traces (`list` of `~plan7.Trace`, optional): A list of
                additional traces to include in the alignment.

        Keyword Arguments:
            trim (`bool`): Trim off any residues that get assigned to
                flanking :math:`N` and :math:`C` states (in profile traces)
                or :math:`I_0` and :math:`I_m` (in core traces).
            digitize (`bool`): If set to `True`, returns a  `DigitalMSA`
                instead of a `TextMSA`.
            all_consensus_cols (`bool`): Force a column to be created for
                every consensus column in the model, even if it means having
                all gap character in a column.

        Returns:
            `~pyhmmer.easel.MSA`: A multiple sequence alignment containing
            the reported hits, either a `~pyhmmer.easel.TextMSA` or a
            `~pyhmmer.easel.DigitalMSA` depending on the value of the
            ``digitize`` argument.

        .. versionadded:: 0.3.0

        .. versionchanged:: 0.6.0
           Added the ``sequences`` and ``traces`` arguments.

        """
        assert self._th != NULL
        assert alphabet._abc != NULL

        cdef int        status
        cdef int        extras
        cdef int        flags  = libhmmer.p7_DEFAULT
        cdef MSA        msa
        cdef Sequence   seq
        cdef Trace      trace
        cdef ESL_SQ**   sqarr  = NULL
        cdef P7_TRACE** trarr  = NULL

        # prepare additional flags
        if trim:
            flags |= libhmmer.p7_TRIM
        if all_consensus_cols:
            flags |= libhmmer.p7_ALL_CONSENSUS_COLS
        if digitize:
            flags |= libhmmer.p7_DIGITIZE
            msa = DigitalMSA.__new__(DigitalMSA, alphabet)
        else:
            msa = TextMSA.__new__(TextMSA)

        # check if any additional sequences/traces are given
        sequences = [] if sequences is None else sequences
        traces = [] if traces is None else traces
        if len(sequences) != len(traces):
            raise ValueError("`sequences` and `traces` must have the same length")
        else:
            extras = len(sequences)

        try:
            # make arrays for the additional sequences/traces
            if extras > 0:
                sqarr = <ESL_SQ**> calloc(sizeof(ESL_SQ*), extras)
                if sqarr == NULL:
                    raise AllocationError("ESL_SQ*", sizeof(ESL_SQ*), extras)
                trarr = <P7_TRACE**> calloc(sizeof(P7_TRACE*), extras)
                if trarr == NULL:
                    raise AllocationError("P7_TRACE*", sizeof(ESL_SQ*), extras)
                for i, (seq, trace) in enumerate(zip(sequences, traces)):
                    sqarr[i] = seq._sq
                    trarr[i] = trace._tr
            # build the alignment
            status = libhmmer.p7_tophits.p7_tophits_Alignment(
                self._th,
                alphabet._abc,
                sqarr,
                trarr,
                extras,
                flags,
                &msa._msa
            )
            if status == libeasel.eslOK:
                return msa
            elif status == libeasel.eslFAIL:
                raise ValueError("No included domains found")
            else:
                raise UnexpectedError(status, "p7_tophits_Alignment")
        finally:
            free(sqarr)
            free(trarr)

    cpdef void write(self, object fh, str format="targets", bint header=True) except *:
        """Write the hits in tabular format to a file-like object.

        Arguments:
            fh (`io.IOBase`): A Python file handle, opened in binary mode.
            format (`str`): The tabular format in which to write the hits.
            header (`bool`): Whether to write a table header. Ignored
                when writing in the ``pfam`` format.

        Hint:
            The hits can be written in one of the following formats:

            ``targets``
              A tabular output format of per-target hits, as obtained
              with the ``--tblout`` output flag of ``hmmsearch`` or
              ``hmmscan``.

            ``domains``
              A tabular output format of per-domain hits, as obtained
              with the ``--domtblout`` output flag of ``hmmsearch`` or
              ``hmmscan``.

            ``pfam``
              A tabular output format suitable for Pfam, merging per-sequence
              and per-domain hits in a single file, with fewer fields and
              sorted by score.

        .. versionadded:: 0.6.1

        """
        cdef FILE* file
        cdef str   fname
        cdef int   status
        cdef bytes qname  = b"-"
        cdef bytes qacc   = b"-"

        if self._query is not None:
            if self._query.name is not None:
                qname = self._query.name
            if self._query.accession is not None:
                qacc = self._query.accession

        file = fopen_obj(fh, "w")
        try:
            if format == "targets":
                fname = "p7_tophits_TabularTargets"
                status = libhmmer.p7_tophits.p7_tophits_TabularTargets(
                    file,
                    qname,
                    qacc,
                    self._th,
                    &self._pli,
                    header
                )
            elif format == "domains":
                fname = "p7_tophits_TabularDomains"
                status = libhmmer.p7_tophits.p7_tophits_TabularDomains(
                    file,
                    qname,
                    qacc,
                    self._th,
                    &self._pli,
                    header
                )
            elif format == "pfam":
                fname = "p7_tophits_TabularXfam"
                status = libhmmer.p7_tophits.p7_tophits_TabularXfam(
                    file,
                    qname,
                    qacc,
                    self._th,
                    &self._pli,
                )
            else:
                raise InvalidParameter("format", format, choices=["targets", "domains", "pfam"])
            if status != libeasel.eslOK:
                _reraise_error()
                raise UnexpectedError(status, fname)
        finally:
            fclose(file)

    def merge(self, *others):
        """Concatenate the hits from this instance and ``others``.

        If the ``Z`` and ``domZ`` values used to compute E-values were
        computed by the `Pipeline` from the number of targets, the returned
        object will update them by summing ``self.Z`` and ``other.Z``. If
        they were set manually, the manual value will be kept, provided
        both values are equal.

        Returns:
            `~pyhmmer.plan7.TopHits`: A new collection of hits containing
            a copy of all the hits from ``self`` and ``other``, sorted
            by key.

        Raises:
            `ValueError`: When trying to merge together several hits
                obtained from different `Pipeline` with incompatible
                parameters.

        Caution:
            This should only be done for hits obtained for the same domain
            on similarly configured pipelines. Some internal checks will be
            done to ensure this is not the case, but the results may not be
            consistent at all.

        Example:
            >>> pli = Pipeline(thioesterase.alphabet)
            >>> hits1 = pli.search_hmm(thioesterase, proteins[:1000])
            >>> hits2 = pli.search_hmm(thioesterase, proteins[1000:2000])
            >>> hits3 = pli.search_hmm(thioesterase, proteins[2000:])
            >>> merged = hits1.merge(hits2, hits3)

        .. versionadded:: 0.5.0

        """
        assert self._th != NULL

        cdef TopHits other
        cdef TopHits other_copy
        cdef TopHits merged     = self.copy()
        cdef int     status     = libeasel.eslOK

        for other in others:
            assert other._th != NULL

            # copy hits (`p7_tophits_Merge` effectively destroys the old storage
            # but because of Python references we cannot be sure that the data is
            # not referenced anywhere else)
            other_copy = other.copy()

            # check that names/accessions are consistent
            if merged._query != other._query:
                raise ValueError("Trying to merge `TopHits` obtained from different queries")

            # just store the copy if merging inside an empty uninitialized `TopHits`
            if merged._th.N == 0:
                merged._query = other._query
                memcpy(&merged._pli, &other_copy._pli, sizeof(P7_PIPELINE))
                merged._th, other_copy._th = other_copy._th, merged._th
                continue

            # check that the parameters are the same
            merged._check_threshold_parameters(&other._pli)

            # merge everything
            with nogil:
                # merge the top hits
                status = libhmmer.p7_tophits.p7_tophits_Merge(merged._th, other_copy._th)
                if status != libeasel.eslOK:
                    raise UnexpectedError(status, "p7_tophits_Merge")
                # merge the pipelines
                status = libhmmer.p7_pipeline.p7_pipeline_Merge(&merged._pli, &other_copy._pli)
                if status != libeasel.eslOK:
                    raise UnexpectedError(status, "p7_pipeline_Merge")

        # Reset nincluded/nreports before thresholding
        # TODO(@althonos, @zdk123): Replace with `p7_tophits_Threshold` as implemented
        #                  in EddyRivasLab/hmmer#307 when formally released.
        for i in range(merged._th.N):
            merged._th.hit[i].flags &= (~p7_hitflags_e.p7_IS_REPORTED)
            merged._th.hit[i].flags &= (~p7_hitflags_e.p7_IS_INCLUDED)
            merged._th.hit[i].nincluded = 0
            merged._th.hit[i].nreported = 0
            for j in range(merged._th.hit[i].ndom):
                merged._th.hit[i].dcl[j].is_reported = False
                merged._th.hit[i].dcl[j].is_included = False

        # threshold the merged hits with new values
        status = libhmmer.p7_tophits.p7_tophits_Threshold(merged._th, &merged._pli)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_tophits_Threshold")

        # return the merged hits
        return merged


@cython.freelist(8)
@cython.no_gc_clear
cdef class Trace:
    """A traceback for the alignment of a model to a sequence.

    .. versionadded:: 0.4.7

    """

    @classmethod
    def from_sequence(cls, Sequence sequence not None):
        """Create a faux trace from a single sequence.

        .. versionadded:: 0.6.0

        """
        cdef int   i
        cdef int   status
        cdef Trace trace  = cls()

        status = libhmmer.p7_trace.p7_trace_Append(trace._tr, p7t_statetype_e.p7T_B, 0, 0)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_trace_Append")

        for i in range(1, sequence._sq.n + 1):
            status = libhmmer.p7_trace.p7_trace_Append(trace._tr, p7t_statetype_e.p7T_M, i, i)
            if status != libeasel.eslOK:
                raise UnexpectedError(status, "p7_trace_Append")

        status = libhmmer.p7_trace.p7_trace_Append(trace._tr, p7t_statetype_e.p7T_E, 0, 0)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_trace_Append")

        trace._tr.M = trace._tr.L = sequence._sq.n
        return trace

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._tr = NULL
        self.traces = None

    def __init__(self, posteriors=False):
        """__init__(self, posteriors=False)\n--\n

        Create a new `Trace` object.

        Arguments:
            posteriors (`bool`): Whether or not to allocate additional memory
                for the storage of posterior probabilties.

        """
        # make `__init__` calllable more than once to avoid bugs
        if self._tr != NULL:
            if self.traces is None:
                libhmmer.p7_trace.p7_trace_Destroy(self._tr)
            self._tr = NULL
            self.traces = None
        if posteriors:
            self._tr = libhmmer.p7_trace.p7_trace_CreateWithPP()
        else:
            self._tr = libhmmer.p7_trace.p7_trace_Create()
        if self._tr == NULL:
            raise AllocationError("P7_TRACE", sizeof(P7_TRACE))

    def __dealloc__(self):
        if self.traces is None:
            libhmmer.p7_trace.p7_trace_Destroy(self._tr)

    def __eq__(self, object other):
        assert self._tr != NULL

        if not isinstance(other, Trace):
            return NotImplemented

        cdef Trace t      = <Trace> other
        cdef int   status = libhmmer.p7_trace.p7_trace_Compare(self._tr, t._tr, 0.0)

        if status == libeasel.eslOK:
            return True
        elif status == libeasel.eslFAIL:
            return False
        else:
            raise UnexpectedError(status, "p7_trace_Compare")

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"<{ty} M={self.M!r} L={self.L!r}>"

    # --- Properties ---------------------------------------------------------

    @property
    def M(self):
        """`int`: The model length.
        """
        assert self._tr != NULL
        return self._tr.M

    @property
    def L(self):
        """`int`: The sequence length.
        """
        assert self._tr != NULL
        return self._tr.L

    @property
    def posterior_probabilities(self):
        """`VectorF` or `None`: The posterior probability of each residue.
        """
        assert self._tr != NULL
        if self._tr.pp == NULL:
            return None
        cdef VectorF pp = VectorF.__new__(VectorF)
        pp._owner = self
        pp._n = pp._shape[0] = self._tr.N
        pp._data = NULL if self._tr.N == 0 else <void*> self._tr.pp
        return pp

    # --- Methods ------------------------------------------------------------

    cpdef float expected_accuracy(self):
        """Returns the sum of the posterior residue decoding probabilities.
        """
        assert self._tr != NULL
        with nogil:
            return libhmmer.p7_trace.p7_trace_GetExpectedAccuracy(self._tr)

    cpdef float score(self, DigitalSequence sequence, Profile profile) except *:
        """Score traceback for target sequence using given profile.

        .. versionadded:: 0.10.5

        """
        assert self._tr != NULL
        assert sequence._sq != NULL
        assert profile._gm != NULL

        cdef int   status
        cdef float score

        with nogil:
            status = libhmmer.p7_trace.p7_trace_Score(self._tr, sequence._sq.dsq, profile._gm, &score)
        if status == libeasel.eslEINVAL:
            error = _recover_error()
            raise ValueError("Invalid trace") from error
        elif status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_trace_Score")

        return score


cdef class Traces:
    """A list of tracebacks obtained by aligning several sequences to a model.

    .. versionadded:: 0.4.7

    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._traces = NULL
        self._ntraces = 0

    def __dealloc__(self):
        libhmmer.p7_trace.p7_trace_DestroyArray(self._traces, self._ntraces)

    def __init__(self):
        """__init__(self)\n--\n

        Create an empty list of traces.

        """
        pass

    def __len__(self):
        return self._ntraces

    def __getitem__(self, ssize_t idx):
        assert self._traces != NULL

        cdef Trace trace

        if idx < 0:
            idx += self._ntraces
        if idx >= <ssize_t> self._ntraces or idx < 0:
            raise IndexError("list index out of range")

        trace = Trace.__new__(Trace)
        trace.traces = self
        trace._tr = self._traces[idx]

        return trace

    def __eq__(self, object other):
        assert self._traces != NULL

        cdef size_t i
        cdef int    status
        cdef Traces other_tr

        if not isinstance(other, Traces):
            return NotImplemented
        other_tr = <Traces> other

        if self._ntraces != other_tr._ntraces:
            return False
        for i in range(self._ntraces):
            status = libhmmer.p7_trace.p7_trace_Compare(self._traces[i], other_tr._traces[i], 0.0)
            if status == libeasel.eslFAIL:
                return False
            elif status != libeasel.eslOK:
                raise UnexpectedError(status, "p7_trace_Compare")

        return True


cdef class TraceAligner:
    """A factory for aligning several sequences to a reference model.

    Example:
        >>> aligner = TraceAligner()
        >>> traces = aligner.compute_traces(thioesterase, proteins[:100])
        >>> msa = aligner.align_traces(thioesterase, proteins[:100], traces)

    .. versionadded:: 0.4.7

    """

    # --- Magic methods ------------------------------------------------------

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"{ty}()"

    # --- Methods ------------------------------------------------------------

    cpdef Traces compute_traces(self, HMM hmm, DigitalSequenceBlock sequences):
        """Compute traces for a collection of sequences relative to an HMM.

        Arguments:
            hmm (`~pyhmmer.plan7.HMM`): The reference HMM to use for the
                alignment.
            sequences (`~pyhmmer.easel.DigitalSequenceBlock`): The target
                sequences to align to the HMM.

        Returns:
            `~pyhmmer.plan7.Traces`: The traces for each sequence.

        Raises:
            `~pyhmmer.errors.AlphabetMismatch`: when the alphabet of any
                of the sequences does not correspond to the HMM alphabet.
            `~pyhmmer.errors.InvalidHMM`: when given a HMM that is not
                valid.

        .. versionchanged:: 0.7.0
            Targets must now be inside a `~pyhmmer.easel.DigitalSequenceBlock`.

        """
        cdef int                 status
        cdef ssize_t             i
        cdef Trace               trace
        cdef Traces              traces = Traces()
        cdef ssize_t             nseq   = len(sequences)
        cdef char[eslERRBUFSIZE] errbuf

        if nseq < 0:
            raise ValueError("Cannot compute traces for a negative number of sequences")
        elif nseq == 0:
            return traces

        # check alphabet
        if not hmm.alphabet._eq(sequences.alphabet):
            raise AlphabetMismatch(hmm.alphabet, sequences.alphabet)
        # check HMM validity, otherwise the function may segfault
        hmm.validate(tolerance=1e-3)

        # allocate the return array of traces and create empty traces
        traces._ntraces = nseq
        traces._traces = <P7_TRACE**> calloc(nseq, sizeof(P7_TRACE*))
        if traces._traces == NULL:
            raise AllocationError("P7_TRACE**", sizeof(P7_TRACE*), nseq)
        for i in range(nseq):
            traces._traces[i] = libhmmer.p7_trace.p7_trace_CreateWithPP()
            if traces._traces[i] == NULL:
                raise AllocationError("P7_TRACE", sizeof(P7_TRACE))

        # compute the traces
        with nogil:
            status = libhmmer.tracealign.p7_tracealign_computeTraces(
                hmm._hmm,
                sequences._refs,
                0,
                nseq,
                traces._traces
            )
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_tracealign_computeTraces")
        return traces

    cpdef MSA align_traces(
        self,
        HMM hmm,
        DigitalSequenceBlock sequences,
        Traces traces,
        bint digitize=False,
        bint trim=False,
        bint all_consensus_cols=False
    ):
        """Compute traces for a collection of sequences relative to an HMM.

        Arguments:
            hmm (`~pyhmmer.plan7.HMM`): The reference HMM to use for the
                alignment.
            sequences (`~pyhmmer.easel.DigitalSequenceBlock`): The target
                sequences to align to the HMM.
            traces (`~pyhmmer.plan7.Traces`): The traces corresponding to the
                alignment of ``sequences`` to ``hmm``, obtained by a previous
                call to `~pyhmmer.plan7.TraceAligner.compute_traces`.
            digitize (`bool`): If set to `True`, returns a `DigitalMSA`
                instead of a `TextMSA`.
            trim (`bool`): Trim off any residues that get assigned to
                flanking :math:`N` and :math:`C` states (in profile traces)
                or :math:`I_0` and :math:`I_m` (in core traces).
            all_consensus_cols (`bool`): Force a column to be created for
                every consensus column in the model, even if it means having
                all gap character in a column. *Note that this is enabled by
                default for* ``hmmalign`` *since HMMER v3.1, but disabled
                here.*

        Returns:
            `~pyhmmer.easel.MSA`: A multiple sequence alignment containing
            the aligned sequences, either a `~pyhmmer.easel.TextMSA` or a
            `~pyhmmer.easel.DigitalMSA` depending on the value of the
            ``digitize`` argument.

        Raises:
            `~pyhmmer.errors.AlphabetMismatch`: when the alphabet of any
                of the sequences does not correspond to the HMM alphabet.
            `~pyhmmer.errors.InvalidHMM`: when given a HMM that is not
                valid.

        .. versionchanged:: 0.7.0
            Targets must now be inside a `~pyhmmer.easel.DigitalSequenceBlock`.

        """
        cdef size_t              i
        cdef MSA                 msa
        cdef int                 status
        cdef char[eslERRBUFSIZE] errbuf
        cdef ssize_t             nseq   = len(sequences)
        cdef ssize_t             ntr    = len(traces)
        cdef int                 flags  = 0

        # check optional flags and prepare the returned MSA
        if trim:
            flags |= libhmmer.p7_TRIM
        if all_consensus_cols:
            flags |= libhmmer.p7_ALL_CONSENSUS_COLS
        if digitize:
            flags |= libhmmer.p7_DIGITIZE
            msa = DigitalMSA.__new__(DigitalMSA, hmm.alphabet)
        else:
            msa = TextMSA.__new__(TextMSA)

        # check number of sequences
        if nseq < 0:
            raise ValueError("Cannot align traces for a negative number of sequences")
        elif nseq != ntr:
            raise ValueError(f"Sequences and traces lengths mismatch ({nseq} sequences, {ntr} traces)")
        elif nseq == 0:
            return msa
        # check alphabet
        if not hmm.alphabet._eq(sequences.alphabet):
            raise AlphabetMismatch(hmm.alphabet, sequences.alphabet)
        # check HMM validity, otherwise the function may segfault
        hmm.validate(tolerance=1e-3)

        # make the alignments
        with nogil:
            status = libhmmer.tracealign.p7_tracealign_Seqs(
                sequences._refs,
                traces._traces,
                nseq,
                hmm._hmm.M,
                flags,
                hmm._hmm,
                &msa._msa
            )
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_tracealign_Seqs")
        return msa


class Transitions(enum.IntEnum):
    """A helper enum for indices of the HMM transition probability matrix.

    The Plan 7 model architecture used in HMMER describes a HMM which has
    3 states and 7 transitions (hence the name) for every node of the model.
    The HMM can transition from a *match* state (:math:`M_n`) to the next
    match stage (:math:`M_{n+1}`), to an *insertion* state (:math:`I_n`)
    or a *deletion* state (:math:`D_{n+1}`). However, there are no
    transitions between a *deletion* and *insertion* state.

    .. versionadded:: 0.8.1

    """
    MM = p7h_transitions_e.p7H_MM
    MI = p7h_transitions_e.p7H_MI
    MD = p7h_transitions_e.p7H_MD
    IM = p7h_transitions_e.p7H_IM
    II = p7h_transitions_e.p7H_II
    DM = p7h_transitions_e.p7H_DM
    DD = p7h_transitions_e.p7H_DD


# --- Module init code -------------------------------------------------------

impl_Init()
p7_FLogsumInit()
