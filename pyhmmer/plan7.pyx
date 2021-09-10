# coding: utf-8
# cython: language_level=3, linetrace=True
"""High-level interface to the Plan7 data model.

Plan7 is the model architecture used by HMMER since HMMER2.

See Also:
    Details about the Plan 7 architecture in the `HMMER documentation
    <http://www.csb.yale.edu/userguides/seq/hmmer/docs/node11.html>`_.

"""

# --- C imports --------------------------------------------------------------

cimport cython
from cpython.bytes cimport PyBytes_FromStringAndSize
from cpython.ref cimport PyObject
from cpython.exc cimport PyErr_Clear
from cpython.unicode cimport PyUnicode_DecodeASCII
from libc.math cimport exp, ceil
from libc.stdlib cimport calloc, malloc, realloc, free
from libc.stdint cimport uint8_t, uint32_t
from libc.stdio cimport fprintf, FILE, stdout, fclose
from libc.string cimport strdup, strndup, strlen, strncpy
from libc.time cimport ctime, strftime, time, time_t, tm, localtime_r
from unicode cimport PyUnicode_DATA, PyUnicode_KIND, PyUnicode_READ, PyUnicode_READY, PyUnicode_GET_LENGTH

cimport libeasel
cimport libeasel.sq
cimport libeasel.alphabet
cimport libeasel.fileparser
cimport libeasel.random
cimport libeasel.getopts
cimport libeasel.vec
cimport libhmmer.modelconfig
cimport libhmmer.p7_hmm
cimport libhmmer.p7_builder
cimport libhmmer.p7_bg
cimport libhmmer.p7_domaindef
cimport libhmmer.p7_hmmfile
cimport libhmmer.p7_pipeline
cimport libhmmer.p7_prior
cimport libhmmer.p7_profile
cimport libhmmer.p7_tophits
from libeasel cimport eslERRBUFSIZE, eslCONST_LOG2R
from libeasel.alphabet cimport ESL_ALPHABET, esl_alphabet_Create, esl_abc_ValidateType
from libeasel.getopts cimport ESL_GETOPTS, ESL_OPTIONS
from libeasel.sq cimport ESL_SQ
from libhmmer cimport p7_LOCAL, p7_EVPARAM_UNSET, p7_CUTOFF_UNSET, p7_NEVPARAM, p7_NCUTOFFS, p7_offsets_e, p7_cutoffs_e, p7_evparams_e
from libhmmer.logsum cimport p7_FLogsumInit
from libhmmer.p7_builder cimport P7_BUILDER, p7_archchoice_e, p7_wgtchoice_e, p7_effnchoice_e
from libhmmer.p7_hmm cimport p7H_NTRANSITIONS, p7H_TC, p7H_GA, p7H_NC, p7H_MAP
from libhmmer.p7_hmmfile cimport p7_hmmfile_formats_e
from libhmmer.p7_tophits cimport p7_hitflags_e
from libhmmer.p7_alidisplay cimport P7_ALIDISPLAY
from libhmmer.p7_pipeline cimport P7_PIPELINE, p7_pipemodes_e, p7_zsetby_e
from libhmmer.p7_profile cimport p7_LOCAL, p7_GLOCAL, p7_UNILOCAL, p7_UNIGLOCAL

IF HMMER_IMPL == "VMX":
    from libhmmer.impl_vmx cimport p7_oprofile, p7_omx, impl_Init
    from libhmmer.impl_vmx.p7_oprofile cimport P7_OPROFILE, p7_oprofile_Dump, p7O_NQB
    from libhmmer.impl_vmx.io cimport p7_oprofile_Write
ELIF HMMER_IMPL == "SSE":
    from libhmmer.impl_sse cimport p7_oprofile, p7_omx, impl_Init, p7_SSVFilter, p7O_EXTRA_SB
    from libhmmer.impl_sse.p7_oprofile cimport P7_OPROFILE, p7_oprofile_Dump, p7O_NQB
    from libhmmer.impl_sse.io cimport p7_oprofile_Write

from .easel cimport (
    Alphabet,
    Sequence,
    DigitalSequence,
    MSA,
    TextMSA,
    DigitalMSA,
    VectorF,
    MatrixF,
    MatrixU8,
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

IF UNAME_SYSNAME == "Linux":
    include "fileobj/linux.pxi"
ELIF UNAME_SYSNAME == "Darwin" or UNAME_SYSNAME.endswith("BSD"):
    include "fileobj/bsd.pxi"

# --- Python imports ---------------------------------------------------------

import array
import collections.abc
import datetime
import errno
import math
import io
import os
import sys
import warnings

from .errors import AllocationError, UnexpectedError, AlphabetMismatch
from .utils import peekable


# --- Cython classes ---------------------------------------------------------


cdef class Alignment:
    """An alignment of a sequence to a profile.

    Attributes:
        domain (`Domain`): The domain this alignment corresponds to.

    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self, Domain domain):
        self.domain = domain
        self._ad = domain._dom.ad

    def __len__(self):
        return self.hmm_to - self.hmm_from

    # --- Properties ---------------------------------------------------------

    @property
    def hmm_from(self):
        """`int`: The start coordinate of the alignment in the query HMM.
        """
        return self._ad.hmmfrom

    @property
    def hmm_to(self):
        """`int`: The end coordinate of the alignment in the query HMM.
        """
        return self._ad.hmmto

    @property
    def hmm_name(self):
        """`bytes`: The name of the query HMM.
        """
        return <bytes> self._ad.hmmname

    @property
    def hmm_accession(self):
        """`bytes`: The accession of the query, or its name if it has none.

        .. versionadded:: 0.1.4

        """
        return <bytes> self._ad.hmmacc

    @property
    def hmm_sequence(self):
        """`str`: The sequence of the query HMM in the alignment.
        """
        return self._ad.model.decode('ascii')

    @property
    def target_from(self):
        """`int`: The start coordinate of the alignment in the target sequence.
        """
        return self._ad.sqfrom

    @property
    def target_name(self):
        """`bytes`: The name of the target sequence.
        """
        return <bytes> self._ad.sqname

    @property
    def target_sequence(self):
        """`str`: The sequence of the target sequence in the alignment.
        """
        return self._ad.aseq.decode('ascii')

    @property
    def target_to(self):
        """`int`: The end coordinate of the alignment in the target sequence.
        """
        return self._ad.sqto

    @property
    def identity_sequence(self):
        """`str`: The identity sequence between the query and the target.
        """
        return self._ad.mline.decode('ascii')


cdef class Background:
    """The null background model of HMMER.

    Attributes:
        alphabet (`~pyhmmer.easel.Alphabet`): The alphabet of the
            backgound model.
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
        """__init__(self, alphabet, uniform=False)\n--

        Create a new background model for the given ``alphabet``.

        Arguments:
          alphabet (`pyhmmer.easel.Alphabet`): The alphabet to create the
              background model with.
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
            raise AllocationError("P7_BG")
        # expose the residue frequencies as the `residue_frequencies` attribute
        self.residue_frequencies = VectorF.__new__(VectorF)
        self.residue_frequencies._data = &(self._bg.f[0])
        self.residue_frequencies._owner = self
        self.residue_frequencies._n = self.residue_frequencies._shape[0] = self.alphabet.K

    def __dealloc__(self):
        libhmmer.p7_bg.p7_bg_Destroy(self._bg)

    def __copy__(self):
        return self.copy()

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
        """copy(self)\n--

        Create a copy of the null model with the same parameters.

        """
        assert self._bg != NULL

        cdef Background new = Background.__new__(Background)
        with nogil:
            new._bg = libhmmer.p7_bg.p7_bg_Clone(self._bg)
        if new._bg == NULL:
            raise AllocationError("P7_BG")

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

    .. versionadded:: 0.2.0

    .. versionchanged:: 0.4.2
       Added the `~Builder.randomness` attribute.

    """

    _ARCHITECTURE_STRATEGY = {
        "fast": p7_archchoice_e.p7_ARCH_FAST,
        "hand": p7_archchoice_e.p7_ARCH_HAND,
    }

    _WEIGHTING_STRATEGY = {
        "pb": p7_wgtchoice_e.p7_WGT_PB,
        "gsc": p7_wgtchoice_e.p7_WGT_GSC,
        "blosum": p7_wgtchoice_e.p7_WGT_BLOSUM,
        "none": p7_wgtchoice_e.p7_WGT_NONE,
        "given": p7_wgtchoice_e.p7_WGT_GIVEN,
    }

    _EFFECTIVE_STRATEGY = {
        "entropy": p7_effnchoice_e.p7_EFFN_ENTROPY,
        "exp": p7_effnchoice_e.p7_EFFN_ENTROPY_EXP,
        "clust": p7_effnchoice_e.p7_EFFN_CLUST,
        "none": p7_effnchoice_e.p7_EFFN_NONE,
    }

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
    ):
        """__init__(self, alphabet, *, architecture="fast", weighting="pb", effective_number="entropy", prior_scheme="alpha", symfrac=0.5, fragthresh=0.5, wid=0.62, esigma=45.0, eid=0.62, EmL=200, EmN=200, EvL=200, EvN=200, EfL=100, EfN=200, Eft=0.04, seed=42, ere=None, popen=None, pextend=None)\n--

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
                parameterizing  from counts, either ``"laplace"`` (the default)
                or ``"alphabet"``.
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
            popen (`float`): The *gap open* probability to use with the score
                system. Default depends on the alphabet: *0.02* for proteins,
                *0.03125* for nucleotides.
            pextend (`float`): The *gap extend* probability to use with the
                score system. Default depends on the alphabet: *0.4* for
                proteins, *0.75* for nucleotides.

        """
        # extract alphabet and create raw P7_BUILDER object
        self.alphabet = alphabet
        abcty = alphabet._abc.type
        with nogil:
            self._bld = libhmmer.p7_builder.p7_builder_Create(NULL, alphabet._abc)
        if self._bld == NULL:
            raise AllocationError("P7_BG")

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
        _arch = self._ARCHITECTURE_STRATEGY.get(architecture)
        if _arch is not None:
            self._bld.arch_strategy = _arch
        else:
            raise ValueError(f"Invalid value for 'architecture': {architecture}")

        # set the weighting strategy
        self.weighting = weighting
        _weighting = self._WEIGHTING_STRATEGY.get(weighting)
        if _weighting is not None:
            self._bld.wgt_strategy = _weighting
        else:
            raise ValueError(f"Invalid value for 'weighting': {weighting}")

        # set the effective sequence number strategy
        self.effective_number = effective_number
        if isinstance(effective_number, (int, float)):
            self._bld.effn_strategy = p7_effnchoice_e.p7_EFFN_SET
            self._bld.eset = effective_number
        elif isinstance(effective_number, str):
            _effn = self._EFFECTIVE_STRATEGY.get(effective_number)
            if _effn is not None:
                self._bld.effn_strategy = _effn
            else:
                raise ValueError(f"Invalid value for 'effective_number': {effective_number}")
        else:
            ty = type(effective_number).__name__
            raise TypeError(f"Invalid type for 'effective_number': {ty}")

        # set the RE target if given one
        self._bld.re_target = -1 if ere is None else ere

        # set the prior scheme
        self.prior_scheme = prior_scheme
        with nogil:
            if prior_scheme is None:
                self._bld.prior = NULL
            elif prior_scheme == "laplace":
                self._bld.prior = libhmmer.p7_prior.p7_prior_CreateLaplace(self.alphabet._abc)
            elif prior_scheme == "alphabet":
                if abcty == libeasel.alphabet.eslAMINO:
                    self._bld.prior = libhmmer.p7_prior.p7_prior_CreateAmino()
                elif abcty == libeasel.alphabet.eslDNA or abcty == libeasel.alphabet.eslRNA:
                    self._bld.prior = libhmmer.p7_prior.p7_prior_CreateNucleic()
                else:
                    self._bld.prior = libhmmer.p7_prior.p7_prior_CreateLaplace(self.alphabet._abc)
            else:
                raise ValueError("Invalid value for 'prior_scheme': {prior_scheme}")

        # set the probability using alphabet-specific defaults or given values
        if abcty == libeasel.alphabet.eslDNA or abcty == libeasel.alphabet.eslRNA:
            self.popen = 0.03125 if popen is None else popen
            self.pextend = 0.75 if pextend is None else pextend
        else:
            self.popen = 0.02 if popen is None else popen
            self.pextend = 0.4 if pextend is None else pextend

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
        self.randomness._seed(seed)

    # --- Methods ------------------------------------------------------------

    cpdef tuple build(
        self,
        DigitalSequence sequence,
        Background background,
    ):
        """build(self, sequence, background)\n--

        Build a new HMM from ``sequence`` using the builder configuration.

        Arguments:
            sequence (`~pyhmmer.easel.DigitalSequence`): A single biological
                sequence in digital mode to build a HMM with.
            background (`pyhmmer.plan7.background`): The background model
                to use to create the HMM.

        Raises:
            `~pyhmmer.errors.AlphabetMismatch`: When either ``sequence`` or
                ``background`` have the wrong alphabet for this builder.

        .. versionchanged:: 0.4.6
           Sets the `HMM.creation_time` attribute with the current time.

        """
        assert self._bld != NULL

        cdef int              status
        cdef HMM              hmm     = HMM.__new__(HMM)
        cdef Profile          profile = Profile.__new__(Profile)
        cdef OptimizedProfile opti    = OptimizedProfile.__new__(OptimizedProfile)
        cdef str              msg

        # use given probabilities
        self._bld.popen = self.popen
        self._bld.pextend = self.pextend

        # reseed RNG used by the builder if needed
        if self._bld.do_reseeding:
            self.randomness._seed(self.seed)

        # check alphabet and assign it to newly created objects
        hmm.alphabet = profile.alphabet = opti.alphabet = self.alphabet
        if not self.alphabet._eq(background.alphabet):
            raise AlphabetMismatch(self.alphabet, background.alphabet)
        if not self.alphabet._eq(sequence.alphabet):
            raise AlphabetMismatch(self.alphabet, sequence.alphabet)

        # load score matrix and build HMM
        # TODO: allow changing from the default scoring matrix
        # TODO: allow caching the parameter values to avoid resetting
        #       everytime `build` is called.
        with nogil:
            status = libhmmer.p7_builder.p7_builder_SetScoreSystem(
                self._bld,
                NULL, # --mxfile
                NULL, # env
                self.popen, # popen
                self.pextend,  # pextend
                background._bg
            )
            if status == libeasel.eslENOTFOUND:
                raise FileNotFoundError("could not open substitution score matrix file")
            elif status == libeasel.eslEINVAL:
                raise ValueError("cannot convert matrix to conditional probabilities")
            elif status != libeasel.eslOK:
                raise UnexpectedError(status, "p7_builder_SetScoreSystem")
            # build HMM and profiles
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
            raise ValueError("Could not build HMM: {}".format(msg))
        else:
            raise UnexpectedError(status, "p7_SingleBuilder")

        # setup attributes
        hmm.evalue_parameters = _EvalueParameters(hmm)
        hmm.cutoffs = _Cutoffs(hmm)
        profile.offsets = _Offsets(profile)
        profile.evalue_parameters = _EvalueParameters(profile)
        profile.cutoffs = _Cutoffs(profile)
        opti.offsets = _Offsets(opti)
        # return newly built HMM, profile and optimized profile
        return hmm, profile, opti

    cpdef tuple build_msa(
        self,
        DigitalMSA msa,
        Background background,
    ):
        """build_msa(self, msa, background)\n--

        Build a new HMM from ``msa`` using the builder configuration.

        Arguments:
            msa (`~pyhmmer.easel.DigitalMSA`): A multiple sequence
                alignment in digital mode to build a HMM with.
            background (`pyhmmer.plan7.background`): The background model
                to use to create the HMM.

        Raises:
            `~pyhmmer.errors.AlphabetMismatch`: When either ``msa`` or
                ``background`` have the wrong alphabet for this builder.

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

        # use given probabilities
        self._bld.popen = self.popen
        self._bld.pextend = self.pextend

        # reseed RNG used by the builder if needed
        if self._bld.do_reseeding:
            self.randomness._seed(self.seed)

        # check alphabet and assign it to newly created objects
        hmm.alphabet = profile.alphabet = opti.alphabet = self.alphabet
        if not self.alphabet._eq(background.alphabet):
            raise AlphabetMismatch(self.alphabet, background.alphabet)
        if not self.alphabet._eq(msa.alphabet):
            raise AlphabetMismatch(self.alphabet, msa.alphabet)

        # load score matrix and build HMM
        with nogil:
            status = libhmmer.p7_builder.p7_builder_SetScoreSystem(
                self._bld,
                NULL, # --mxfile
                NULL, # env
                self.popen, # popen
                self.pextend,  # pextend
                background._bg
            )
            if status == libeasel.eslENOTFOUND:
                raise FileNotFoundError("could not open substitution score matrix file")
            elif status == libeasel.eslEINVAL:
                raise ValueError("cannot convert matrix to conditional probabilities")
            elif status != libeasel.eslOK:
                raise UnexpectedError(status, "p7_builder_SetScoreSystem")
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
            raise ValueError("Could not build HMM: {}".format(msg))
        else:
            raise UnexpectedError(status, "p7_Builder")

        # setup attributes
        hmm.evalue_parameters = _EvalueParameters(hmm)
        hmm.cutoffs = _Cutoffs(hmm)
        profile.offsets = _Offsets(profile)
        profile.evalue_parameters = _EvalueParameters(profile)
        profile.cutoffs = _Cutoffs(profile)
        opti.offsets = _Offsets(opti)
        # return newly built HMM, profile and optimized profile
        return hmm, profile, opti

    cpdef Builder copy(self):
        """copy(self)\n--

        Create a duplicate `Builder` instance with the same arguments.

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
            pextend=self.pextend
        )


cdef class _Cutoffs:
    """A read-only view over the Pfam score cutoffs of a `HMM` or a `Profile`.

    .. versionadded:: 0.4.6

    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._owner = None
        self._flags = NULL
        self._cutoffs = NULL
        self._is_profile = True

    def __init__(self, object owner):
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
            raise TypeError("expected Profile or HMM, found {ty}")

    def __copy__(self):
        cdef _Cutoffs c
        c = _Cutoffs.__new__(_Cutoffs)
        c._owner = self._owner
        c._cutoffs = self._cutoffs
        c._flags = self._flags
        c._owner = self._owner
        return c

    def __str__(self):
        ty = type(self).__name__
        return "<Pfam score cutoffs of {!r}>".format(
            self._owner,
        )

    def __eq__(self, object other):
        if isinstance(other, _Cutoffs):
            return self.as_vector() == other.as_vector()
        return NotImplemented

    # --- Properties ---------------------------------------------------------

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
        """gathering_available(self)\n--

        Check whether the gathering thresholds are available.

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
        """trusted_available(self)\n--

        Check whether the trusted cutoffs are available.

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
        """noise_available(self)\n--

        Check whether the noise cutoffs are available.

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
        """as_vector(self)\n --

        Return a view over the score cutoffs as a `~VectorF`.

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

    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self, Hit hit, size_t index):
        self._dom = &hit._hit.dcl[index]
        self.hit = hit
        self.alignment = Alignment(self)

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
            return exp(self._dom.lnP) * self.hit.hits.domZ

    @property
    def i_evalue(self):
        """`float`: The independent e-value for the domain.
        """
        assert self._dom != NULL
        if self.hit.hits.long_targets:
            return exp(self._dom.lnP)
        else:
            return exp(self._dom.lnP) * self.hit.hits.Z

    @property
    def pvalue(self):
        """`float`: The p-value of the domain bitscore.
        """
        assert self._dom != NULL
        return exp(self._dom.lnP)


cdef class Domains:
    """A read-only view over the domains of a single `~pyhmmer.plan7.Hit`.

    Attributes:
        hit (`~pyhmmer.plan7.Hit`): The target hit these domains belong hit.

    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self, Hit hit):
        self.hit = hit

    def __len__(self):
        return self.hit._hit.ndom

    def __getitem__(self, int index):
        if index < 0:
            index += self.hit._hit.ndom
        if index >= self.hit._hit.ndom or index < 0:
            raise IndexError("list index out of range")
        return Domain(self.hit, <size_t> index)


cdef class _EvalueParameters:
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
        cdef str ty
        if isinstance(owner, Profile):
            self._evparams = &(<Profile> owner)._gm.evparam
            self._owner = owner
        elif isinstance(owner, HMM):
            self._evparams = &(<HMM> owner)._hmm.evparam
            self._owner = owner
        else:
            ty = type(owner).__name__
            raise TypeError("expected Profile or HMM, found {ty}")

    def __copy__(self):
        cdef _EvalueParameters ev
        ev = _EvalueParameters.__new__(_EvalueParameters)
        ev._owner = self._owner
        ev._evparams = self._evparams
        return ev

    def __str__(self):
        ty = type(self).__name__
        return "<e-value parameters of {!r}>".format(
            self._owner,
        )

    def __eq__(self, object other):
        if isinstance(other, _EvalueParameters):
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
        """as_vector(self)\n--

        Return a view over the e-value parameters as a `~VectorF`.

        """
        assert self._evparams != NULL

        cdef VectorF new

        new = VectorF.__new__(VectorF)
        new._owner = self
        new._n = new._shape[0] = p7_NEVPARAM
        new._data = <float*> self._evparams
        return new


cdef class Hit:
    """A high-scoring database hit found by the comparison pipeline.
    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self, TopHits hits, size_t index):
        assert hits._th != NULL
        assert index < hits._th.N
        self.hits = hits
        self._hit = hits._th.hit[index]

    # --- Properties ---------------------------------------------------------

    @property
    def name(self):
        """`bytes`: The name of the database hit.
        """
        assert self._hit != NULL
        assert self._hit.name != NULL
        return <bytes> self._hit.name

    @property
    def accession(self):
        """`bytes` or `None`: The accession of the database hit, if any.
        """
        assert self._hit != NULL
        if self._hit.acc == NULL:
            return None
        return <bytes> self._hit.acc

    @property
    def description(self):
        """`bytes` or `None`: The description of the database hit, if any.
        """
        assert self._hit != NULL
        if self._hit.desc == NULL:
            return None
        return <bytes> self._hit.acc

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
        if self.hits.long_targets:
            return exp(self._hit.lnP)
        else:
            return exp(self._hit.lnP) * self.hits.Z

    @property
    def pvalue(self):
        """`float`: The p-value of the bitscore.

        .. versionadded:: 0.4.2

        """
        assert self._hit != NULL
        return exp(self._hit.lnP)

    # --- Methods ------------------------------------------------------------

    cpdef bint is_included(self):
        """is_included(self)\n--

        Check if the hit should be included with respect to the thresholds.

        """
        return self._hit.flags & p7_hitflags_e.p7_IS_INCLUDED != 0

    cpdef bint is_reported(self):
        """is_reported(self)\n--

        Check if the hit should be reported with respect to the thresholds.

        """
        return self._hit.flags & p7_hitflags_e.p7_IS_REPORTED != 0

    cpdef bint is_new(self):
        """is_new(self)\n--

        Check if the hit is a new hit.

        """
        return self._hit.flags & p7_hitflags_e.p7_IS_NEW != 0

    cpdef bint is_dropped(self):
        """is_dropped(self)\n--

        Check if the hit was dropped.

        """
        return self._hit.flags & p7_hitflags_e.p7_IS_DROPPED != 0

    cpdef bint is_duplicate(self):
        """is_duplicate(self)\n--

        Check if the hit is a duplicate.

        """
        return self._hit.flags & p7_hitflags_e.p7_IS_DUPLICATE


@cython.freelist(8)
cdef class HMM:
    """A data structure storing the Plan7 Hidden Markov Model.

    Attributes:
        alphabet (`~pyhmmer.easel.Alphabet`): The alphabet of the model.
        evalue_parameters (`~pyhmmer.plan7._EvalueParameters`): The e-value
            parameters for this HMM.
        cutoffs (`~pyhmmer.plan7._Cutoffs`): The Pfam score cutoffs for this
            HMM, if any.

    .. versionchanged:: 0.4.6
       Added the `~HMM.evalue_parameters` and `~HMM.cutoffs` attributes.

    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self.alphabet = None
        self.cutoffs = None
        self.evalue_parameters = None
        self._hmm = NULL

    def __init__(self, int M, Alphabet alphabet):
        """__init__(self, M, alphabet)\n--

        Create a new HMM from scratch.

        Arguments:
            M (`int`): The length of the model (i.e. the number of nodes).
            alphabet (`~pyhmmer.easel.Alphabet`): The alphabet of the model.

        """
        # store the alphabet so it's not deallocated
        self.alphabet = alphabet
        # create a new HMM suitable for at least M nodes
        with nogil:
            self._hmm = libhmmer.p7_hmm.p7_hmm_Create(M, alphabet._abc)
        if not self._hmm:
            raise AllocationError("P7_HMM")
        # FIXME(@althonos): Remove following block when
        # https://github.com/EddyRivasLab/hmmer/pull/236
        # is merged and released in a new HMMER version
        cdef int status = libhmmer.p7_hmm.p7_hmm_SetConsensus(self._hmm, NULL)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_hmm_SetConsensus")
        self._hmm.flags &= ~libhmmer.p7_hmm.p7H_CONS
        # expose the e-value parameters as the `evalue_parameters` attribute
        self.evalue_parameters = _EvalueParameters.__new__(_EvalueParameters)
        self.evalue_parameters._evparams = &self._hmm.evparam
        self.evalue_parameters._owner = self
        # expose the score cutoffs as the `cutoffs` attribute
        self.cutoffs = _Cutoffs.__new__(_Cutoffs)
        self.cutoffs._cutoffs = &self._hmm.cutoff
        self.cutoffs._flags = &self._hmm.flags
        self.cutoffs._is_profile = False

    def __dealloc__(self):
        libhmmer.p7_hmm.p7_hmm_Destroy(self._hmm)

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
        return HMM, (self.M, self.alphabet), self.__getstate__()

    cpdef dict __getstate__(self):
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
            "flags": self._hmm.flags
        }

        # copy alignment map if available
        if self._hmm.flags & libhmmer.p7_hmm.p7H_MAP:
            assert self._hmm.map != NULL
            state["map"] = map_ = array.array('i')
            for i in range(self._hmm.M + 1):
                map_.append(self._hmm.map[i])
        # add optional annotations
        if self._hmm.flags & libhmmer.p7_hmm.p7H_RF:
            assert self._hmm.rf != NULL
            state["rf"] = PyBytes_FromStringAndSize(self._hmm.rf, self._hmm.M + 1)
        if self._hmm.flags & libhmmer.p7_hmm.p7H_MM:
            assert self._hmm.mm != NULL
            state["mm"] = PyBytes_FromStringAndSize(self._hmm.mm, self._hmm.M + 1)
        if self._hmm.flags & libhmmer.p7_hmm.p7H_CS:
            assert self._hmm.cs != NULL
            state["cs"] = PyBytes_FromStringAndSize(self._hmm.cs, self._hmm.M + 1)
        if self._hmm.flags & libhmmer.p7_hmm.p7H_CA:
            assert self._hmm.ca != NULL
            state["ca"] = PyBytes_FromStringAndSize(self._hmm.ca, self._hmm.M + 1)
        if self._hmm.flags & libhmmer.p7_hmm.p7H_CONS:
            assert self._hmm.consensus != NULL
            state["consensus"] = PyBytes_FromStringAndSize(self._hmm.consensus, self._hmm.M + 1)

        return state

    cpdef object __setstate__(self, dict state):
        assert self._hmm != NULL

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

        # strings that must be set manually
        if self._hmm.rf != NULL:
            free(self._hmm.rf)
            self._hmm.rf = NULL
        if self._hmm.flags & libhmmer.p7_hmm.p7H_RF:
            self._hmm.rf = strndup(<const char*> state["rf"], M+2)
            if self._hmm.rf == NULL:
                raise AllocationError("char*")
        if self._hmm.mm != NULL:
            free(self._hmm.mm)
            self._hmm.mm = NULL
        if self._hmm.flags & libhmmer.p7_hmm.p7H_MM:
            self._hmm.mm = strndup(<const char*> state["mm"], M+2)
            if self._hmm.mm == NULL:
                raise AllocationError("char*")
        if self._hmm.consensus != NULL:
            free(self._hmm.consensus)
            self._hmm.consensus = NULL
        if self._hmm.flags & libhmmer.p7_hmm.p7H_CONS:
            self._hmm.consensus = strndup(<const char*> state["consensus"], M+2)
            if self._hmm.consensus == NULL:
                raise AllocationError("char*")
        if self._hmm.cs != NULL:
            free(self._hmm.cs)
            self._hmm.cs = NULL
        if self._hmm.flags & libhmmer.p7_hmm.p7H_CS:
            self._hmm.cs = strndup(<const char*> state["cs"], M+2)
            if self._hmm.cs == NULL:
                raise AllocationError("char*")
        if self._hmm.ca != NULL:
            free(self._hmm.ca)
            self._hmm.ca = NULL
        if self._hmm.flags & libhmmer.p7_hmm.p7H_CA:
            self._hmm.ca = strndup(<const char*> state["ca"], M+2)
            if self._hmm.ca == NULL:
                raise AllocationError("char*")

        # numerical values that can be set directly on the C struct
        self._hmm.nseq = state["nseq"]
        self._hmm.eff_nseq = state["eff_nseq"]
        self._hmm.max_length = state["max_length"]
        self._hmm.checksum = state["checksum"]
        self._hmm.offset = state["offset"]

        # vector data that must be copied
        assert self._hmm.ins != NULL and self._hmm.ins[0] != NULL
        assert self._hmm.mat != NULL and self._hmm.mat[0] != NULL
        assert self._hmm.t != NULL and self._hmm.t[0] != NULL
        assert ins.ndim == 2 and ins.shape[0] == M+1 and ins.shape[1] == K
        assert mat.ndim == 2 and mat.shape[0] == M+1 and mat.shape[1] == K
        assert t.ndim == 2 and t.shape[0] == M+1 and t.shape[1] == p7H_NTRANSITIONS
        memcpy(self._hmm.ins[0],  &ins[0][0], (M+1) * K                * sizeof(float))
        memcpy(self._hmm.mat[0],  &mat[0][0], (M+1) * K                * sizeof(float))
        memcpy(self._hmm.t[0],    &t[0][0],   (M+1) * p7H_NTRANSITIONS * sizeof(float))

        # arrays that must be copied (compo may be None in the stat dictionary)
        assert evparam.ndim == 1 and evparam.shape[0] == p7_NEVPARAM
        assert cutoff.ndim == 1 and cutoff.shape[0] == p7_NCUTOFFS
        memcpy(&self._hmm.evparam[0], &evparam[0], p7_NEVPARAM * sizeof(float))
        memcpy(&self._hmm.cutoff[0], &cutoff[0], p7_NCUTOFFS * sizeof(float))
        if self._hmm.flags & libhmmer.p7_hmm.p7H_COMPO:
            compo = state["compo"]
            assert compo.ndim == 1 and compo.shape[0] == K
            memcpy(&self._hmm.cutoff[0], &compo[0], K * sizeof(float))

        # copy alignment map only if it is available
        if self._hmm.flags & p7H_MAP == 0:
            if self._hmm.map != NULL:
                free(self._hmm.map)
                self._hmm.map = NULL
        else:
            map_ = state["map"]
            if self._hmm.map != NULL:
                self._hmm.map = <int*> realloc(self._hmm.map, (M + 1) * sizeof(int))
            else:
                self._hmm.map = <int*> calloc(M + 1, sizeof(int))
            if self._hmm.map == NULL:
                raise AllocationError("int*")
            assert map_.ndim == 1 and map_.shape[0] == M + 1
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
        """`bytes` or `None`: The name of the HMM, if any.
        """
        assert self._hmm != NULL
        return None if self._hmm.name == NULL else <bytes> self._hmm.name

    @name.setter
    def name(self, bytes name):
        assert self._hmm != NULL

        cdef char* name_ = NULL if name is None else <char*> name
        cdef int err = libhmmer.p7_hmm.p7_hmm_SetName(self._hmm, name_)

        if err == libeasel.eslEMEM:
            raise AllocationError("char*")
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

        cdef char* acc = NULL if accession is None else <char*> accession
        cdef int err = libhmmer.p7_hmm.p7_hmm_SetAccession(self._hmm, acc)

        if err == libeasel.eslEMEM:
            raise AllocationError("char*")
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

                >>> dna = easel.Alphabet.dna()  # dna.K=4
                >>> hmm = plan7.HMM(100, dna)
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
        return (&self._hmm.consensus[1]).decode("ascii")

    @property
    def consensus_structure(self):
        """`str` or `None`: The consensus structure of the HMM, if any.

        .. versionadded:: 0.3.1

        """
        assert self._hmm != NULL
        if not (self._hmm.flags & libhmmer.p7_hmm.p7H_CS):
            return None
        assert self._hmm.cs != NULL
        return (&self._hmm.cs[1]).decode("ascii")

    @property
    def consensus_accessibility(self):
        """`str` or `None`: The consensus accessibility of the HMM, if any.

        .. versionadded:: 0.3.1

        """
        assert self._hmm != NULL
        if not (self._hmm.flags & libhmmer.p7_hmm.p7H_CA):
            return None
        assert self._hmm.ca != NULL
        return (&self._hmm.ca[1]).decode("ascii")

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
        return (&self._hmm.rf[1]).decode("ascii")

    @property
    def model_mask(self):
        """`str` or `None`: The model mask line from the alignment, if any.

        .. versionadded:: 0.3.1

        """
        assert self._hmm != NULL
        if not (self._hmm.flags & libhmmer.p7_hmm.p7H_MM):
            return None
        assert self._hmm.mm != NULL
        return (&self._hmm.mm[1]).decode("ascii")

    @property
    def description(self):
        """`bytes` or `None`: The description of the HMM, if any.
        """
        assert self._hmm != NULL
        return None if self._hmm.desc == NULL else <bytes> self._hmm.desc

    @description.setter
    def description(self, bytes description):
        assert self._hmm != NULL

        cdef char* desc = NULL if description is None else <char*> description
        cdef int err = libhmmer.p7_hmm.p7_hmm_SetDescription(self._hmm, desc)

        if err == libeasel.eslEMEM:
            raise AllocationError("char*")
        elif err != libeasel.eslOK:
            raise UnexpectedError(err, "p7_hmm_SetDescription")

    @property
    def transition_probabilities(self):
        r"""`~pyhmmer.easel.MatrixF`: The transition probabilities of the model.

        The property exposes a matrix of shape :math:`(M+1, 7)`, with one row
        per node (plus one extra row for the entry probabilities), and one
        column per transition.

        Columns indices correspond to the following: 0 for :math:`M_n \to M_{n+1}``,
        1 for :math:`M_n \to I_{n+1}`, 2 for :math:`M_n \to M_{n+1}`,
        3 for :math:`I_n \to M_{n+1}`, 4 for :math:`I_n \to I_{n+1}`,
        5 for :math:`D_n \to D_{n+1}`, 6 for :math:`D_n \to D_{n+1}`.

        Example:
            >>> t = thioesterase.transition_probabilities
            >>> t[0, 5]  # TDM, 1 by convention
            1.0

        Caution:
            If editing this matrix manually, note that some invariants need
            to hold for the HMM to be valid: :math:`I_n`, :math:`M_n` and
            :math:`D_n` transition probabilities should only contain
            probabilities between 0 and 1, and sum to 1::

                >>> t = thioesterase.transition_probabilities
                >>> t[50, 0] + t[50, 1] + t[50, 2]  # M_n probabilities
                1.000...
                >>> t[50, 3] + t[50, 4]  # I_n probabilities
                1.000...
                >>> t[50, 5] + t[50, 6]  # D_n probabilties
                1.000...

        .. versionadded:: 0.3.1

        .. versionchanged:: 0.4.0
            This property is now a `~pyhmmer.easel.MatrixF`.

        """
        assert self._hmm != NULL
        cdef MatrixF mat = MatrixF.__new__(MatrixF)
        mat._m = mat._shape[0] = self._hmm.M + 1
        mat._n = mat._shape[1] = libhmmer.p7_hmm.p7H_NTRANSITIONS
        mat._owner = self
        mat._data = self._hmm.t
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

                >>> hmm = HMM(100, alphabet=easel.Alphabet.dna())
                >>> hmm.match_emissions[0]
                VectorF([1.0, 0.0, 0.0, 0.0])

        Caution:
            If editing this matrix manually, note that rows must contain
            valid probabilities for the HMM to be valid: each row must
            contains values between 0 and 1, and sum to 1.

        .. versionadded:: 0.3.1

        .. versionchanged:: 0.4.0
            This property is now a `~pyhmmer.easel.MatrixF`, and stores row 0.

        """
        assert self._hmm != NULL
        cdef MatrixF mat = MatrixF.__new__(MatrixF)
        mat._m = mat._shape[0] = self._hmm.M + 1
        mat._n = mat._shape[1] = self.alphabet.K
        mat._owner = self
        mat._data = self._hmm.mat
        return mat

    @property
    def insert_emissions(self):
        """`memoryview` of `float`: The insert emissions of the model.

        The property exposes a matrix of shape :math:`(M+1, K)`, with one
        row per node and one column per alphabet symbol.

        Caution:
            If editing this matrix manually, note that rows must contain
            valid probabilities for the HMM to be valid: each row must
            contains values between 0 and 1, and sum to 1.

        .. versionadded:: 0.3.1

        .. versionchanged:: 0.4.0
            This property is now a `~pyhmmer.easel.MatrixF`, and stores row 0.

        """
        assert self._hmm != NULL
        cdef MatrixF mat = MatrixF.__new__(MatrixF)
        mat._m = mat._shape[0] = self._hmm.M + 1
        mat._n = mat._shape[1] = self.alphabet.K
        mat._owner = self
        mat._data = self._hmm.ins
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
            raise AllocationError("char*")

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
            raise AllocationError("char*")
        if self._hmm.ctime != NULL:
            strncpy(self._hmm.ctime, <const char*> formatted, n + 1)

    # --- Methods ------------------------------------------------------------

    cpdef void write(self, object fh, bint binary=False) except *:
        """write(self, fh, binary=False)\n--

        Write the HMM to a file handle.

        Arguments:
            fh (`io.IOBase`): A Python file handle, opened in binary mode
                (this must be the case even with ``binary=False``, since
                the C code will emit bytes in either case).
            binary (`bool`): Pass ``False`` to emit the file in ASCII mode
                using the latest supported HMMER format, or ``True`` to use
                the binary HMMER3 format.

        """
        cdef int     status
        cdef FILE*   file
        cdef P7_HMM* hm     = self._hmm

        file = fopen_obj(fh, mode="w")

        if binary:
            status = libhmmer.p7_hmmfile.p7_hmmfile_WriteBinary(file, -1, hm)
        else:
            status = libhmmer.p7_hmmfile.p7_hmmfile_WriteASCII(file, -1, hm)

        if status == libeasel.eslOK:
            fclose(file)
        else:
            raise UnexpectedError(status, "p7_hmmfile_WriteASCII")

    cpdef void zero(self):
        """zero(self)\n--

        Set all parameters to zero, including model composition.

        """
        assert self._hmm != NULL
        with nogil:
            libhmmer.p7_hmm.p7_hmm_Zero(self._hmm)

    cpdef HMM copy(self):
        """copy(self)\n--

        Return a copy of the HMM with the exact same configuration.

        .. versionadded:: 0.3.0

        """
        assert self._hmm != NULL

        cdef HMM new = HMM.__new__(HMM)
        new.alphabet = self.alphabet
        new.evalue_parameters = _EvalueParameters(new)
        new.cutoffs = _Cutoffs(new)

        with nogil:
            new._hmm = libhmmer.p7_hmm.p7_hmm_Clone(self._hmm)
        if not new._hmm:
            raise AllocationError("P7_HMM")
        return new

    cpdef void renormalize(self):
        """renormalize(self)\n--

        Renormalize all parameter vectors (emissions and transitions).

        .. versionadded:: 0.4.0

        """
        assert self._hmm != NULL

        cdef int status
        with nogil:
            status = libhmmer.p7_hmm.p7_hmm_Renormalize(self._hmm)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_hmm_Renormalize")

    cpdef void scale(self, double scale, bint exponential=False):
        """scale(self, scale, exponential=False)\n--

        In a model containing counts, rescale counts by a factor.

        This method only affects core probability model emissions and
        transitions.

        Arguments:
            scale (`float`): The scaling factor to use (:math:`1.0` for no
                scaling). Often computed using the ratio of effective
                sequences (:math:`\\frac{n_{eff}}{n_{seq}}`)
            exponential (`bool`): When set to `True`, use ``scale`` as an
                exponential factor (:math:`C_i = C_i \exp{s}`) instead of
                a multiplicative factor (:math: `C_i = C_i \\times s`),
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

    cpdef void set_composition(self):
        """set_composition(self)\n--

        Calculate and set the model composition.

        .. versionadded:: 0.4.0

        """
        assert self._hmm != NULL

        cdef int status
        with nogil:
            status = libhmmer.p7_hmm.p7_hmm_SetComposition(self._hmm)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_hmm_SetComposition")


cdef class HMMFile:
    """A wrapper around a file (or database), storing serialized HMMs.
    """

    _FORMATS = {
        "2.0": p7_hmmfile_formats_e.p7_HMMFILE_20,
        "3/a": p7_hmmfile_formats_e.p7_HMMFILE_3a,
        "3/b": p7_hmmfile_formats_e.p7_HMMFILE_3b,
        "3/c": p7_hmmfile_formats_e.p7_HMMFILE_3c,
        "3/d": p7_hmmfile_formats_e.p7_HMMFILE_3d,
        "3/e": p7_hmmfile_formats_e.p7_HMMFILE_3e,
        "3/f": p7_hmmfile_formats_e.p7_HMMFILE_3f,
    }

    _MAGIC = {
        v3a_magic: p7_hmmfile_formats_e.p7_HMMFILE_3a,
        v3b_magic: p7_hmmfile_formats_e.p7_HMMFILE_3b,
        v3c_magic: p7_hmmfile_formats_e.p7_HMMFILE_3c,
        v3d_magic: p7_hmmfile_formats_e.p7_HMMFILE_3d,
        v3e_magic: p7_hmmfile_formats_e.p7_HMMFILE_3e,
        v3f_magic: p7_hmmfile_formats_e.p7_HMMFILE_3f,
    }

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._alphabet = None
        self._hfp = NULL

    def __init__(self, object file, bint db = True):
        """__init__(self, file, db=True)\n--

        Create a new HMM reader from the given file.

        Arguments:
            file (`str` or file-like object): Either the path to a file
                containing the HMMs to read, or a file-like object opened in
                binary-mode.
            db (`bool`): Set to `False` to force the parser to ignore the
                pressed HMM database if it finds one. Defaults to `False`.

        """
        cdef int       status
        cdef bytes     fspath
        cdef bytearray errbuf = bytearray(eslERRBUFSIZE)

        try:
            fspath = os.fsencode(file)
            if db:
                function = "p7_hmmfile_OpenE"
                status = libhmmer.p7_hmmfile.p7_hmmfile_OpenE(fspath, NULL, &self._hfp, errbuf)
            else:
                function = "p7_hmmfile_OpenENoDB"
                status = libhmmer.p7_hmmfile.p7_hmmfile_OpenENoDB(fspath, NULL, &self._hfp, errbuf)
        except TypeError:
            self._hfp = self._open_fileobj(file)
            status    = libeasel.eslOK

        if status == libeasel.eslENOTFOUND:
            raise FileNotFoundError(errno.ENOENT, "no such file or directory: {!r}".format(file))
        elif status == libeasel.eslEFORMAT:
            if os.stat(file).st_size:
                raise ValueError("format not recognized by HMMER")
            else:
                raise EOFError("HMM file is empty")
        elif status != libeasel.eslOK:
            raise UnexpectedError(status, function)

        self._alphabet = Alphabet.__new__(Alphabet)
        self._alphabet._abc = NULL

    def __dealloc__(self):
        if self._hfp:
            warnings.warn("unclosed HMM file", ResourceWarning)
            self.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def __iter__(self):
        return self

    def __next__(self):
        cdef int status
        cdef HMM py_hmm
        cdef P7_HMM* hmm = NULL

        if self._hfp == NULL:
            raise ValueError("I/O operation on closed file.")

        # with nogil:
        status = libhmmer.p7_hmmfile.p7_hmmfile_Read(self._hfp, &self._alphabet._abc, &hmm)

        if status == libeasel.eslOK:
            py_hmm = HMM.__new__(HMM)
            py_hmm.alphabet = self._alphabet # keep a reference to the alphabet
            py_hmm._hmm = hmm
            py_hmm.evalue_parameters = _EvalueParameters(py_hmm)
            py_hmm.cutoffs = _Cutoffs(py_hmm)
            return py_hmm
        elif status == libeasel.eslEOF:
            raise StopIteration()
        elif status == libeasel.eslEMEM:
            raise AllocationError("P7_HMM")
        elif status == libeasel.eslESYS:
            raise OSError(self._hfp.errbuf.decode("utf-8", "replace"))
        elif status == libeasel.eslEFORMAT:
            raise ValueError("Invalid format in file: {}".format(self._hfp.errbuf.decode("utf-8", "replace")))
        elif status == libeasel.eslEINCOMPAT:
            alphabet = libeasel.alphabet.esl_abc_DecodeType(self._alphabet.type)
            raise ValueError("HMM is not in the expected {} alphabet".format(alphabet))
        else:
            raise UnexpectedError(status, "p7_hmmfile_Read")

    # --- Methods ------------------------------------------------------------

    cpdef void close(self):
        """close(self)\n--

        Close the HMM file and free resources.

        This method has no effect if the file is already closed. It is called
        automatically if the `HMMFile` was used in a context::

            >>> with HMMFile("tests/data/hmms/bin/PKSI-AT.h3m") as hmm_file:
            ...     hmm = next(hmm_file)
        """
        if self._hfp:
            libhmmer.p7_hmmfile.p7_hmmfile_Close(self._hfp)
            self._hfp = NULL

    # --- Utils --------------------------------------------------------------

    cdef P7_HMMFILE* _open_fileobj(self, object fh) except *:
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
            raise AllocationError("P7_HMMFILE")

        # store options
        hfp.f            = fopen_obj(fh_)
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
        if hasattr(fh, "name"):
            filename = fh.name.encode()
            hfp.fname = strdup(filename)
            if hfp.fname == NULL:
                raise AllocationError("char*")

        # check if the parser is in binary format,
        magic = int.from_bytes(fh_.peek(4)[:4], sys.byteorder)
        if magic in self._MAGIC:
            hfp.format = self._MAGIC[magic]
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
            raise AllocationError("ESL_FILEPARSER")
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
            if format in self._FORMATS:
                hfp.format = self._FORMATS[format]
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


cdef class OptimizedProfile:
    """An optimized profile that uses platform-specific instructions.

    Attributes:
        alphabet (`~pyhmmer.easel.Alphabet`): The alphabet for which this
            optimized profile is configured.
        offsets (`~pyhmmer.plan7._Offsets`): The disk offsets for this
            optimized profile, if it was loaded from a pressed HMM file.

    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._om = NULL
        self.alphabet = None

    def __init__(self, int M, Alphabet alphabet):
        """__init__(self, M, alphabet)\n--

        Create a new optimized profile from scratch.

        Optimized profiles use platform-specific code to accelerate the
        various algorithms. Although you can allocate an optimized profile
        yourself, the only way to obtain a fully configured profile is to
        create it with the `Profile.optimized` method, after having
        configured the profile for a given HMM with `Profile.configure`.

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
            raise AllocationError("P7_OPROFILE")
        # expose the disk offsets as the `offsets` attribute
        self.offsets = _Offsets.__new__(_Offsets)
        self.offsets._offs = &self._om.offs
        self.offsets._owner = self

    def __dealloc__(self):
        p7_oprofile.p7_oprofile_Destroy(self._om)

    def __copy__(self):
        return self.copy()

    # --- Properties ---------------------------------------------------------

    @property
    def M(self):
        """`int`: The number of nodes in the model.

        .. versionadded:: v0.4.0

        """
        assert self._om != NULL
        return self._om.M

    @property
    def L(self):
        """`int`: The currently configured target sequence length.

        .. versionadded:: v0.4.0

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
    def tbm(self):
        r"""`int`: The constant cost for a :math:`B \to M_k` transition.

        .. versionadded:: v0.4.0

        """
        assert self._om != NULL
        return self._om.tbm_b

    @property
    def tec(self):
        r"""`int`: The constant cost for a :math:`E \to C` transition.

        .. versionadded:: v0.4.0

        """
        assert self._om != NULL
        return self._om.tec_b

    @property
    def tjb(self):
        """`int`: The constant cost for a :math:`NJC` move.

        .. versionadded:: v0.4.0

        """
        assert self._om != NULL
        return self._om.tjb_b

    @property
    def base(self):
        assert self._om != NULL
        return self._om.base_b

    @property
    def bias(self):
        """`int`: The positive bias to emission scores.

        .. versionadded:: v0.4.0

        """
        assert self._om != NULL
        return self._om.bias_b

    @property
    def sbv(self):
        """`~pyhmmer.easel.MatrixU8`: The match scores for the SSV filter.

        .. versionadded:: v0.4.0

        """
        assert self._om != NULL

        cdef int nqb = p7O_NQB(self._om.M)
        cdef int nqs = nqb + p7O_EXTRA_SB

        cdef MatrixU8 mat = MatrixU8.__new__(MatrixU8)
        mat._m = mat._shape[0] = self.alphabet.Kp
        mat._n = mat._shape[1] = 16 * nqs
        mat._owner = self
        mat._data = <uint8_t**> self._om.sbv
        return mat

    @property
    def rbv(self):
        """`~pyhmmer.easel.MatrixU8`: The match scores for the MSV filter.
        """
        assert self._om != NULL

        cdef MatrixU8 mat = MatrixU8.__new__(MatrixU8)
        mat._m = mat._shape[0] = self.alphabet.Kp
        mat._n = mat._shape[1] = 16 * p7O_NQB(self._om.M)
        mat._owner = self
        mat._data = <uint8_t**> self._om.rbv
        return mat

    # --- Methods ------------------------------------------------------------

    cpdef bint is_local(self):
        """is_local(self)\n--

        Return whether or not the profile is in a local alignment mode.

        """
        assert self._om != NULL
        return p7_oprofile.p7_oprofile_IsLocal(self._om)

    cpdef OptimizedProfile copy(self):
        """copy(self)\n--

        Create an exact copy of the optimized profile.

        """
        assert self._om != NULL
        cdef OptimizedProfile new = OptimizedProfile.__new__(OptimizedProfile)
        new.alphabet = self.alphabet
        new.offsets = _Offsets(new)
        with nogil:
            new._om = p7_oprofile.p7_oprofile_Clone(self._om)
        if new._om == NULL:
            raise AllocationError("P7_OPROFILE")
        return new

    cpdef void write(self, object fh_filter, object fh_profile):
        """write(self, fh_filter, fh_profile)\n--

        Write an optimized profile to two separate files.

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

        pfp = fopen_obj(fh_profile, mode="w")
        ffp = fopen_obj(fh_filter, mode="w")
        status = p7_oprofile_Write(ffp, pfp, self._om)
        if status == libeasel.eslOK:
            fclose(ffp)
            fclose(pfp)
        else:
            raise UnexpectedError(status, "p7_oprofile_Write")

    def convert(self, Profile profile):
        """convert(self, profile)\n--

        Store the given profile into ``self`` as a platform-specific profile.

        Use this method to obtained an optimized profile from a `Profile`
        while recycling the internal vector buffers.

        See Also:
            The `Profile.optimized` method, which allows getting an
            `OptimizedProfile` directly from a profile without having to
            allocate first.

        """
        assert self._om != NULL
        assert profile._gm != NULL

        with nogil:
            self._convert(profile._gm)

    cdef int _convert(self, P7_PROFILE* gm) nogil except 1:
        cdef int status
        if self._om.allocM < gm.M:
            raise ValueError("Optimized profile is too small to hold profile")
        status = p7_oprofile.p7_oprofile_Convert(gm, self._om)
        if status == libeasel.eslOK:
            return 0
        elif status == libeasel.eslEINVAL:
            raise ValueError("Standard and optimized profiles are not compatible.")
        elif status == libeasel.eslEMEM:
            raise AllocationError("P7_OPROFILE")
        elif status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_oprofile_Convert")

    cpdef object ssv_filter(self, DigitalSequence seq):
        """ssv_filter(self, seq)\n--

        Compute the SSV filter score for the given sequence.

        Arguments:
            seq (`~pyhmmer.easel.DigitalSequence`):

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
            `~pyhmmer.errors.AlphabetMismatch`: when the alphabet of the
                sequence does not correspond to the profile alphabet.

        Caution:
            This method is not available on the PowerPC platform (calling
            it will raise a `NotImplementedError`).

        .. versionadded:: v0.4.0

        """
        IF HMMER_IMPL == "SSE":
            assert self._om != NULL

            if self.alphabet != seq.alphabet:
                raise AlphabetMismatch(self.alphabet, seq.alphabet)

            cdef float score
            cdef int status

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
        ELSE:
            raise NotImplementedError("p7_SSVFilter is not available on VMX platforms")


cdef class _Offsets:
    """A mutable view over the disk offsets of a profile.
    """

    def __cinit__(self):
        self._owner = None
        self._offs  = NULL

    def __init__(self, object owner):
        cdef str ty

        if isinstance(owner, Profile):
            self._offs = &(<Profile> owner)._gm.offs
            self._owner = owner
        elif isinstance(owner, OptimizedProfile):
            self._offs = &(<OptimizedProfile> owner)._om.offs
            self._owner = owner
        else:
            ty = type(owner).__name__
            raise TypeError("expected Profile or OptimizedProfile, found {ty}")

    def __copy__(self):
        assert self._offs != NULL
        cdef _Offsets copy = _Offsets.__new__(_Offsets)
        copy._offs = self._offs
        copy._owner = self._owner
        return copy

    def __str__(self):
        ty = type(self).__name__
        return "<offsets of {!r} model={!r} filter={!r} profile={!r}>".format(
            self._owner,
            self.model,
            self.filter,
            self.profile,
        )

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
        self._offs[0][<int> p7_offsets_e.p7_MOFFSET] = -1 if profile is None else profile


cdef class Pipeline:
    """An HMMER3 accelerated sequence/profile comparison pipeline.

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
    LONG_TARGETS = False # always False for now

    _BIT_CUTOFFS = {
        "gathering": libhmmer.p7_hmm.p7H_GA,
        "noise": libhmmer.p7_hmm.p7H_NC,
        "trusted": libhmmer.p7_hmm.p7H_TC,
    }

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._refs = NULL
        self._nref = 0
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
        uint32_t seed=42,
        object Z=None,
        object domZ=None,
        double F1=0.02,
        double F2=1e-3,
        double F3=1e-5,
        double E=10.0,
        double domE=10.0,
        double incE=0.01,
        double incdomE=0.01,
        str bit_cutoffs=None,
    ):
        """__init__(self, alphabet, background=None, *, bias_filter=True, null2=True, seed=42, Z=None, domZ=None, F1=0.02, F2=1e-3, F3=1e-5, E=10.0, domE=10.0, incE=0.01, incdomE=0.01, bit_cutoffs=None)\n--

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
            domE (`float`): The per-domain E-value threshold for reporting
                a domain hit.
            incE (`float`): The per-target E-value threshold for including
                a hit in the resulting `TopHits`.
            incdomE (`float`): The per-domain E-value threshold for including
                a domain in the resulting `TopHits`.
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
        cdef int  m_hint       = self.M_HINT
        cdef int  l_hint       = self.L_HINT
        cdef bint long_targets = self.LONG_TARGETS

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
                long_targets,
                p7_pipemodes_e.p7_SEARCH_SEQS
            )
        if self._pli == NULL:
            raise AllocationError("P7_PIPELINE")

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
        self.domE = domE
        self.incE = incE
        self.incdomE = incdomE

        # setup the model-specific reporting cutoffs
        self._save_cutoff_parameters()
        self.bit_cutoffs = bit_cutoffs

    def __dealloc__(self):
        libhmmer.p7_pipeline.p7_pipeline_Destroy(self._pli)
        free(self._refs)

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
        self.randomness._seed(seed)

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
    def bit_cutoffs(self):
        """`str` or `None`: The model-specific thresholding option, if any.

        .. versionadded:: 0.4.6

        """
        assert self._pli != NULL
        return next(
            (k for k,v in self._BIT_CUTOFFS.items() if v == self._pli.use_bit_cutoffs),
            None
        )

    @bit_cutoffs.setter
    def bit_cutoffs(self, str bit_cutoffs):
        assert self._pli != NULL
        if bit_cutoffs is not None:
            #
            flag = self._BIT_CUTOFFS.get(bit_cutoffs)
            if flag is None:
                raise ValueError(f"invalid bit cutoff: {bit_cutoffs}")
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

    # --- Methods ------------------------------------------------------------

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

    cpdef void clear(self):
        """clear(self)\n--

        Reset the pipeline configuration to its default state.

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
        self.randomness._seed(self._seed)
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
        self._pli.pos_past_fwd    = 0;
        self._pli.pos_output      = 0
        self._pli.strands         = 0
        self._pli.W               = 0
        self._pli.hfp             = NULL
        self._pli.errbuf[0]       = b'\0'

    cpdef TopHits search_hmm(self, HMM query, object sequences):
        """search_hmm(self, query, sequences)\n--

        Run the pipeline using a query HMM against a sequence database.

        Arguments:
            query (`~pyhmmer.plan7.HMM`): The HMM object to use to query the
                sequence database.
            sequences (collection of `~pyhmmer.easel.DigitalSequence`): The
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

        .. versionadded:: 0.2.0

        """
        assert self._pli != NULL

        cdef ssize_t              nseqs
        cdef int                  status
        cdef int                  allocM
        cdef DigitalSequence      seq
        cdef TopHits              hits = TopHits()

        # check the pipeline was configured with the same alphabet
        if not self.alphabet._eq(query.alphabet):
            raise AlphabetMismatch(self.alphabet, query.alphabet)
        # check we can rewind the sequences
        if not isinstance(sequences, collections.abc.Sequence):
            raise TypeError("`sequences` cannot be an iterator, expected an iterable")
        # allow peeking the sequences
        nseqs = len(sequences)
        sequences = peekable(sequences)
        seq = sequences.peek()

        # collect sequences: we prepare an array of pointer to sequences
        # so that the C backend can iter through them without having to
        # acquire the GIL between each iteration.
        if self._nref <= nseqs:
            self._nref = nseqs + 1
            self._refs = <void**> realloc(self._refs, sizeof(void*) * (self._nref))
            if self._refs == NULL:
                raise AllocationError("void**")
        for i, seq in enumerate(sequences):
            # check alphabet
            if not self.alphabet._eq(seq.alphabet):
                raise AlphabetMismatch(self.alphabet, seq.alphabet)
            # check length
            if len(seq) > 100000:
                raise ValueError("sequence length over comparison pipeline limit (100,000)")
            # record the pointer to the raw C struct
            self._refs[i] = <void*> seq._sq
        self._refs[nseqs] = NULL

        with nogil:
            # make sure the pipeline is set to search mode and ready for a new HMM
            self._pli.mode = p7_pipemodes_e.p7_SEARCH_SEQS
            # reallocate the profile if it is too small, otherwise just clear it
            if self.profile._gm.allocM < query._hmm.M:
                libhmmer.p7_profile.p7_profile_Destroy(self.profile._gm)
                self.profile._gm = libhmmer.p7_profile.p7_profile_Create(query._hmm.M, self.alphabet._abc)
                if self.profile._gm == NULL:
                    raise AllocationError("P7_PROFILE")
            else:
                self.profile._clear()
            # configure the profile from the query HMM
            self.profile._configure(query, self.background, seq._sq.L)
            # reallocate the optimized profile if it is too small
            if self.opt._om.allocM < query._hmm.M:
                p7_oprofile.p7_oprofile_Destroy(self.opt._om)
                self.opt._om = p7_oprofile.p7_oprofile_Create(query._hmm.M, self.alphabet._abc)
                if self.opt._om == NULL:
                    raise AllocationError("P7_OPROFILE")
            # convert the profile to an optimized one
            self.opt._convert(self.profile._gm)
            # run the search loop on all database sequences while recycling memory
            Pipeline._search_loop(
                self._pli,
                self.opt._om,
                self.background._bg,
                <ESL_SQ**> self._refs,
                hits._th,
            )
            # sort hits
            hits._sort_by_key()
            hits._threshold(self)

        # return the hits
        return hits

    cpdef TopHits search_msa(
        self,
        DigitalMSA query,
        object sequences,
        Builder builder = None,
    ):
        """search_msa(self, query, sequences, builder=None)\n--

        Run the pipeline using a query alignment against a sequence database.

        Arguments:
            query (`~pyhmmer.easel.DigitalMSA`): The multiple sequence
                alignment to use to query the sequence database.
            sequences (collection of `~pyhmmer.easel.DigitalSequence`): The
                sequences to query.
            builder (`~pyhmmer.plan7.Builder`, optional): A HMM builder to
                use to convert the query to a `~pyhmmer.plan7.HMM`. If
                `None` is given, it will use a default one.

        Returns:
            `~pyhmmer.plan7.TopHits`: the hits found in the sequence database.

        Raises:
            `~pyhmmer.errors.AlphabetMismatch`: When the alphabet of the
                current pipeline does not match the alphabet of the given
                query.

        .. versionadded:: 0.3.0

        """
        assert self._pli != NULL

        cdef int                  status
        cdef int                  allocM
        cdef DigitalSequence      seq
        cdef Profile              profile
        cdef HMM                  hmm
        cdef OptimizedProfile     opt
        cdef TopHits              hits    = TopHits()

        # check the pipeline was configured with the same alphabet
        if not self.alphabet._eq(query.alphabet):
            raise AlphabetMismatch(self.alphabet, query.alphabet)
        # check we can rewind the sequences
        if not isinstance(sequences, collections.abc.Sequence):
            raise TypeError("`sequences` cannot be an iterator, expected an iterable")

        # make sure the pipeline is set to search mode and ready for a new HMM
        self._pli.mode = p7_pipemodes_e.p7_SEARCH_SEQS

        # use a default HMM builder if none was given
        builder = Builder(self.alphabet, seed=self.seed) if builder is None else builder

        # build the HMM and the profile from the query MSA
        hmm, profile, opt = builder.build_msa(query, self.background)

        # prepare an array to pass the sequences to the C function
        if self._nref < len(sequences) + 1:
            self._nref = len(sequences) + 1
            self._refs = <void**> realloc(self._refs, sizeof(void*) * (len(sequences)+1))
            if self._refs == NULL:
                raise AllocationError("void**")

        # collect sequences: we prepare an array of pointer to sequences
        # so that the C backend can iter through them without having to
        # acquire the GIL between each iteration.

        for i, seq in enumerate(sequences):
            # check alphabet
            if not self.alphabet._eq(seq.alphabet):
                raise AlphabetMismatch(self.alphabet, seq.alphabet)
            # check length
            if len(seq) > 100000:
                raise ValueError("sequence length over comparison pipeline limit (100,000)")
            # record the pointer to the raw C struct
            self._refs[i] = <void*> seq._sq
        self._refs[len(sequences)] = NULL

        with nogil:
            # run the search loop on all database sequences while recycling memory
            Pipeline._search_loop(
                self._pli,
                opt._om,
                self.background._bg,
                <ESL_SQ**> self._refs,
                hits._th,
            )
            # threshold hits
            hits._sort_by_key()
            hits._threshold(self)
        # return the hits
        return hits

    cpdef TopHits search_seq(
        self,
        DigitalSequence query,
        object sequences,
        Builder builder = None,
    ):
        """search_seq(self, query, sequences, builder=None)\n--

        Run the pipeline using a query sequence against a sequence database.

        Arguments:
            query (`~pyhmmer.easel.DigitalSequence`): The sequence object to
                use to query the sequence database.
            sequences (collection of `~pyhmmer.easel.DigitalSequence`): The
                sequences to query.
            builder (`~pyhmmer.plan7.Builder`, optional): A HMM builder to
                use to convert the query to a `~pyhmmer.plan7.HMM`. If
                `None` is given, it will use a default one.

        Returns:
            `~pyhmmer.plan7.TopHits`: the hits found in the sequence database.

        Raises:
            `~pyhmmer.errors.AlphabetMismatch`: When the alphabet of the
                current pipeline does not match the alphabet of the given
                query.

        .. versionadded:: 0.2.0

        """
        assert self._pli != NULL

        cdef int                  allocM
        cdef DigitalSequence      seq
        cdef Profile              profile
        cdef HMM                  hmm
        cdef OptimizedProfile     opt
        cdef TopHits              hits    = TopHits()
        cdef ESL_SQ**             seqs

        # check the pipeline was configure with the same alphabet
        if not self.alphabet._eq(query.alphabet):
            raise AlphabetMismatch(self.alphabet, query.alphabet)
        # check we can rewind the sequences
        if not isinstance(sequences, collections.abc.Sequence):
            raise TypeError("`sequences` cannot be an iterator, expected an iterable")

        # make sure the pipeline is set to search mode and ready for a new HMM
        self._pli.mode = p7_pipemodes_e.p7_SEARCH_SEQS

        # use a default HMM builder if none was given
        builder = Builder(self.alphabet, seed=self.seed) if builder is None else builder

        # build the HMM and the profile from the query sequence, using the first
        # as a hint for the model configuration, or reuse the profile we
        # cached if it is large enough to hold the new HMM
        hmm, profile, opt = builder.build(query, self.background)

        # prepare an array to pass the sequences to the C function
        if self._nref <= len(sequences):
            self._nref = len(sequences) + 1
            self._refs = <void**> realloc(self._refs, sizeof(void*) * (len(sequences)+1))
            if self._refs == NULL:
                raise AllocationError("void**")

        # collect sequences: we prepare an array of pointer to sequences
        # so that the C backend can iter through them without having to
        # acquire the GIL between each iteration.
        for i, seq in enumerate(sequences):
            # check alphabet
            if not self.alphabet._eq(seq.alphabet):
                raise AlphabetMismatch(self.alphabet, seq.alphabet)
            # check length
            if len(seq) > 100000:
                raise ValueError("sequence length over comparison pipeline limit (100,000)")
            # record the pointer to the raw C struct
            self._refs[i] = <void*> seq._sq
        self._refs[len(sequences)] = NULL

        with nogil:
            # run the search loop on all database sequences while recycling memory
            Pipeline._search_loop(
                self._pli,
                opt._om,
                self.background._bg,
                <ESL_SQ**> self._refs,
                hits._th,
            )
            # threshold hits
            hits._sort_by_key()
            hits._threshold(self)
        # return the hits
        return hits

    @staticmethod
    cdef int _search_loop(
        P7_PIPELINE* pli,
        P7_OPROFILE* om,
        P7_BG*       bg,
        ESL_SQ**     sq,
        P7_TOPHITS*  th,
    ) nogil except 1:
        cdef int status

        # configure the pipeline for the current HMM
        status = libhmmer.p7_pipeline.p7_pli_NewModel(pli, om, bg)
        if status == libeasel.eslEINVAL:
            raise ValueError("model does not have bit score thresholds expected by the pipeline")
        elif status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_pli_NewModel")

        # run the inner loop on all sequences
        while sq[0] != NULL:

            # configure the profile, background and pipeline for the new sequence
            status = libhmmer.p7_pipeline.p7_pli_NewSeq(pli, sq[0])
            if status != libeasel.eslOK:
                raise UnexpectedError(status, "p7_pli_NewSeq")
            status = libhmmer.p7_bg.p7_bg_SetLength(bg, sq[0].n)
            if status != libeasel.eslOK:
                raise UnexpectedError(status, "p7_bg_SetLength")
            status = p7_oprofile.p7_oprofile_ReconfigLength(om, sq[0].n)
            if status != libeasel.eslOK:
                raise UnexpectedError(status, "p7_oprofile_ReconfigLength")

            # run the pipeline on the sequence
            status = libhmmer.p7_pipeline.p7_Pipeline(pli, om, bg, sq[0], NULL, th)
            if status == libeasel.eslEINVAL:
                raise ValueError("model does not have bit score thresholds expected by the pipeline")
            elif status == libeasel.eslERANGE:
                raise OverflowError("numerical overflow in the optimized vector implementation")
            elif status != libeasel.eslOK:
                raise UnexpectedError(status, "p7_Pipeline")

            # clear pipeline for reuse for next target
            libhmmer.p7_pipeline.p7_pipeline_Reuse(pli)

            # advance to next sequence
            sq += 1 #&(sq[1])

        # Return 0 to indicate success
        return 0

    cpdef TopHits scan_seq(
        self,
        DigitalSequence query,
        object hmms,
    ):
        """scan_seq(self, query, hmms)\n--

        Run the pipeline using a query sequence against a profile database.

        Arguments:
            query (`~pyhmmer.easel.DigitalSequence`): The sequence object to
                use to query the profile database.
            hmms (iterable of `~pyhmmer.easel.DigitalSequence`): The
                HMM profiles to query. Pass a `~pyhmmer.easel.HMMFile`
                instance to read from disk iteratively.

        Returns:
            `~pyhmmer.plan7.TopHits`: the hits found in the profile database.

        Raises:
            `~pyhmmer.errors.AlphabetMismatch`: When the alphabet of the
                current pipeline does not match the alphabet of the given
                query or profile.

        Caution:
            In the current version, this method is not optimized to use
            the *pressed* database, even if it exists. This will cause the
            MSV and SSV filters to be rebuilt at each iteration, which could
            be slow. Consider at least pre-fetching the HMM database if
            calling this method several times in a row.

        .. versionadded:: v0.4.0

        """
        cdef int                  allocM
        cdef Profile              profile
        cdef HMM                  hmm
        cdef OptimizedProfile     opt
        cdef TopHits              hits    = TopHits()

        assert self._pli != NULL

        # check the pipeline was configure with the same alphabet
        if not self.alphabet._eq(query.alphabet):
            raise AlphabetMismatch(self.alphabet, query.alphabet)

        # make sure the pipeline is set to scan mode
        self._pli.mode = p7_pipemodes_e.p7_SCAN_MODELS

        # get an iterator over the input hmms, with an early return if
        # no sequences were given, and extract the raw pointer to the
        # HMM from the Python object
        hmms_iter = iter(hmms)
        hmm = next(hmms_iter, None)
        if hmm is None:
            return hits
        if not self.alphabet._eq(hmm.alphabet):
            raise AlphabetMismatch(self.alphabet, hmm.alphabet)

        # run the search loop on all database sequences while recycling memory
        self._scan_loop(
            self._pli,
            query._sq,
            self.background._bg,
            hmm._hmm,
            hits._th,
            hmms_iter,
            hmm.alphabet
        )

        # threshold hits
        hits._sort_by_key()
        hits._threshold(self)

        # return the hits
        return hits

    cdef int _scan_loop(
        self,
        P7_PIPELINE* pli,
        ESL_SQ*      sq,
        P7_BG*       bg,
        P7_HMM*      hm,
        P7_TOPHITS*  th,
        object       hmm_iter,
        Alphabet     hmm_alphabet
    ) except 1:
        cdef int              status
        cdef Alphabet         self_alphabet = self.alphabet
        cdef Profile          profile       = Profile(hm.M, hmm_alphabet)
        cdef OptimizedProfile optimized     = OptimizedProfile(hm.M, hmm_alphabet)
        cdef P7_PROFILE*      gm            = profile._gm
        cdef P7_OPROFILE*     om            = optimized._om
        cdef HMM              hmm

        with nogil:
            # verify the alphabet
            if not self_alphabet._eq(hmm_alphabet):
                raise AlphabetMismatch(self_alphabet, hmm_alphabet)
            # verify the length
            if sq.n > 100000:
                raise ValueError("sequence length over comparison pipeline limit (100,000)")

            # configure the pipeline for the current sequence
            status = libhmmer.p7_pipeline.p7_pli_NewSeq(pli, sq)
            if status != libeasel.eslOK:
                raise UnexpectedError(status, "p7_pli_NewSeq")

            # run the inner loop on all HMMs
            while hm != NULL:
                # reallocate if too small
                if gm.M < hm.M:
                    with gil:
                        profile = Profile(hm.M, hmm_alphabet)
                        optimized = OptimizedProfile(hm.M, hmm_alphabet)
                        gm, om = profile._gm, optimized._om

                # configure the new profile
                status = libhmmer.modelconfig.p7_ProfileConfig(hm, bg, gm, sq.L, p7_LOCAL)
                if status != libeasel.eslOK:
                    raise UnexpectedError(status, "p7_ProfileConfig")

                # convert the new profile to an optimized one
                status = p7_oprofile.p7_oprofile_Convert(gm, om)
                if status != libeasel.eslOK:
                    raise UnexpectedError(status, "p7_oprofile_Convert")

                # configure the profile, background and pipeline for the new HMM
                status = libhmmer.p7_pipeline.p7_pli_NewModel(pli, om, bg)
                if status != libeasel.eslOK:
                    raise UnexpectedError(status, "p7_pli_NewModel")
                status = libhmmer.p7_bg.p7_bg_SetLength(bg, sq.n)
                if status != libeasel.eslOK:
                    raise UnexpectedError(status, "p7_bg_SetLength")
                status = p7_oprofile.p7_oprofile_ReconfigLength(om, sq.n)
                if status != libeasel.eslOK:
                    raise UnexpectedError(status, "p7_oprofile_ReconfigLength")

                # run the pipeline on the sequence
                status = libhmmer.p7_pipeline.p7_Pipeline(pli, om, bg, sq, NULL, th)
                if status == libeasel.eslEINVAL:
                    raise ValueError("model does not have bit score thresholds expected by the pipeline")
                elif status == libeasel.eslERANGE:
                    raise OverflowError("numerical overflow in the optimized vector implementation")
                elif status != libeasel.eslOK:
                    raise UnexpectedError(status, "p7_Pipeline")

                # clear pipeline for reuse for next target
                libhmmer.p7_pipeline.p7_pipeline_Reuse(pli)

                # acquire the GIL just long enough to get the next HMM
                with gil:
                    hmm = next(hmm_iter, None)
                    if hmm is None:
                        hm = NULL
                    else:
                        hm = hmm._hmm
                        hmm_alphabet = hmm.alphabet
                        # verify the alphabet
                        if not self_alphabet._eq(hmm_alphabet):
                            raise AlphabetMismatch(self_alphabet, hmm_alphabet)

        # Return 0 to indicate success
        return 0


cdef class Profile:
    """A Plan7 search profile.

    Attributes:
        alphabet (`~pyhmmer.easel.Alphabet`): The alphabet for which this
            profile was configured.
        offsets (`~pyhmmer.plan7._Offsets`): The disk offsets for this
            profile, if it was loaded from a pressed HMM file.
        evalue_parameters (`~pyhmmer.plan7._EvalueParameters`): The e-value
            parameters for this profile.
        cutoffs (`~pyhmmer.plan7._Cutoffs`): The Pfam score cutoffs for this
            profile, if any.

    .. versionchanged:: 0.4.6
       Added the `~Profile.evalue_parameters` and `~Profile.cutoffs` attributes.

    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._gm = NULL
        self.alphabet = None
        self.offsets = None
        self.evalue_parameters = None
        self.cutoffs = None

    def __init__(self, int M, Alphabet alphabet):
        """__init__(self, M, alphabet)\n--

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
            raise AllocationError("P7_PROFILE")
        # expose the disk offsets as the `offsets` attribute
        self.offsets = _Offsets.__new__(_Offsets)
        self.offsets._offs = &self._gm.offs
        self.offsets._owner = self
        # expose the e-value parameters as the `evalue_parameters` attribute
        self.evalue_parameters = _EvalueParameters.__new__(_EvalueParameters)
        self.evalue_parameters._evparams = &self._gm.evparam
        self.evalue_parameters._owner = self
        # expose the score cutoffs as the `cutoffs` attribute
        self.cutoffs = _Cutoffs.__new__(_Cutoffs)
        self.cutoffs._cutoffs = &self._gm.cutoff
        self.cutoffs._flags = NULL
        self.cutoffs._is_profile = True

    def __dealloc__(self):
        libhmmer.p7_profile.p7_profile_Destroy(self._gm)

    def __copy__(self):
        return self.copy()

    def __eq__(self, object other):
        assert self._gm != NULL

        if not isinstance(other, Profile):
            return False

        cdef Profile      p = <Profile> other
        cdef int     status = libhmmer.p7_profile.p7_profile_Compare(self._gm, p._gm, 0.0)

        if status == libeasel.eslOK:
            return True
        elif status == libeasel.eslFAIL:
            return False
        else:
            raise UnexpectedError(status, "p7_profile_Compare")

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

    # --- Methods ------------------------------------------------------------

    cdef int _clear(self) nogil except 1:
        cdef int status
        status = libhmmer.p7_profile.p7_profile_Reuse(self._gm)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_profile_Reuse")
        return 0

    def clear(self):
        """clear(self)\n--

        Clear internal buffers to reuse the profile without reallocation.

        """
        assert self._gm != NULL
        self._clear()

    cdef int _configure(self, HMM hmm, Background background, int L, bint multihit=True, bint local=True) nogil except 1:
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
        status = libhmmer.modelconfig.p7_ProfileConfig(hm, bg, gm, L, mode)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_ProfileConfig")

        return 0

    def configure(self, HMM hmm, Background background, int L, bint multihit=True, bint local=True):
        """configure(self, hmm, background, L, multihit=True, local=True)\n--

        Configure a search profile using the given models.

        Arguments:
            hmm (`pyhmmer.plan7.HMM`): The model HMM with core probabilities.
            bg (`pyhmmer.plan7.Background`): The null background model.
            L (`int`): The expected target sequence length.
            multihit (`bool`): Whether or not to use multihit modes.
            local (`bool`): Whether or not to use non-local modes.

        """
        assert self._gm != NULL
        with nogil:
            self._configure(hmm, background, L, multihit, local)

    cpdef Profile copy(self):
        """copy(self)\n--

        Return a copy of the profile with the exact same configuration.

        """
        assert self._gm != NULL

        cdef int status
        cdef Profile new = Profile.__new__(Profile)
        new.alphabet = self.alphabet
        new.offsets = _Offsets(new)
        new.evalue_parameters = _EvalueParameters(new)
        new.cutoffs = _Cutoffs(new)

        # allocate a new profile
        with nogil:
            new._gm = libhmmer.p7_profile.p7_profile_Create(self._gm.allocM, self.alphabet._abc)
        if not new._gm:
            raise AllocationError("P7_PROFILE")

        # copy the current profile to the new profile
        with nogil:
            status = libhmmer.p7_profile.p7_profile_Copy(self._gm, new._gm)
        if status == libeasel.eslOK:
            return new
        else:
            raise UnexpectedError(status, "p7_profile_Copy")

    cpdef bint is_local(self):
        """is_local(self)\n--

        Return whether or not the profile is in a local alignment mode.

        """
        assert self._gm != NULL
        return libhmmer.p7_profile.p7_profile_IsLocal(self._gm)

    cpdef bint is_multihit(self):
        """is_multihit(self)\n--

        Returns whether or not the profile is in a multihit alignment mode.

        """
        assert self._gm != NULL
        return libhmmer.p7_profile.p7_profile_IsMultihit(self._gm)

    cpdef OptimizedProfile optimized(self):
        """optimized(self)\n--

        Convert the profile to a platform-specific optimized profile.

        Returns:
            `OptimizedProfile`: The platform-specific optimized profile built
            using the configuration of this profile.

        """
        assert self._gm != NULL

        cdef int              status
        cdef OptimizedProfile opt

        opt = OptimizedProfile(self._gm.M, self.alphabet)
        with nogil:
            opt._convert(self._gm)

        return opt


cdef class TopHits:
    """A ranked list of top-scoring hits.

    `TopHits` are thresholded using the parameters from the pipeline, and are
    sorted by key when you obtain them from a `Pipeline` instance::

        >>> abc = thioesterase.alphabet
        >>> hits = Pipeline(abc).search_hmm(thioesterase, proteins)
        >>> hits.is_sorted()
        True

    Use `len` to query the number of top hits, and the usual indexing notation
    to extract a particular `Hit`::

        >>> len(hits)
        1
        >>> hits[0].name
        b'938293.PRJEB85.HG003687_113'

    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._th = NULL

    def __init__(self):
        """__init__(self)\n--

        Create an empty `TopHits` instance.

        """
        if self._th != NULL:
            libhmmer.p7_tophits.p7_tophits_Destroy(self._th)
        with nogil:
            self._th = libhmmer.p7_tophits.p7_tophits_Create()
        if self._th == NULL:
            raise AllocationError("P7_TOPHITS")

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

    # --- Properties ---------------------------------------------------------

    @property
    def reported(self):
        """`int`: The number of hits that are above the reporting threshold.
        """
        return self._th.nreported

    @property
    def included(self):
        """`int`: The number of hits that are above the inclusion threshold.
        """
        return self._th.nincluded

    # --- Methods ------------------------------------------------------------

    cdef int _threshold(self, Pipeline pipeline) nogil except 1:
        cdef int status = libhmmer.p7_tophits.p7_tophits_Threshold(self._th, pipeline._pli)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_tophits_Threshold")
        self.Z = pipeline._pli.Z
        self.domZ = pipeline._pli.domZ
        self.long_targets = pipeline._pli.long_targets
        return 0

    def threshold(self, Pipeline pipeline):
        """threshold(self, pipeline)\n--

        Apply score and e-value thresholds using pipeline parameters.

        This function is automatically called in `Pipeline.search_hmm` or
        `Pipeline.search_seq`, so it should have limited usefulness at the
        Python level.

        """
        assert self._th != NULL
        assert pipeline._pli != NULL
        with nogil:
            self._threshold(pipeline)

    cdef int _sort_by_key(self) nogil except 1:
        cdef int status = libhmmer.p7_tophits.p7_tophits_SortBySortkey(self._th)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_tophits_SortBySortkey")
        return 0

    cdef int _sort_by_seqidx(self) nogil except 1:
        cdef int status = libhmmer.p7_tophits.p7_tophits_SortBySeqidxAndAlipos(self._th)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_tophits_SortBySeqidxAndAlipos")
        return 0

    def sort(self, str by="key"):
        """sort(self, by="key")\n--

        Sort hits in the current instance using the given method.

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
            raise ValueError("Invalid value for `by` argument: {!r}".format(by))

    def is_sorted(self, str by="key"):
        """is_sorted(self, by="key")\n--

        Check whether or not the hits are sorted with the given method.

        See `~pyhmmer.plan7.TopHits.sort` for a list of allowed values for
        the ``by`` argument.

        """
        assert self._th != NULL
        if by == "key":
            return self._th.is_sorted_by_sortkey
        elif by == "seqidx":
            return self._th.is_sorted_by_seqidx
        else:
            raise ValueError("Invalid value for `by` argument: {!r}".format(by))

    cpdef MSA to_msa(self, Alphabet alphabet, bint trim=False, bint digitize=False, bint all_consensus_cols=False):
        """to_msa(self, alphabet, trim=False, digitize=False, all_consensus_cols=False)\n--

        Create multiple alignment of all included domains.

        Arguments:
            alphabet (`~pyhmmer.easel.Alphabet`): The alphabet of the
                HMM this `TopHits` was obtained from. It is required to
                convert back hits to single sequences.
            trim (`bool`): Trim off any residues that get assigned to
                flanking :math:`N` and :math:`C` states (in profile traces)
                or :math:`I_0` and :math:`I_m` (in core traces).
            digitize (`bool`): If set to `True`, returns a `DigitalMSA`
                instead of a `TextMSA`.
            all_consensus_cols (`bool`): Force a column to be created for
                every consensus column in the model, even if it means having
                all gap character in a column.

        Returns:
            `~pyhmmer.easel.MSA`: A multiple sequence alignment containing
            the reported hits, either a `TextMSA` or a `DigitalMSA`
            depending on the value of the ``digitize`` argument.

        .. versionadded:: 0.3.0

        """
        assert self._th != NULL
        assert alphabet._abc != NULL

        cdef int status
        cdef int flags  = libhmmer.p7_DEFAULT
        cdef MSA msa

        if trim:
            flags |= libhmmer.p7_TRIM
        if all_consensus_cols:
            flags |= libhmmer.p7_ALL_CONSENSUS_COLS
        if digitize:
            flags |= libhmmer.p7_DIGITIZE
            msa = DigitalMSA.__new__(DigitalMSA, alphabet)
        else:
            msa = TextMSA.__new__(TextMSA)

        status = libhmmer.p7_tophits.p7_tophits_Alignment(
            self._th,
            alphabet._abc,
            NULL, # inc_sqarr
            NULL, # inc_trarr
            0,    # inc_n
            flags,
            &msa._msa
        )

        if status == libeasel.eslOK:
            return msa
        elif status == libeasel.eslFAIL:
            raise ValueError("No included domains found")
        else:
            raise UnexpectedError(status, "p7_tophits_Alignment")


# --- Module init code -------------------------------------------------------

impl_Init()
p7_FLogsumInit()

include "exceptions.pxi"
