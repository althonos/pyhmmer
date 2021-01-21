# coding: utf-8
# cython: language_level=3, linetrace=True
"""High-level interface to the Plan7 data model.

Plan7 is the model architecture used by HMMER since HMMER2.

See Also:
    Details about the Plan 7 architecture in the `HMMER documentation
    <http://www.csb.yale.edu/userguides/seq/hmmer/docs/node11.html>`_.

"""

# --- C imports --------------------------------------------------------------

from cpython.ref cimport PyObject
from cpython.exc cimport PyErr_Clear
from libc.stdlib cimport malloc, realloc, free
from libc.stdint cimport uint32_t
from libc.stdio cimport fprintf, FILE, stdout, fclose
from libc.string cimport strdup
from libc.math cimport exp, ceil

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
from libhmmer cimport p7_LOCAL, p7_offsets_e
from libhmmer.logsum cimport p7_FLogsumInit
from libhmmer.p7_builder cimport P7_BUILDER, p7_archchoice_e, p7_wgtchoice_e, p7_effnchoice_e
from libhmmer.p7_hmmfile cimport p7_hmmfile_formats_e
from libhmmer.p7_tophits cimport p7_hitflags_e
from libhmmer.p7_alidisplay cimport P7_ALIDISPLAY
from libhmmer.p7_pipeline cimport P7_PIPELINE, p7_pipemodes_e, p7_zsetby_e
from libhmmer.p7_profile cimport p7_LOCAL, p7_GLOCAL, p7_UNILOCAL, p7_UNIGLOCAL

IF HMMER_IMPL == "VMX":
    from libhmmer.impl_vmx cimport p7_oprofile, p7_omx, impl_Init
    from libhmmer.impl_vmx.p7_oprofile cimport P7_OPROFILE
    from libhmmer.impl_vmx.io cimport p7_oprofile_Write
ELIF HMMER_IMPL == "SSE":
    from libhmmer.impl_sse cimport p7_oprofile, p7_omx, impl_Init
    from libhmmer.impl_sse.p7_oprofile cimport P7_OPROFILE
    from libhmmer.impl_sse.io cimport p7_oprofile_Write

from .easel cimport Alphabet, Sequence, DigitalSequence
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

import collections.abc
import errno
import io
import os
import sys
import warnings

from .errors import AllocationError, UnexpectedError


# --- Cython classes ---------------------------------------------------------


cdef class Alignment:
    """A single alignment of a sequence to a profile.
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
    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._bg = NULL
        self.alphabet = None
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
        self.alphabet = alphabet
        self.uniform = uniform
        with nogil:
            if uniform:
                self._bg = libhmmer.p7_bg.p7_bg_CreateUniform(alphabet._abc)
            else:
                self._bg = libhmmer.p7_bg.p7_bg_Create(alphabet._abc)
        if self._bg == NULL:
            raise AllocationError("P7_BG")

    def __dealloc__(self):
        libhmmer.p7_bg.p7_bg_Destroy(self._bg)

    def __copy__(self):
        return self.copy()

    # --- Magic methods ------------------------------------------------------

    @property
    def L(self):
        """`int`: The mean of the null model length distribution, in residues.
        """
        assert self._bg != NULL
        return <int> ceil(self._bg.p1 / (1 - self._bg.p1))

    @L.setter
    def L(self, int L):
        cdef int    status
        cdef P7_BG* bg     = self._bg

        assert self._bg != NULL

        with nogil:
            status = libhmmer.p7_bg.p7_bg_SetLength(bg, L)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_bg_SetLength")

    # --- Methods ------------------------------------------------------------

    cpdef Background copy(self):
        """copy(self)\n--

        Create a copy of the null model with the same parameters.

        """
        cdef Background new = Background.__new__(Background)
        new.alphabet = self.alphabet
        new.uniform = self.uniform
        with nogil:
            new._bg = libhmmer.p7_bg.p7_bg_Clone(self._bg)
        if new._bg == NULL:
            raise AllocationError("P7_BG")
        return new


cdef class Builder:
    """A factory for constructing new HMMs from raw sequences.
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
        str effective_number="entropy",
        str prior_scheme="alphabet",
        float symfrac=0.5,
        float fragthresh=0.5,
        float wid=0.62,
        float esigma=45.0,
        float eid=0.62,
        int EmL=200,
        int EmN=200,
        int EvL=200,
        int EvN=200,
        int EfL=100,
        int EfN=200,
        float Eft=0.04,
        int seed=42,
        object ere=None,
    ):
        """__init__(self, alphabet, *, architecture="fast", weighting="pb", effective_number="entropy", prior_scheme="alpha", symfrac=0.5, fragthresh=0.5, wid=0.62, esigma=45.0, eid=0.62, EmL=200, EmN=200, EvL=200, EvN=200, EfL=100, EfN=200, Eft=0.04, seed=42, ere=None)\n--

        Create a new sequence builder with the given configuration.

        Arguments:
            alphabet (`~pyhmmer.easel.Alphabet`): The alphabet the builder
                expects the sequences to be in.

        """
        # extract alphabet and create raw P7_BUILDER object
        self.alphabet = alphabet
        with nogil:
            self._bld = libhmmer.p7_builder.p7_builder_Create(NULL, alphabet._abc)
        if self._bld == NULL:
            raise AllocationError("P7_BG")

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

        # reset the random number generator
        self.seed = seed

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
        else:
            _effn = self._EFFECTIVE_STRATEGY.get(effective_number)
            if _effn is not None:
                self._bld.effn_strategy = _effn
            else:
                raise ValueError(f"Invalid value for 'effective_number': {effective_number}")

        # set the RE target if given one
        if ere is not None:
            self._bld.re_target = ere

        # set the prior scheme
        self.prior_scheme = prior_scheme
        with nogil:
            if prior_scheme is None:
                self._bld.prior = NULL
            elif prior_scheme == "laplace":
                self._bld.prior = libhmmer.p7_prior.p7_prior_CreateLaplace(self.alphabet._abc)
            elif prior_scheme != "alphabet":
                raise ValueError("Invalid value for 'prior_scheme': {prior_scheme}")

    def __dealloc__(self):
        libhmmer.p7_builder.p7_builder_Destroy(self._bld)

    def __copy__(self):
        return self.copy()

    # --- Properties ---------------------------------------------------------

    @property
    def seed(self):
        """`int`: The seed used by the internal random number generator.

        Setting the seed will effectively reinitialize the internal RNG. In
        the special case the seed is *0*, a one-time arbitrary seed will be
        chosen and the RNG will no be reseeded for reproducibility.
        """
        return libeasel.random.esl_randomness_GetSeed(self._bld.r)

    @seed.setter
    def seed(self, uint32_t seed):
        with nogil:
            status = libeasel.random.esl_randomness_Init(self._bld.r, seed)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_randomness_Init")
        self._bld.do_reseeding = seed != 0

    # --- Methods ------------------------------------------------------------

    cpdef tuple build(
        self,
        DigitalSequence sequence,
        Background background,
        float popen = 0.02,
        float pextend = 0.4
    ):
        """build(self, sequence, background, popen=0.02, pextend=0.4)\n--

        Build a new HMM from ``sequence`` using the builder configuration.

        Arguments:
            sequence (`~pyhmmer.easel.DigitalSequence`): A single biological
                sequence in digital mode to build a HMM with.
            background (`pyhmmer.plan7.background`): The background model to
                use to create the HMM.
            popen (`float`): The *gap open* probability to use with the score
                system. Defaults to *0.02*.
            pextend (`float`): The *gap extend* probability to use with the
                score system. Defaults to *0.4*.

        Raises:
            `ValueError`: When either ``sequence`` or ``background`` have the
                wrong alphabet for this builder.

        """
        cdef int              status
        cdef HMM              hmm     = HMM.__new__(HMM)
        cdef Profile          profile = Profile.__new__(Profile)
        cdef OptimizedProfile opti    = OptimizedProfile.__new__(OptimizedProfile)

        # check alphabet and assign it to newly created objects
        hmm.alphabet = profile.alphabet = opti.alphabet = self.alphabet
        if background.alphabet != self.alphabet:
            raise ValueError("background does not have the right alphabet")

        # load score matrix
        # TODO: allow changing from the default scoring matrix
        # TODO: allow caching the parameter values to avoid resetting
        #       everytime `build` is called.
        with nogil:
            status = libhmmer.p7_builder.p7_builder_SetScoreSystem(
                self._bld,
                NULL, # --mxfile
                NULL, # env
                popen, # popen
                pextend,  # pextend
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
            if status != libeasel.eslOK:
                raise UnexpectedError(status, "p7_SingleBuilder")
        return (hmm, profile, opti)

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
            EfL=self._bld.EfN,
            Eft=self._bld.Eft,
            seed=self.seed,
            ere=self._bld.re_target,
        )


cdef class Domain:
    """A single domain in a query `~pyhmmer.plan7.Hit`.
    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self, Hit hit, size_t index):
        self._dom = &hit._hit.dcl[index]
        self.hit = hit
        self.alignment = Alignment(self)

    # --- Properties ---------------------------------------------------------

    @property
    def env_from(self):
        assert self._dom != NULL
        return self._dom.ienv

    @property
    def env_to(self):
        assert self._dom != NULL
        return self._dom.jenv

    @property
    def score(self):
        """`float`: The overall score in bits, null-corrected.
        """
        assert self._dom != NULL
        return self._dom.bitscore

    @property
    def bias(self):
        assert self._dom != NULL
        return self._dom.dombias * libeasel.eslCONST_LOG2R

    @property
    def correction(self):
        assert self._dom != NULL
        return self._dom.domcorrection * libeasel.eslCONST_LOG2R

    @property
    def envelope_score(self):
        assert self._dom != NULL
        return self._dom.envsc * libeasel.eslCONST_LOG2R

    @property
    def c_evalue(self):
        """`float`: The independent e-value for the domain.
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


cdef class Domains:
    """A sequence of domains corresponding to a single `~pyhmmer.plan7.Hit`.
    """

    def __cinit__(self, Hit hit):
        self.hit = hit

    def __len__(self):
        return self.hit._hit.ndom

    def __getitem__(self, index):
        if index < 0:
            index += self.hit._hit.ndom
        if index >= self.hit._hit.ndom or index < 0:
            raise IndexError("list index out of range")
        return Domain(self.hit, <size_t> index)


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
    def bias(self):
        assert self._hit != NULL
        return self._hit.pre_score - self._hit.score

    @property
    def domains(self):
        """`Domains`: The list of domains aligned to this hit.
        """
        return Domains(self)

    @property
    def evalue(self):
        """`float`: The e-value of the hit.
        """
        assert self._hit != NULL
        if self.hits.long_targets:
            return exp(self._hit.lnP)
        else:
            return exp(self._hit.lnP) * self.hits.Z

    # --- Methods ------------------------------------------------------------

    cpdef bint is_included(self):
        return self._hit.flags & p7_hitflags_e.p7_IS_INCLUDED != 0

    cpdef bint is_reported(self):
        return self._hit.flags & p7_hitflags_e.p7_IS_REPORTED != 0

    cpdef bint is_new(self):
        return self._hit.flags & p7_hitflags_e.p7_IS_NEW != 0

    cpdef bint is_dropped(self):
        return self._hit.flags & p7_hitflags_e.p7_IS_DROPPED != 0

    cpdef bint is_duplicate(self):
        return self._hit.flags & p7_hitflags_e.p7_IS_DUPLICATE


cdef class HMM:
    """A data structure storing the Plan7 Hidden Markov Model.
    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self.alphabet = None
        self._hmm = NULL

    def __init__(self, int M, Alphabet alphabet):
        """__init__(self, M, alphabet)\n--

        Create a new HMM from scratch.

        Arguments:
            M (`int`): The length of the model (i.e. the number of nodes).
            alphabet (`~pyhmmer.easel.Alphabet`): The alphabet of the model.

        """
        self.alphabet = alphabet
        with nogil:
            self._hmm = libhmmer.p7_hmm.p7_hmm_Create(M, alphabet._abc)
        if not self._hmm:
            raise AllocationError("P7_HMM")

    def __dealloc__(self):
        libhmmer.p7_hmm.p7_hmm_Destroy(self._hmm)

    # --- Properties ---------------------------------------------------------

    @property
    def M(self):
        """`int`: The length of the model (i.e. the number of nodes).
        """
        return self._hmm.M

    @property
    def name(self):
        """`bytes` or `None`: The name of the HMM, if any.
        """
        return None if self._hmm.name == NULL else <bytes> self._hmm.name

    @name.setter
    def name(self, bytes name):
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
        return None if self._hmm.acc == NULL else <bytes> self._hmm.acc

    @accession.setter
    def accession(self, bytes accession):
        cdef char* acc = NULL if accession is None else <char*> accession
        cdef int err = libhmmer.p7_hmm.p7_hmm_SetAccession(self._hmm, acc)
        if err == libeasel.eslEMEM:
            raise AllocationError("char*")
        elif err != libeasel.eslOK:
            raise UnexpectedError(err, "p7_hmm_SetAccession")

    @property
    def description(self):
        """`bytes` or `None`: The description of the HMM, if any.
        """
        return None if self._hmm.desc == NULL else <bytes> self._hmm.desc

    @description.setter
    def description(self, bytes description):
        cdef char* desc = NULL if description is None else <char*> description
        cdef int err = libhmmer.p7_hmm.p7_hmm_SetDescription(self._hmm, desc)
        if err == libeasel.eslEMEM:
            raise AllocationError("char*")
        elif err != libeasel.eslOK:
            raise UnexpectedError(err, "p7_hmm_SetDescription")

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
        with nogil:
            libhmmer.p7_hmm.p7_hmm_Zero(self._hmm)


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
            return py_hmm
        elif status == libeasel.eslEOF:
            raise StopIteration()
        elif status == libeasel.eslEMEM:
            raise AllocationError("P7_HMM")
        elif status == libeasel.eslESYS:
            raise OSError(self._hfp.errbuf.decode('ascii'))
        elif status == libeasel.eslEFORMAT:
            raise ValueError("Invalid format in file: {}".format(self._hfp.errbuf.decode('ascii')))
        elif status == libeasel.eslEINCOMPAT:
            alphabet = libeasel.alphabet.esl_abc_DecodeType(self._alphabet.type)
            raise ValueError("HMM is not in the expected {} alphabet".format(alphabet))
        else:
            raise UnexpectedError(status, "p7_hmmfile_Read")


    # --- Methods ------------------------------------------------------------

    cpdef void close(self):
        """close(self)\n--

        Close the HMM file and free resources.

        This method has no effect if the file is already closed. It called
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
            format = token[5:].decode("ascii", errors="replace")
            if format in self._FORMATS:
                hfp.format = self._FORMATS[format]
            else:
                hfp.parser = NULL
        elif token.startswith(b"HMMER2.0"):
            hfp.parser = read_asc20hmm
            hfp.format = p7_hmmfile_formats_e.p7_HMMFILE_20

        # check the format tag was recognized
        if hfp.parser == NULL:
            text = token.decode(encoding='ascii', errors='replace')
            libhmmer.p7_hmmfile.p7_hmmfile_Close(hfp)
            raise ValueError("Unrecognized format tag in HMM file: {!r}".format(text))

        # return the finalized P7_HMMFILE*
        return hfp


cdef class OptimizedProfile:
    """An optimized profile that uses platform-specific instructions.
    """

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
        self.alphabet = alphabet
        with nogil:
            self._om = p7_oprofile.p7_oprofile_Create(M, alphabet._abc)
        if self._om == NULL:
            raise AllocationError("P7_OPROFILE")

    def __dealloc__(self):
        p7_oprofile.p7_oprofile_Destroy(self._om)

    def __copy__(self):
        return self.copy()

    @property
    def offsets(self):
        return _Offsets.__new__(_Offsets, self)

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


cdef class _Offsets:
    """A view over the disk offsets of an optimized profile.
    """

    def __cinit__(self, OptimizedProfile opt):
        self.opt = opt
        self._offs = &opt._om.offs

    def __copy__(self):
        return _Offsets.__new__(_Offsets, self.om)

    def __str__(self):
        ty = type(self).__name__
        return "<offsets of {!r} model={!r} filter={!r} profile={!r}>".format(
            self.opt,
            self.model,
            self.filter,
            self.profile,
        )

    @property
    def model(self):
        cdef off_t model = self._offs[0][<int> p7_offsets_e.p7_MOFFSET]
        return None if model == -1 else model

    @model.setter
    def model(self, object model):
        self._offs[0][<int> p7_offsets_e.p7_MOFFSET] = -1 if model is None else model

    @property
    def filter(self):
        cdef off_t filter = self._offs[0][<int> p7_offsets_e.p7_FOFFSET]
        return None if filter == -1 else filter

    @filter.setter
    def filter(self, object filter):
        self._offs[0][<int> p7_offsets_e.p7_FOFFSET] = -1 if filter is None else filter

    @property
    def profile(self):
        cdef off_t profile = self._offs[0][<int> p7_offsets_e.p7_POFFSET]
        return None if profile == -1 else profile

    @profile.setter
    def profile(self, object profile):
        self._offs[0][<int> p7_offsets_e.p7_MOFFSET] = -1 if profile is None else profile


cdef class Pipeline:
    """An HMMER3 accelerated sequence/profile comparison pipeline.
    """

    M_HINT = 100         # default model size
    L_HINT = 100         # default sequence size
    LONG_TARGETS = False

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._pli = NULL
        self.profile = None
        self.background = None

    def __init__(
        self,
        Alphabet alphabet,
        Background background = None,
        *,
        bint bias_filter=True,
        float report_e=10.0,
        bint null2=True,
        uint32_t seed=42,
        object Z=None,
        object domZ=None,
    ):
        """__init__(self, alphabet, background=None, *, bias_filter=True, report_e=10.0, null2=True, seed=42, Z=None, domZ=None)\n--

        Instantiate and configure a new accelerated comparison pipeline.

        Arguments:
            alphabet (`~pyhmmer.easel.Alphabet`): The biological alphabet the
                of the HMMs and sequences that are going to be compared. Used
                to build the background model.
            background (`~pyhmmer.plan7.Background`, optional): The background
                model to use with the pipeline, or ``None`` to create and use
                a default one. *The pipeline needs ownership of the background
                model, so make sure to use* `Background.copy` *if passing a
                custom background model here.*

        Keyword Arguments:
            bias_filter (`bool`): Whether or not to enable composition bias
                filter. Defaults to ``True``.
            null2 (`bool`): Whether or not to compute biased composition score
                corrections. Defaults to ``True``.
            report_e (`float`): Report hits with e-value lower than or
                equal to this threshold in output. Defaults to ``10.0``.
            seed (`int`, optional): The seed to use with the random number
                generator. Pass *0* to use a one-time arbitrary seed, or
                `None` to keep the default seed from HMMER.

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

        # configure the pipeline with the additional keyword arguments
        self.null2 = null2
        self.bias_filter = bias_filter
        self.report_e = report_e
        self.Z = Z
        self.domZ = domZ
        self.seed = seed

    def __dealloc__(self):
        libhmmer.p7_pipeline.p7_pipeline_Destroy(self._pli)

    # --- Properties ---------------------------------------------------------

    @property
    def Z(self):
        """`float` or `None`: The number of effective targets searched.

        It is used to compute the independent e-value for each domain, and
        for an entire hit. If `None`, the parameter number will be set
        automatically after all the comparisons have been done. Otherwise, it
        can be set to an arbitrary number.

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
        `None`, the parameter number will be set automatically after all the
        comparisons have been done, and all hits have been thresholded.
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
        """`int`: The seed used by the internal random number generator.
        """
        return libeasel.random.esl_randomness_GetSeed(self._pli.r)

    @seed.setter
    def seed(self, uint32_t seed):
        cdef int status
        with nogil:
            status = libeasel.random.esl_randomness_Init(self._pli.r, seed)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_randomness_Init")
        self._pli.do_reseeding = self._pli.ddef.do_reseeding = seed != 0

    # --- Methods ------------------------------------------------------------

    cpdef void clear(self):
        """clear(self)\n--

        Reset the pipeline configuration to its default state.

        """
        cdef int      status
        cdef uint32_t seed   = self.seed

        # reset the Z and domZ values from the CLI
        self.Z = self._Z
        self.domZ = self._domZ

        #
        self._pli.do_alignment_score_calc = 0
        self._pli.long_targets = self.LONG_TARGETS

        with nogil:
            # reinitialize the random number generator
            if self._pli.do_reseeding:
                status = libeasel.random.esl_randomness_Init(self._pli.r, seed)
                if status != libeasel.eslOK:
                    raise UnexpectedError(status, "esl_randomness_Init")
            # reinitialize the domaindef
            libhmmer.p7_domaindef.p7_domaindef_Reuse(self._pli.ddef)

        # Configure reporting thresholds
        self._pli.by_E            = True
        self._pli.E               = 10.0
        self._pli.T               = 0.0
        self._pli.dom_by_E        = True
        self._pli.domE            = 10.0
        self._pli.domT            = 0.0
        self._pli.use_bit_cutoffs = False
        # if (go && esl_opt_IsOn(go, "-T"))
        #   {
        #     pli->T    = esl_opt_GetReal(go, "-T");
        #     pli->by_E = FALSE;
        #   }
        # if (go && esl_opt_IsOn(go, "--domT"))
        #   {
        #     pli->domT     = esl_opt_GetReal(go, "--domT");
        #     pli->dom_by_E = FALSE;
        #   }

        # Configure inclusion thresholds
        self._pli.inc_by_E           = True
        self._pli.incE               = 0.01
        self._pli.incT               = 0.0
        self._pli.incdom_by_E        = True
        self._pli.incdomE            = 0.01
        self._pli.incdomT            = 0.0
        # if (go && esl_opt_IsOn(go, "--incT"))
        #   {
        #     pli->incT     = esl_opt_GetReal(go, "--incT");
        #     pli->inc_by_E = FALSE;
        #   }
        # if (go && esl_opt_IsOn(go, "--incdomT"))
        #   {
        #     pli->incdomT     = esl_opt_GetReal(go, "--incdomT");
        #     pli->incdom_by_E = FALSE;
        #   }

        # Configure for one of the model-specific thresholding options
        # if (go && esl_opt_GetBoolean(go, "--cut_ga"))
        #   {
        #     pli->T        = pli->domT        = 0.0;
        #     pli->by_E     = pli->dom_by_E    = FALSE;
        #     pli->incT     = pli->incdomT     = 0.0;
        #     pli->inc_by_E = pli->incdom_by_E = FALSE;
        #     pli->use_bit_cutoffs = p7H_GA;
        #   }
        # if (go && esl_opt_GetBoolean(go, "--cut_nc"))
        #   {
        #     pli->T        = pli->domT        = 0.0;
        #     pli->by_E     = pli->dom_by_E    = FALSE;
        #     pli->incT     = pli->incdomT     = 0.0;
        #     pli->inc_by_E = pli->incdom_by_E = FALSE;
        #     pli->use_bit_cutoffs = p7H_NC;
        #   }
        # if (go && esl_opt_GetBoolean(go, "--cut_tc"))
        #   {
        #     pli->T        = pli->domT        = 0.0;
        #     pli->by_E     = pli->dom_by_E    = FALSE;
        #     pli->incT     = pli->incdomT     = 0.0;
        #     pli->inc_by_E = pli->incdom_by_E = FALSE;
        #     pli->use_bit_cutoffs = p7H_TC;
        #   }

        # Configure search space sizes for E value calculations
        self._pli.Z       = self._pli.domZ       = 0.0
        self._pli.Z_setby = self._pli.domZ_setby = p7_zsetby_e.p7_ZSETBY_NTARGETS
        # if (go && esl_opt_IsOn(go, "-Z"))
        #   {
        #     pli->Z_setby = p7_ZSETBY_OPTION;
        #     pli->Z       = esl_opt_GetReal(go, "-Z");
        #   }
        # if (go && esl_opt_IsOn(go, "--domZ"))
        #   {
        #     pli->domZ_setby = p7_ZSETBY_OPTION;
        #     pli->domZ       = esl_opt_GetReal(go, "--domZ");
        #   }

        # Configure acceleration pipeline thresholds
        self._pli.do_max        = False
        self._pli.do_biasfilter = self.bias_filter
        self._pli.do_null2      = self.null2
        self._pli.F1            = 0.02
        self._pli.F2            = 1e-3
        self._pli.F3            = 1e-5
        if self._pli.long_targets:
            self._pli.B1     = 100
            self._pli.B2     = 240
            self._pli.B3     = 1000
        else:
            self._pli.B1 = self._pli.B2 = self._pli.B3 = -1

        # if (go && esl_opt_GetBoolean(go, "--max"))
        #   {
        #     pli->do_max        = TRUE;
        #     pli->do_biasfilter = FALSE;
        #
        #     pli->F2 = pli->F3 = 1.0;
        #     pli->F1 = (pli->long_targets ? 0.3 : 1.0); // need to set some threshold for F1 even on long targets. Should this be tighter?
        #   }

        # Accounting as we collect results
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
        self._pli.show_accessions = True #(go && esl_opt_GetBoolean(go, "--acc")   ? TRUE  : FALSE);
        self._pli.show_alignments = False #(go && esl_opt_GetBoolean(go, "--noali") ? FALSE : TRUE);
        self._pli.hfp             = NULL
        self._pli.errbuf[0]       = b'\0'

    cpdef TopHits search_hmm(self, HMM query, object sequences):
        """search_hmm(self, query, sequences)\n--

        Run the pipeline using a query HMM against a sequence database.

        Arguments:
            query (`~pyhmmer.plan7.HMM`): The HMM object to use to query the
                sequence database.
            sequences (iterable of `~pyhmmer.easel.DigitalSequence`): The
                sequences to query with the HMM. For instance, pass a
                `~pyhmmer.easel.SequenceFile` in digital mode to read from
                disk iteratively.

        Returns:
            `~pyhmmer.plan7.TopHits`: the hits found in the sequence database.

        Raises:
            `ValueError`: When the alphabet of the current pipeline does not
                match the alphabet of the given HMM.

        """
        cdef int                  status
        cdef int                  allocM
        cdef DigitalSequence      seq
        cdef Profile              profile = self.profile
        cdef OptimizedProfile     opt
        cdef TopHits              hits    = TopHits()

        assert self._pli != NULL

        # check the pipeline was configure with the same alphabet
        if query.alphabet != self.alphabet:
            raise ValueError("Wrong alphabet in input HMM: expected {!r}, found {!r}".format(
                self.alphabet,
                query.alphabet
            ))

        # make sure the pipeline is set to search mode and ready for a new HMM
        self._pli.mode = p7_pipemodes_e.p7_SEARCH_SEQS

        # get an iterator over the input sequences, with an early return if
        # no sequences were given, and extract the raw pointer to the sequence
        # from the Python object
        seqs_iter = iter(sequences)
        seq = next(seqs_iter, None)
        if seq is None:
            return hits
        if seq.alphabet != self.alphabet:
            raise ValueError("Wrong alphabet in input sequence: expected {!r}, found {!r}".format(
                self.alphabet,
                seq.alphabet,
            ))

        # build the profile from the HMM, using the first sequence length
        # as a hint for the model configuration, or reuse the profile we
        # cached if it is large enough to hold the new HMM
        if profile is None or profile._gm.allocM < query.M:
            profile = self.profile = Profile(query.M, self.alphabet)
        else:
            profile.clear()
        profile.configure(query, self.background, len(seq))

        # build the optimized model from the profile
        opt = profile.optimized()

        # run the search loop on all database sequences while recycling memory
        self._search_loop(
            self._pli,
            opt._om,
            self.background._bg,
            seq._sq,
            hits._th,
            seqs_iter
        )

        # threshold hits
        hits.sort(by="key")
        hits.threshold(self)

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
            query (`~pyhmmer.plan7.DigitalSequence`): The sequence object to
                use to query the sequence database.
            sequences (iterable of `~pyhmmer.easel.DigitalSequence`): The
                sequences to query. Pass a `~pyhmmer.easel.SequenceFile`
                instance in digital mode to read from disk iteratively.
            builder (`~pyhmmer.plan7.Builder`, optional): A HMM builder to
                use to convert the query to a `~pyhmmer.plan7.HMM`. If `None`
                is given, it will use a default one.

        Returns:
            `~pyhmmer.plan7.TopHits`: the hits found in the sequence database.

        Raises:
            `ValueError`: When the alphabet of the current pipeline does not
                match the alphabet of the given query.

        """
        cdef int                  allocM
        cdef DigitalSequence      seq
        cdef Profile              profile
        cdef HMM                  hmm
        cdef OptimizedProfile     opt
        cdef TopHits              hits    = TopHits()

        assert self._pli != NULL

        # check the pipeline was configure with the same alphabet
        if query.alphabet != self.alphabet:
            raise ValueError("Wrong alphabet in query sequence: expected {!r}, found {!r}".format(
                self.alphabet,
                query.alphabet
            ))

        # make sure the pipeline is set to search mode and ready for a new HMM
        self._pli.mode = p7_pipemodes_e.p7_SEARCH_SEQS

        # use a default HMM builder if none was given
        builder = Builder(self.alphabet) if builder is None else builder

        # get an iterator over the input sequences, with an early return if
        # no sequences were given, and extract the raw pointer to the sequence
        # from the Python object
        seqs_iter = iter(sequences)
        seq = next(seqs_iter, None)
        if seq is None:
            return hits
        if seq.alphabet != self.alphabet:
            raise ValueError("Wrong alphabet in database sequence: expected {!r}, found {!r}".format(
                self.alphabet,
                seq.alphabet,
            ))

        # build the HMM and the profile from the query sequence, using the first
        # as a hint for the model configuration, or reuse the profile we
        # cached if it is large enough to hold the new HMM
        hmm, profile, opt = builder.build(query, self.background)

        # run the search loop on all database sequences while recycling memory
        self._search_loop(
            self._pli,
            opt._om,
            self.background._bg,
            seq._sq,
            hits._th,
            seqs_iter
        )

        # threshold hits
        hits.sort(by="key")
        hits.threshold(self)

        # return the hits
        return hits

    cdef  void    _search_loop(
        self,
        P7_PIPELINE* pli,
        P7_OPROFILE* om,
        P7_BG*       bg,
        ESL_SQ*      sq,
        P7_TOPHITS*  th,
        object       seqs_iter,
    ):
        cdef int      status
        cdef Sequence seq

        with nogil:
            # configure the pipeline for the current HMM
            status = libhmmer.p7_pipeline.p7_pli_NewModel(pli, om, bg)
            if status == libeasel.eslEINVAL:
                raise ValueError("model does not have bit score thresholds expected by the pipeline")
            elif status != libeasel.eslOK:
                raise UnexpectedError(status, "p7_pli_NewModel")

            # run the inner loop on all sequences
            while sq != NULL:
                # configure the profile, background and pipeline for the new sequence
                status = libhmmer.p7_pipeline.p7_pli_NewSeq(pli, sq)
                if status != libeasel.eslOK:
                    raise UnexpectedError(status, "p7_pli_NewSeq")
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

                # acquire the GIL just long enough to get the next sequence
                with gil:
                    seq = next(seqs_iter, None)
                    if seq is not None:
                        sq = seq._sq
                        if seq.alphabet != self.alphabet:
                            raise ValueError("Wrong alphabet in input sequence: expected {!r}, found {!r}".format(
                                self.alphabet,
                                seq.alphabet,
                            ))
                    else:
                        sq = NULL


cdef class Profile:
    """A Plan7 search profile.
    """

    def __cinit__(self):
        self._gm = NULL
        self.alphabet = None

    def __init__(self, int M, Alphabet alphabet):
        """__init__(self, M, alphabet)\n--

        Create a new profile for the given ``alphabet``.

        Arguments:
            M (`int`): The length of the profile, i.e. the number of nodes.
            alphabet (`~pyhmmer.easel.Alphabet`): The alphabet to use with
                this profile.

        """
        self.alphabet = alphabet
        with nogil:
            self._gm = libhmmer.p7_profile.p7_profile_Create(M, alphabet._abc)
        if not self._gm:
            raise AllocationError("P7_PROFILE")

    def __dealloc__(self):
        libhmmer.p7_profile.p7_profile_Destroy(self._gm)

    def __copy__(self):
        return self.copy()

    cpdef void clear(self):
        """clear(self)\n--

        Clear internal buffers to reuse the profile without reallocation.

        """
        assert self._gm != NULL
        cdef int status
        with nogil:
            status = libhmmer.p7_profile.p7_profile_Reuse(self._gm)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_profile_Reuse")

    cpdef void configure(self, HMM hmm, Background background, int L, bint multihit=True, bint local=True):
        """configure(self, hmm, background, L, multihit=True, local=True)\n--

        Configure a search profile using the given models.

        Arguments:
            hmm (`pyhmmer.plan7.HMM`): The model HMM with core probabilities.
            bg (`pyhmmer.plan7.Background`): The null background model.
            L (`int`): The expected target sequence length.
            multihit (`bool`): Whether or not to use multihit modes.
            local (`bool`): Whether or not to use non-local modes.

        """
        cdef int         mode
        cdef int         status
        cdef P7_HMM*     hm     = hmm._hmm
        cdef P7_BG*      bg     = background._bg
        cdef P7_PROFILE* gm     = self._gm

        if hmm.alphabet != self.alphabet:
            raise ValueError("HMM and profile alphabet don't match")
        elif hm.M > self._gm.allocM:
            raise ValueError("Profile is too small to hold HMM")

        if multihit:
            mode = p7_LOCAL if local else p7_GLOCAL
        else:
            mode = p7_UNILOCAL if local else p7_UNIGLOCAL

        with nogil:
            status = libhmmer.modelconfig.p7_ProfileConfig(hm, bg, gm, L, mode)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_ProfileConfig")

    cpdef Profile copy(self):
        """copy(self)\n--

        Return a copy of the profile with the exact same configuration.

        """
        assert self._gm != NULL

        cdef int status
        cdef Profile new = Profile.__new__(Profile)
        new.alphabet = self.alphabet
        with nogil:
            new._gm = libhmmer.p7_profile.p7_profile_Create(self._gm.allocM, self.alphabet._abc)
        if not new._gm:
            raise AllocationError("P7_PROFILE")

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
        cdef int              status
        cdef OptimizedProfile opt
        cdef P7_PROFILE*      gm     = self._gm

        assert self._gm != NULL

        opt = OptimizedProfile.__new__(OptimizedProfile)
        OptimizedProfile.__init__(opt, self._gm.M, self.alphabet)

        with nogil:
            status = p7_oprofile.p7_oprofile_Convert(gm, opt._om)
        if status == libeasel.eslOK:
            return opt
        elif status == libeasel.eslEINVAL:
            raise ValueError("Standard and optimized profiles are not compatible.")
        elif status == libeasel.eslEMEM:
            raise AllocationError("P7_OPROFILE")
        elif status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_oprofile_Convert")


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

    def __init__(self):
        """__init__(self)\n--

        Create an empty `TopHits` instance.

        """
        assert self._th == NULL, "called TopHits.__init__ more than once"
        with nogil:
            self._th = libhmmer.p7_tophits.p7_tophits_Create()
        if self._th == NULL:
            raise AllocationError("P7_TOPHITS")

    def __cinit__(self):
        self._th = NULL

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

    @property
    def reported(self):
        return self._th.nreported

    @property
    def included(self):
        return self._th.nincluded

    cdef void threshold(self, Pipeline pipeline):
        """threshold(self, pipeline)\n--

        Apply score and e-value thresholds using pipeline parameters.

        This function is automatically called in `Pipeline.search_hmm` or
        `Pipeline.search_seq`, and is therefore not exposed in the Python API.

        """
        cdef int          status
        cdef P7_TOPHITS*  th     = self._th
        cdef P7_PIPELINE* pli    = pipeline._pli

        assert self._th != NULL

        with nogil:
            status = libhmmer.p7_tophits.p7_tophits_Threshold(th, pli)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_tophits_Threshold")

        self.Z = pli.Z
        self.domZ = pli.domZ
        self.long_targets = pli.long_targets

    cpdef void clear(self):
        """clear(self)\n--

        Free internals to allow reusing for a new pipeline run.

        """
        assert self._th != NULL
        cdef int status = p7_tophits_Reuse(self._th)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_tophits_Reuse")

    cpdef void sort(self, str by="key"):
        """sort(self, by="key")\n--

        Sort hits in the current instance using the given method.

        Arguments:
            by (`str`): The comparison method to use to compare hits.
                Allowed values are: ``key`` (the default) to sort by key, or
                ``seqidx`` to sort by sequence index and alignment position.

        """
        assert self._th != NULL

        cdef P7_TOPHITS* th = self._th

        if by == "key":
            function = "p7_tophits_SortBySortkey"
            with nogil:
                status = libhmmer.p7_tophits.p7_tophits_SortBySortkey(th)
        elif by == "seqidx":
            function = "p7_tophits_SortBySeqidxAndAlipos"
            with nogil:
                status = libhmmer.p7_tophits.p7_tophits_SortBySeqidxAndAlipos(th)
        # elif by == "modelname"
        #     status = libhmmer.p7_tophits.p7_tophits_SortByModelnameAndAlipos(self._th)
        #     function = "p7_tophits_SortByModelnameAndAlipos"
        else:
            raise ValueError("Invalid value for `by` argument: {!r}".format(by))

        if status != libeasel.eslOK:
            raise UnexpectedError(status, function)

    cpdef bint is_sorted(self, str by="key"):
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

        raise ValueError("Invalid value for `by` argument: {!r}".format(by))


# --- Module init code -------------------------------------------------------

impl_Init()
p7_FLogsumInit()

include "exceptions.pxi"
