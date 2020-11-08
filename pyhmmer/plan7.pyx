# coding: utf-8
# cython: language_level=3, linetrace=True

# --- C imports --------------------------------------------------------------

from libc.stdio cimport fprintf, FILE, stdout
from libc.math cimport exp

cimport libeasel
cimport libeasel.sq
cimport libeasel.alphabet
cimport libeasel.getopts
cimport libhmmer.impl_sse.p7_oprofile
cimport libhmmer.modelconfig
cimport libhmmer.p7_hmm
cimport libhmmer.p7_bg
cimport libhmmer.p7_hmmfile
cimport libhmmer.p7_pipeline
cimport libhmmer.p7_profile
cimport libhmmer.p7_tophits
from libeasel cimport eslERRBUFSIZE, eslCONST_LOG2R
from libeasel.alphabet cimport ESL_ALPHABET, esl_alphabet_Create, esl_abc_ValidateType
from libeasel.getopts cimport ESL_GETOPTS, ESL_OPTIONS
from libeasel.sq cimport ESL_SQ
from libhmmer cimport p7_LOCAL
from libhmmer.logsum cimport p7_FLogsumInit
from libhmmer.impl_sse cimport impl_Init
from libhmmer.impl_sse.p7_oprofile cimport P7_OPROFILE
from libhmmer.p7_pipeline cimport P7_PIPELINE, p7_pipemodes_e

from .easel cimport Alphabet, Sequence


# --- Python imports ---------------------------------------------------------

import errno
import os
import io
import warnings

from .errors import AllocationError, UnexpectedError


# --- Cython classes ---------------------------------------------------------

cdef class Domain:

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self.hit = None
        self._dom = NULL

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
    def ali_from(self):
        assert self._dom != NULL
        return self._dom.iali

    @property
    def ali_to(self):
        assert self._dom != NULL
        return self._dom.jali

    @property
    def hmm_from(self):
        assert self._dom != NULL
        return self._dom.iorf

    @property
    def hmm_to(self):
        assert self._dom != NULL
        return self._dom.jorf

    @property
    def score(self):
        """`float`: The overall score in BITS, null-corrected.
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
        assert self._dom != NULL
        return exp(self._dom.lnP)

    # @property
    # def c_evalue(self):
    #     assert self._dom != NULL
    #     return exp(self._dom.lnP) * pli->domZ,
    #
    # @property
    # def i_evalue(self):
    #     assert self._dom != NULL
    #     exp(th->hit[h]->dcl[d].lnP) * pli->Z,


cdef class Domains:

    def __cinit__(self, Hit hit):
        self.hit = hit

    def __len__(self):
        return self.hit._hit.ndom

    def __getitem__(self, index):
        if index < 0:
            index += self.hit._hit.ndom
        if index >= self.hit._hit.ndom or index < 0:
            raise IndexError("list index out of range")

        cdef Domain dom = Domain.__new__(Domain)
        dom.hit = self.hit
        dom._dom = &self.hit._hit.dcl[<size_t> index]

        return dom


cdef class Hit:

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self, TopHits hits, size_t index):
        assert hits._th != NULL
        assert index < hits._th.N

        self.hits = hits
        self._hit = hits._th.hit[index]

    # --- Properties ---------------------------------------------------------

    @property
    def name(self):
        assert self._hit != NULL
        assert self._hit.name != NULL
        return <bytes> self._hit.name

    @property
    def accession(self):
        assert self._hit != NULL
        if self._hit.acc == NULL:
            return None
        return <bytes> self._hit.acc

    @property
    def description(self):
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
        """`float`: Bit score of the sequence before null2 correction.
        """
        assert self._hit != NULL
        return self._hit.pre_score

    @property
    def bias(self):
        assert self._hit != NULL
        return self._hit.pre_score - self._hit.score

    @property
    def lnP(self):
        assert self._hit != NULL
        return self._hit.lnP

    @property
    def domains(self):
        return Domains(self)


cdef class HMM:

    # --- Magic methods ------------------------------------------------------

    def __init__(self, int m, Alphabet alphabet):
        self.alphabet = alphabet
        self._hmm = libhmmer.p7_hmm.p7_hmm_Create(m, alphabet._abc)
        if not self.hmm:
           raise MemoryError("could not allocate C object")

    def __cinit__(self):
        self.alphabet = None
        self._hmm = NULL

    def __dealloc__(self):
        libhmmer.p7_hmm.p7_hmm_Destroy(self._hmm)

    # --- Properties ---------------------------------------------------------

    @property
    def name(self):
        """`bytes`: The name of the HMM.
        """
        return <bytes> self._hmm.name

    @name.setter
    def name(self, bytes name):
        cdef int err = libhmmer.p7_hmm.p7_hmm_SetName(self._hmm, name)
        if err == libeasel.eslEMEM:
            raise MemoryError("could not allocate C string")
        elif err != libeasel.eslOK:
            raise RuntimeError("unexpected error in p7_hmm_SetName: {}".format(err))

    @property
    def accession(self):
        """`bytes`: The accession of the HMM.
        """
        return <bytes> self._hmm.acc

    @accession.setter
    def accession(self, bytes accession):
        cdef int err = libhmmer.p7_hmm.p7_hmm_SetAccession(self._hmm, accession)
        if err == libeasel.eslEMEM:
            raise MemoryError("could not allocate C string")
        elif err != libeasel.eslOK:
            raise RuntimeError("unexpected error in p7_hmm_SetAccession: {}".format(err))

    @property
    def description(self):
        """`bytes`: The description of the HMM.
        """
        return <bytes> self._hmm.desc

    @description.setter
    def description(self, bytes description):
        cdef int err = libhmmer.p7_hmm.p7_hmm_SetDescription(self._hmm, description)
        if err == libeasel.eslEMEM:
            raise MemoryError("could not allocate C string")
        elif err != libeasel.eslOK:
            raise RuntimeError("unexpected error in p7_hmm_SetDescription: {}".format(err))

    # --- Methods ------------------------------------------------------------

    cpdef void zero(self):
        """Set all parameters to zero (including model composition).
        """
        libhmmer.p7_hmm.p7_hmm_Zero(self._hmm)


cdef class HMMFile:

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._alphabet = None
        self._hfp = NULL

    def __init__(self, str filename):
        cdef bytes fspath = os.fsencode(filename)
        cdef errbuf = bytearray(eslERRBUFSIZE)

        cdef err = libhmmer.p7_hmmfile.p7_hmmfile_OpenE(fspath, NULL, &self._hfp, errbuf)
        if err == libeasel.eslENOTFOUND:
            raise FileNotFoundError(errno.ENOENT, "no such file or directory: {!r}".format(filename))
        elif err == libeasel.eslEFORMAT:
            raise ValueError("format not recognized by HMMER")

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

        with nogil:
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


    # --- Utils --------------------------------------------------------------

    cpdef void close(self):
        libhmmer.p7_hmmfile.p7_hmmfile_Close(self._hfp)
        self._hfp = NULL


cdef class Pipeline:

    def __cinit__(self):
        self._pli = NULL

    def __init__(self, Alphabet alphabet):

        M_hint = 100
        L_hint = 400
        long_targets = False

        self.alphabet = alphabet

        self._pli = libhmmer.p7_pipeline.p7_pipeline_Create(NULL, M_hint, L_hint, long_targets, p7_pipemodes_e.p7_SEARCH_SEQS)
        if self._pli == NULL:
            raise AllocationError("P7_PIPELINE")

        self._bg = libhmmer.p7_bg.p7_bg_Create(alphabet._abc)
        if self._bg == NULL:
            raise AllocationError("P7_BG")

    def __dealloc__(self):
        libhmmer.p7_pipeline.p7_pipeline_Destroy(self._pli)
        libhmmer.p7_bg.p7_bg_Destroy(self._bg)

    cpdef TopHits search(self, HMM hmm, object seqs, TopHits hits = None):
        cdef int           status
        cdef Sequence      seq
        cdef ESL_ALPHABET* abc     = self.alphabet._abc
        cdef P7_BG*        bg      = self._bg
        cdef P7_HMM*       hm      = hmm._hmm
        cdef P7_PROFILE*   gm
        cdef P7_OPROFILE*  om
        cdef P7_TOPHITS*   th
        cdef P7_PIPELINE*  pli     = self._pli
        cdef ESL_SQ*       sq      = NULL

        assert self._pli != NULL

        # check the pipeline was configure with the same alphabet
        if hmm.alphabet != self.alphabet:
            raise ValueError("Wrong alphabet in input HMM: expected {!r}, found {!r}".format(
              self.alphabet,
              hmm.alphabet
            ))

        # make sure the pipeline is set to search mode
        self._pli.mode = p7_pipemodes_e.p7_SEARCH_SEQS

        # get a pointer to the P7_TOPHITS struct to use
        hits = TopHits() if hits is None else hits
        th = hits._th

        # get an iterator over the input sequences, with an early return if
        # no sequences were given, and extract the raw pointer to the sequence
        # from the Python object
        seqs_iter = iter(seqs)
        seq = next(seqs_iter, None)
        if seq is None:
            return hits
        sq = seq._sq

        # release the GIL for as long as possible
        with nogil:

            # build the profile from the HMM, using the first sequence length
            # as a hint for the model configuraiton
            gm = libhmmer.p7_profile.p7_profile_Create(hm.M, abc)
            if libhmmer.modelconfig.p7_ProfileConfig(hm, bg, gm, sq.L, p7_LOCAL):
                raise RuntimeError()

            # build the optimized model from the HMM
            om = libhmmer.impl_sse.p7_oprofile.p7_oprofile_Create(hm.M, abc)
            if libhmmer.impl_sse.p7_oprofile.p7_oprofile_Convert(gm, om):
                raise RuntimeError()

            # configure the pipeline for the current HMM
            if libhmmer.p7_pipeline.p7_pli_NewModel(pli, om, bg):
                raise RuntimeError()

            # run the inner loop on all sequences
            while sq != NULL:

                # digitize the sequence if needed
                # if libeasel.sq.esl_sq_Copy(seq._sq, dbsq):
                #     raise RuntimeError()
                if not libeasel.sq.esl_sq_IsDigital(sq):
                    if libeasel.sq.esl_sq_Digitize(abc, sq):
                        raise RuntimeError()

                # configure the profile, background and pipeline for the new sequence
                if libhmmer.p7_pipeline.p7_pli_NewSeq(pli, sq):
                    raise RuntimeError()
                if libhmmer.p7_bg.p7_bg_SetLength(bg, sq.n):
                    raise RuntimeError()
                if libhmmer.impl_sse.p7_oprofile.p7_oprofile_ReconfigLength(om, sq.n):
                    raise RuntimeError()

                # run the pipeline on the sequence
                status = libhmmer.p7_pipeline.p7_Pipeline(pli, om, bg, sq, NULL, th);
                if status != libeasel.eslOK:
                    raise RuntimeError()

                # clear pipeline for reuse for next target
                libhmmer.p7_pipeline.p7_pipeline_Reuse(pli)

                # acquire the GIL just long enough to get the next sequence
                with gil:
                    seq = next(seqs_iter, None)
                    sq = NULL if seq is None else seq._sq

            # deallocate the profile and optimized model
            libhmmer.p7_profile.p7_profile_Destroy(gm)
            libhmmer.impl_sse.p7_oprofile.p7_oprofile_Destroy(om)

        #
        return hits


cdef class Profile:

    def __cinit__(self):
        self._gm = NULL

    def __init__(self, int m, Alphabet alphabet):
        self.alphabet = alphabet
        self._gm = libhmmer.p7_profile.p7_profile_Create(m, alphabet._abc)
        if not self._gm:
            raise AllocationError("P7_PROFILE")

    def __dealloc__(self):
        libhmmer.p7_profile.p7_profile_Destroy(self._gm)

    def __copy__(self):
        return self.copy()

    cpdef void clear(self):
        assert self._gm != NULL
        cdef int status = libhmmer.p7_profile.p7_profile_Reuse(self._gm)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_profile_Reuse")

    cpdef Profile copy(self):
        assert self._gm != NULL

        cdef Profile new = Profile.__new__(Profile)
        new.alphabet = self.alphabet
        new._gm = libhmmer.p7_profile.p7_profile_Create(self._gm.allocM, self.alphabet._abc)
        if not new._gm:
            raise AllocationError("P7_PROFILE")

        cdef int status = libhmmer.p7_profile.p7_profile_Copy(self._gm, new._gm)
        if status == libeasel.eslOK:
            return new
        else:
            raise UnexpectedError(status, "p7_profile_Copy")

    cpdef bint is_local(self):
        """Returns whether or not the profile is in a local alignment mode.
        """
        assert self._gm != NULL
        return libhmmer.p7_profile.p7_profile_IsLocal(self._gm)

    cpdef bint is_multihit(self):
        """Returns whether or not the profile is in a multihit alignment mode.
        """
        assert self._gm != NULL
        return libhmmer.p7_profile.p7_profile_IsMultihit(self._gm)


cdef class TopHits:

    def __init__(self):
        assert self._th == NULL, "called TopHits.__init__ more than once"
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
        cdef size_t index_
        if index < 0:
            index += self._th.N
        if index >= self._th.N or index < 0:
            raise IndexError("list index out of range")
        return Hit(self, index)

    def __iadd__(self, other):
        assert self._th != NULL

        if not isinstance(other, TopHits):
            self_ty = type(self).__name__
            other_ty = type(other).__name__
            return TypeError("Cannot merge {!r} object into a {!r} instance".format(other_ty, self_ty))

        cdef int status = libhmmer.p7_tophits.p7_tophits_Merge(self._th, (<TopHits> other)._th)
        if status == libeasel.eslOK:
            libhmmer.p7_tophits.p7_tophits_Reuse((<TopHits> other)._th)
            return self
        elif status == libeasel.eslEMEM:
            raise AllocationError("P7_TOPHITS")
        else:
            raise UnexpectedError(status, "p7_tophits_Merge")

    cpdef void sort_by_sortkey(self):
        assert self._th != NULL
        libhmmer.p7_tophits.p7_tophits_SortBySortkey(self._th)

    cpdef void threshold(self, Pipeline pipeline):
        assert self._th != NULL
        libhmmer.p7_tophits.p7_tophits_Threshold(self._th, pipeline._pli)

    cpdef void clear(self):
        """Free internals to allow reusing the top hits for something else.
        """
        assert self._th != NULL
        cdef int status = libhmmer.p7_tophits.p7_tophits_Reuse(self._th)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_tophits_Reuse")

    cpdef void sort(self, str by="key"):
        assert self._th != NULL

        if by == "key":
            status = libhmmer.p7_tophits.p7_tophits_SortBySortkey(self._th)
            function = "p7_tophits_SortBySortkey"
        elif by == "seqidx":
            status = libhmmer.p7_tophits.p7_tophits_SortBySeqidxAndAlipos(self._th)
            function = "p7_tophits_SortBySeqidxAndAlipos"
        # elif by == "modelname"
        #     status = libhmmer.p7_tophits.p7_tophits_SortByModelnameAndAlipos(self._th)
        #     function = "p7_tophits_SortByModelnameAndAlipos"
        else:
            raise ValueError("Invalid value for `by` argument: {!r}".format(by))

        if status != libeasel.eslOK:
            raise UnexpectedError(status, function)

    cpdef bint is_sorted(self, str by="key"):
        assert self._th != NULL

        if by == "key":
            return self._th.is_sorted_by_sortkey
        elif by == "seqidx":
            return self._th.is_sorted_by_seqidx

        raise ValueError("Invalid value for `by` argument: {!r}".format(by))




    # cpdef void print_targets(self, Pipeline pipeline):
    #     assert self._th != NULL
    #     libhmmer.p7_tophits.p7_tophits_Targets(stdout, self._th, pipeline._pli, 0)
    #
    # cpdef void print_domains(self, Pipeline pipeline):
    #     assert self._th != NULL
    #     libhmmer.p7_tophits.p7_tophits_Domains(stdout, self._th, pipeline._pli, 0)


# --- Module init code -------------------------------------------------------

impl_Init()
p7_FLogsumInit()
