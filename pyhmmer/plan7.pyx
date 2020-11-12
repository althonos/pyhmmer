# coding: utf-8
# cython: language_level=3, linetrace=True
"""High-level interface to the Plan7 data model.

Plan7 is the model architecture used by HMMER since HMMER2.

See Also:
    Details about the Plan 7 architecture in the `HMMER documentation
    <http://www.csb.yale.edu/userguides/seq/hmmer/docs/node11.html>`_.

"""

# --- C imports --------------------------------------------------------------

from libc.stdint cimport uint32_t
from libc.stdio cimport fprintf, FILE, stdout
from libc.math cimport exp

cimport libeasel
cimport libeasel.sq
cimport libeasel.alphabet
cimport libeasel.random
cimport libeasel.getopts
cimport libhmmer.modelconfig
cimport libhmmer.p7_hmm
cimport libhmmer.p7_bg
cimport libhmmer.p7_domaindef
cimport libhmmer.p7_hmmfile
cimport libhmmer.p7_pipeline
cimport libhmmer.p7_profile
cimport libhmmer.p7_tophits
from libeasel cimport eslERRBUFSIZE, eslCONST_LOG2R
from libeasel.alphabet cimport ESL_ALPHABET, esl_alphabet_Create, esl_abc_ValidateType
from libeasel.getopts cimport ESL_GETOPTS, ESL_OPTIONS
from libeasel.sq cimport ESL_SQ
from libhmmer cimport p7_LOCAL
from libhmmer.p7_alidisplay cimport P7_ALIDISPLAY
from libhmmer.logsum cimport p7_FLogsumInit
from libhmmer.p7_pipeline cimport P7_PIPELINE, p7_pipemodes_e

IF HMMER_IMPL == "VMX":
    from libhmmer.impl_vmx cimport p7_oprofile
    from libhmmer.impl_vmx cimport impl_Init
    from libhmmer.impl_vmx.p7_oprofile cimport P7_OPROFILE
ELIF HMMER_IMPL == "SSE":
    from libhmmer.impl_sse cimport p7_oprofile
    from libhmmer.impl_sse cimport impl_Init
    from libhmmer.impl_sse.p7_oprofile cimport P7_OPROFILE

from .easel cimport Alphabet, Sequence


# --- Python imports ---------------------------------------------------------

import errno
import os
import io
import warnings
import collections.abc

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
        return self._ad.hmmfrom

    @property
    def hmm_to(self):
        return self._ad.hmmto

    @property
    def hmm_name(self):
        return <bytes> self._ad.hmmname

    @property
    def hmm_sequence(self):
        return self._ad.model.decode('ascii')

    @property
    def target_from(self):
        return self._ad.sqfrom

    @property
    def target_name(self):
        return <bytes> self._ad.sqname

    @property
    def target_sequence(self):
        return self._ad.aseq.decode('ascii')

    @property
    def target_to(self):
        return self._ad.sqto

    @property
    def identity_sequence(self):
        return self._ad.mline.decode('ascii')


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
        assert self._dom != NULL
        return exp(self._dom.lnP)

    # @property
    # def i_evalue(self):
    #     assert self._dom != NULL
    #     return exp(self._dom.lnP) * pli->domZ,


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
    """A data structure storing the Plan7 Hidden Markov Model.
    """

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
    """A wrapper around a file (or database), storing serialized HMMs.
    """

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

        if self._hfp == NULL:
            raise ValueError("I/O operation on closed file.")

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
        """Close the HMM file and free resources.

        This method has no effect if the file is already closed.
        """
        libhmmer.p7_hmmfile.p7_hmmfile_Close(self._hfp)
        self._hfp = NULL


cdef class Pipeline:
    """An HMMER3 accelerated sequence/profile comparison pipeline.
    """

    M_HINT = 100         # default model size
    L_HINT = 400         # default sequence size
    LONG_TARGETS = False

    def __cinit__(self):
        self._pli = NULL

    def __init__(
        self,
        Alphabet alphabet,
        *,
        bint bias_filter=True,
        float report_e=10.0,
        bint null2=True,
        object seed=None,
    ):
        """Instantiate and configure a new accelerated comparison pipeline.

        Arguments:
            alphabet (`~pyhmmer.easel.Alphabet`): The biological alphabet the
                of the HMMs and sequences that are going to be compared. Used
                to build the background model.

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
        # store a reference to the alphabet to avoid deallocation
        self.alphabet = alphabet

        # create the background model and the pipeline
        self._bg = libhmmer.p7_bg.p7_bg_Create(alphabet._abc)
        if self._bg == NULL:
            raise AllocationError("P7_BG")
        self._pli = libhmmer.p7_pipeline.p7_pipeline_Create(
            NULL,
            self.M_HINT,
            self.L_HINT,
            self.LONG_TARGETS,
            p7_pipemodes_e.p7_SEARCH_SEQS
        )
        if self._pli == NULL:
            raise AllocationError("P7_PIPELINE")

        # configure the pipeline with the additional keyword arguments
        self._pli.do_null2 = null2
        self._pli.do_biasfilter = bias_filter
        self._pli.E = report_e
        if seed is not None:  # we need to destroy and reallocate if given a seed
            libeasel.random.esl_randomness_Destroy(self._pli.r)
            self._pli.r = libeasel.random.esl_randomness_Create(<uint32_t> seed)
            if self._pli.r == NULL:
                raise AllocationError("ESL_RANDOMNESS")
            libhmmer.p7_domaindef.p7_domaindef_Destroy(self._pli.ddef)
            self._pli.ddef = libhmmer.p7_domaindef.p7_domaindef_Create(self._pli.r)
            if self._pli.ddef == NULL:
                raise AllocationError("P7_DOMAINDEF")
            self._pli.do_reseeding = self._pli.ddef.do_reseeding = seed == 0


    def __dealloc__(self):
        libhmmer.p7_pipeline.p7_pipeline_Destroy(self._pli)
        libhmmer.p7_bg.p7_bg_Destroy(self._bg)

    cpdef TopHits search(self, HMM hmm, object sequences, TopHits hits = None):
        """Run the pipeline using a query HMM against a sequence database.

        Arguments:
            hmm (`~pyhmmer.plan7.HMM`): The HMM object to use to query the
                sequence database.
            sequences (`Iterable` of `~pyhmmer.easel.Sequence`): The sequences
                to query with the HMMs. Pass a `~pyhmmer.easel.SequenceFile`
                instance to iteratively read from disk.
            hits (`~pyhmmer.plan7.TopHits`, optional): A hit collection to
                store results in, if any. If `None`, a new collection will
                be allocated. Use a non-`None` value if you are running
                several HMMs sequentially and want to aggregate the results,
                or if you are able to recycle a `~pyhmmer.plan7.TopHits`
                instance with `TopHits.clear`.

        Returns:
            `~pyhmmer.plan7.TopHits`: the hits found in the sequence database.
            If an instance was passed via the ``hits`` argument, it is safe
            to ignore the return value and to access it by reference.

        Raises:
            `ValueError`: When the alphabet of the current pipeline does not
                match the alphabet of the given HMM.

        """

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
        seqs_iter = iter(sequences)
        seq = next(seqs_iter, None)
        if seq is None:
            return hits
        sq = seq._sq

        # release the GIL for as long as possible to allow several search
        # pipelines to run in true parallel mode even in Python threads
        with nogil:

            # build the profile from the HMM, using the first sequence length
            # as a hint for the model configuraiton
            gm = libhmmer.p7_profile.p7_profile_Create(hm.M, abc)
            if libhmmer.modelconfig.p7_ProfileConfig(hm, bg, gm, sq.L, p7_LOCAL):
                raise RuntimeError()

            # build the optimized model from the HMM
            om = p7_oprofile.p7_oprofile_Create(hm.M, abc)
            if p7_oprofile.p7_oprofile_Convert(gm, om):
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
                if p7_oprofile.p7_oprofile_ReconfigLength(om, sq.n):
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
            p7_oprofile.p7_oprofile_Destroy(om)

        # return the hits
        return hits


cdef class Profile:
    """A Plan7 search profile.
    """

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
        """Clear internal buffers to reuse the profile without reallocation.
        """
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
    """A ranked list of top-scoring hits.
    """

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

    cpdef void threshold(self, Pipeline pipeline):
        """Apply score and e-value thresholds using ``pipeline`` parameters.

        This function will mark the targets and domains satisfying the
        reporting thresholds set for ``pipeline`` as significant, and update
        the total number of reported and included targets.

        """
        assert self._th != NULL
        libhmmer.p7_tophits.p7_tophits_Threshold(self._th, pipeline._pli)

    cpdef void clear(self):
        """Free internals to allow reusing for a new pipeline run.
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
