# coding: utf-8
# cython: language_level=3, linetrace=True

from libc.stdio cimport fprintf, FILE, stdout

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
from libeasel cimport eslERRBUFSIZE
from libeasel.alphabet cimport ESL_ALPHABET, esl_alphabet_Create, esl_abc_ValidateType
from libeasel.getopts cimport ESL_GETOPTS, ESL_OPTIONS
from libeasel.sq cimport ESL_SQ
from libhmmer cimport p7_LOCAL
from libhmmer.impl_sse.p7_oprofile cimport P7_OPROFILE
from libhmmer.p7_bg cimport P7_BG
from libhmmer.p7_hmm cimport P7_HMM
from libhmmer.p7_hmmfile cimport P7_HMMFILE
from libhmmer.p7_pipeline cimport P7_PIPELINE, p7_pipemodes_e
from libhmmer.p7_profile cimport P7_PROFILE
from libhmmer.p7_tophits cimport P7_TOPHITS

from .easel cimport Alphabet, Sequence

import errno
import os
import io


cdef extern from "hmmsearch.c":
    cdef ESL_OPTIONS* options


cdef class P7Profile:

    cdef P7_PROFILE* _gm

    def __cinit__(self):
        self._gm = NULL

    def __dealloc__(self):
        libhmmer.p7_profile.p7_profile_Destroy(self._gm)


cdef class P7HMM:

    # keep a reference to the Alphabet Python object to avoid deallocation of
    # the inner ESL_ALPHABET; the Python object provides reference counting
    # for free.
    cdef readonly Alphabet alphabet
    cdef P7_HMM* _hmm

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

    cpdef zero(self):
        """Set all parameters to zero (including model composition).
        """
        libhmmer.p7_hmm.p7_hmm_Zero(self._hmm)


cdef class P7HMMFile:

    cdef P7_HMMFILE* _hfp
    cdef Alphabet _alphabet

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self, str filename):
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
        libhmmer.p7_hmmfile.p7_hmmfile_Close(self._hfp)

    def __iter__(self):
        # cdef int err = p7_hmmfile.p7_hmmfile_Position(self._hfp, 0)
        # if err == libeasel.eslEINVAL:
        #     raise io.UnsupportedOperation("file is not seekable")
        return self

    def __next__(self):
        cdef P7HMM py_hmm
        cdef P7_HMM* hmm = NULL

        cdef int err = libhmmer.p7_hmmfile.p7_hmmfile_Read(self._hfp, &self._alphabet._abc, &hmm)
        if err == libeasel.eslOK:
            py_hmm = P7HMM.__new__(P7HMM)
            py_hmm.alphabet = self._alphabet # keep a reference to the alphabet
            py_hmm._hmm = hmm
            return py_hmm
        elif err == libeasel.eslEOF:
            raise StopIteration()
        elif err == libeasel.eslEMEM:
            raise MemoryError("could not allocate C object")
        elif err == libeasel.eslESYS:
            raise OSError(self._hfp.errbuf.decode('ascii'))
        elif err == libeasel.eslEFORMAT:
            raise ValueError("invalid format in file: {}".format(self._hfp.errbuf.decode('ascii')))
        elif err == libeasel.eslEINCOMPAT:
            alphabet = libeasel.alphabet.esl_abc_DecodeType(self._alphabet.type)
            raise ValueError("HMM is not in the expected {} alphabet".format(alphabet))
        else:
            raise RuntimeError("unexpected error ({}): {}".format(err, self._hfp.errbuf.decode('ascii')))


cpdef void hmmsearch(P7HMM hmm, Sequence seq):

    cdef int          status
    cdef FILE*        ofp     = stdout
    cdef P7_BG*       bg
    cdef P7_PROFILE*  gm
    cdef P7_OPROFILE* om
    cdef P7_TOPHITS*  th
    cdef P7_PIPELINE* pli
    cdef ESL_GETOPTS* go      = libeasel.getopts.esl_getopts_Create(options)
    cdef ESL_SQ*      dbsq    = libeasel.sq.esl_sq_Create()

    # built the top hits list
    th = libhmmer.p7_tophits.p7_tophits_Create();

    # copy and then digitize the sequence
    if libeasel.sq.esl_sq_Copy(seq._sq, dbsq):
        raise RuntimeError()
    if libeasel.sq.esl_sq_Digitize(hmm.alphabet._abc, dbsq):
        raise RuntimeError()

    # configure the profile and the background for the sequence
    # create the background model from the HMM alphabet
    bg = libhmmer.p7_bg.p7_bg_Create(hmm.alphabet._abc);
    if libhmmer.p7_bg.p7_bg_SetLength(bg, dbsq.n):
        raise RuntimeError()

    # build optimized model
    gm = libhmmer.p7_profile.p7_profile_Create(hmm._hmm.M, hmm.alphabet._abc)
    om = libhmmer.impl_sse.p7_oprofile.p7_oprofile_Create(hmm._hmm.M, hmm.alphabet._abc)
    if libhmmer.modelconfig.p7_ProfileConfig(hmm._hmm, bg, gm, 100, p7_LOCAL):
        raise RuntimeError()
    if libhmmer.impl_sse.p7_oprofile.p7_oprofile_Convert(gm, om):
        raise RuntimeError()
    if libhmmer.impl_sse.p7_oprofile.p7_oprofile_ReconfigLength(om, dbsq.n):
        raise RuntimeError()

    # build the processing pipeline
    pli = libhmmer.p7_pipeline.p7_pipeline_Create(go, om.M, 100, False, p7_pipemodes_e.p7_SEARCH_SEQS) # use same dummy as `hmmseach.c` code
    if libhmmer.p7_pipeline.p7_pli_NewSeq(pli, dbsq):
        raise RuntimeError()
    if libhmmer.p7_pipeline.p7_pli_NewModel(pli, om, bg):
        raise RuntimeError()

    # run the pipeline
    status = libhmmer.p7_pipeline.p7_Pipeline(pli, om, bg, dbsq, NULL, th);
    if status != libeasel.eslOK:
        raise RuntimeError()

    # sort and print results
    libhmmer.p7_tophits.p7_tophits_SortBySortkey(th)
    libhmmer.p7_tophits.p7_tophits_Threshold(th, pli)
    libhmmer.p7_tophits.p7_tophits_Domains(ofp, th, pli, True)
    fprintf(stdout, "\n\n");
