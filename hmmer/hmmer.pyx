cimport libeasel
cimport libeasel.alphabet
from libeasel cimport eslERRBUFSIZE
from libeasel.alphabet cimport ESL_ALPHABET, esl_alphabet_Create, esl_abc_ValidateType
from libhmmer cimport p7_hmm, p7_hmmfile
from libhmmer.h2_io cimport p7_h2io_WriteASCII
from libhmmer.p7_hmm cimport P7_HMM
from libhmmer.p7_hmmfile cimport P7_HMMFILE

from .easel cimport Alphabet

import errno
import os
import io


cdef class Plan7HMM:

    # keep a reference to the Alphabet Python object to avoid deallocation of
    # the inner ESL_ALPHABET; the Python object provides reference counting
    # for free.
    cdef readonly Alphabet alphabet
    cdef P7_HMM* _hmm

    # --- Magic methods ------------------------------------------------------

    def __init__(self, int m, Alphabet alphabet):
        self.alphabet = alphabet
        self._hmm = p7_hmm.p7_hmm_Create(m, <const ESL_ALPHABET*> &alphabet._alphabet)
        if not self.hmm:
           raise MemoryError("could not allocate C object")

    def __cinit__(self):
        self.alphabet = None
        self._hmm = NULL

    def __dealloc__(self):
        p7_hmm.p7_hmm_Destroy(self._hmm)

    # --- Properties ---------------------------------------------------------

    @property
    def name(self):
        """`bytes`: The name of the HMM.
        """
        return <bytes> self._hmm.name

    @name.setter
    def name(self, bytes name):
        cdef int err = p7_hmm.p7_hmm_SetName(self._hmm, name)
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
        cdef int err = p7_hmm.p7_hmm_SetAccession(self._hmm, accession)
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
        cdef int err = p7_hmm.p7_hmm_SetDescription(self._hmm, description)
        if err == libeasel.eslEMEM:
            raise MemoryError("could not allocate C string")
        elif err != libeasel.eslOK:
            raise RuntimeError("unexpected error in p7_hmm_SetDescription: {}".format(err))

    # --- Methods ------------------------------------------------------------

    cpdef zero(self):
        """Set all parameters to zero (including model composition).
        """
        p7_hmm.p7_hmm_Zero(self._hmm)


cdef class Plam7HMMFile:

    cdef P7_HMMFILE* _hfp
    cdef Alphabet _alphabet

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self, str filename):
        cdef bytes fspath = os.fsencode(filename)
        cdef errbuf = bytearray(eslERRBUFSIZE)

        cdef err = p7_hmmfile.p7_hmmfile_OpenE(fspath, NULL, &self._hfp, errbuf)
        if err == libeasel.eslENOTFOUND:
            raise FileNotFoundError(errno.ENOENT, "no such file or directory: {!r}".format(filename))
        elif err == libeasel.eslEFORMAT:
            raise ValueError("format not recognized by HMMER")

        self._alphabet = Alphabet.__new__(Alphabet)
        self._alphabet._alphabet = NULL

    def __dealloc__(self):
        p7_hmmfile.p7_hmmfile_Close(self._hfp)

    def __iter__(self):
        # cdef int err = p7_hmmfile.p7_hmmfile_Position(self._hfp, 0)
        # if err == libeasel.eslEINVAL:
        #     raise io.UnsupportedOperation("file is not seekable")
        return self

    def __next__(self):
        cdef Plan7HMM py_hmm
        cdef P7_HMM* hmm = NULL

        cdef int err = p7_hmmfile.p7_hmmfile_Read(self._hfp, &self._alphabet._alphabet, &hmm)
        if err == libeasel.eslOK:
            py_hmm = Plan7HMM.__new__(Plan7HMM)
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
