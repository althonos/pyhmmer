from libc.stdio cimport FILE
from posix.types cimport off_t

from libeasel.alphabet cimport ESL_ALPHABET
from libeasel.fileparser cimport ESL_FILEPARSER
from libeasel.ssi cimport ESL_SSI
from libhmmer.p7_hmm cimport P7_HMM


cdef extern from "easel.h" nogil:

    DEF eslERRBUFSIZE = 128


cdef extern from "hmmer.h" nogil:


    cdef enum p7_hmmfile_formats_e:
        p7_HMMFILE_20 = 0
        p7_HMMFILE_3a = 1
        p7_HMMFILE_3b = 2
        p7_HMMFILE_3c = 3
        p7_HMMFILE_3d = 4
        p7_HMMFILE_3e = 5
        p7_HMMFILE_3f = 6


    ctypedef p7_hmmfile_s P7_HMMFILE
    cdef struct p7_hmmfile_s:
        FILE* f
        char* fname
        ESL_SSI* ssi

        bint do_gzip
        bint do_stdin
        bint newly_opened
        bint is_pressed

        int format
        int (*parser)(p7_hmmfile_s*, ESL_ALPHABET**, P7_HMM**)
        ESL_FILEPARSER* efp

        FILE* ffp
        FILE* pfp

        # int syncRead
        # pthread_mutex_t readMutex

        char[eslERRBUFSIZE] errbuf


    int  p7_hmmfile_OpenE    (const char *filename, char *env, P7_HMMFILE **ret_hfp, char *errbuf)
    int  p7_hmmfile_OpenENoDB(const char *filename, char *env, P7_HMMFILE **ret_hfp, char *errbuf)
    int  p7_hmmfile_Open     (const char *filename, char *env, P7_HMMFILE **ret_hfp) # Deprecated
    int  p7_hmmfile_OpenNoDB (const char *filename, char *env, P7_HMMFILE **ret_hfp) # Deprecated
    int  p7_hmmfile_OpenBuffer(const char *buffer, int size, P7_HMMFILE **ret_hfp)
    void p7_hmmfile_Close(P7_HMMFILE *hfp)

    int  p7_hmmfile_WriteBinary(FILE *fp, int format, P7_HMM *hmm) except *
    int  p7_hmmfile_WriteASCII (FILE *fp, int format, P7_HMM *hmm) except *
    int  p7_hmmfile_WriteToString (char **s, int format, P7_HMM *hmm)

    int  p7_hmmfile_Read(P7_HMMFILE *hfp, ESL_ALPHABET **ret_abc,  P7_HMM **opt_hmm) except *
    int  p7_hmmfile_PositionByKey(P7_HMMFILE *hfp, const char *key)
    int  p7_hmmfile_Position(P7_HMMFILE *hfp, const off_t offset)
