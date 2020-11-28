from libc.stdint cimport uint32_t

from libeasel.alphabet cimport ESL_ALPHABET
from libhmmer.p7_hmm cimport P7_HMM
from libhmmer.p7_hmmfile cimport P7_HMMFILE

cdef extern from "reexports/p7_hmmfile.h" nogil:
    cdef uint32_t v3a_magic
    cdef uint32_t v3b_magic
    cdef uint32_t v3c_magic
    cdef uint32_t v3d_magic
    cdef uint32_t v3e_magic
    cdef uint32_t v3f_magic

    int read_asc20hmm(P7_HMMFILE* hfp, ESL_ALPHABET** ret_abc, P7_HMM** opt_hmm)
    int read_asc30hmm(P7_HMMFILE* hfp, ESL_ALPHABET** ret_abc, P7_HMM** opt_hmm)
    int read_bin30hmm(P7_HMMFILE* hfp, ESL_ALPHABET** ret_abc, P7_HMM** opt_hmm)
