from libc.stdio cimport FILE

from libhmmer.p7_hmm cimport P7_HMM


cdef extern from "hmmer.h":

    int p7_h2io_WriteASCII(FILE* fp, P7_HMM* hmm)
