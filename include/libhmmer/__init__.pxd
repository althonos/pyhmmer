from libc.stdint cimport uint8_t


cdef extern from "hmmer.h" nogil:

    # Search modes
    cdef enum:
        p7_NO_MODE = 0
        p7_LOCAL = 1
        p7_GLOCAL = 2
        p7_UNILOCAL = 3
        p7_UNIGLOCAL = 4
