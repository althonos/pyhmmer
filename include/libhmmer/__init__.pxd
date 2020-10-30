from libc.stdint cimport uint8_t


cdef extern from "hmmer.h" nogil:

    # Search modes
    DEF p7_NO_MODE = 0
    DEF p7_LOCAL = 1
    DEF p7_GLOCAL = 2
    DEF p7_UNILOCAL = 3
    DEF p7_UNIGLOCAL = 4
