from libeasel.alphabet cimport ESL_ALPHABET
from libeasel.sq cimport ESL_SQ

cdef class Alphabet:
    cdef ESL_ALPHABET* _abc
    cdef _init_default(self, int ty)

cdef class Sequence:
    cdef ESL_SQ* _sq
