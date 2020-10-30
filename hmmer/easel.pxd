from libeasel.alphabet cimport ESL_ALPHABET

cdef class Alphabet:
    cdef ESL_ALPHABET* _alphabet
    cdef _init_default(self, int ty)
