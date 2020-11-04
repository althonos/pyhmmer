# coding: utf-8
# cython: language_level=3, linetrace=True

# --- C imports --------------------------------------------------------------

from libc.stdint cimport uint32_t

cimport libeasel.sq
from libeasel.alphabet cimport ESL_ALPHABET
from libeasel.sq cimport ESL_SQ
from libeasel.sqio cimport ESL_SQFILE


# --- Cython classes ---------------------------------------------------------

cdef class Alphabet:
    cdef ESL_ALPHABET* _abc
    cdef _init_default(self, int ty)


cdef class Sequence:
    cdef ESL_SQ* _sq

    cpdef void digitize(self, Alphabet alphabet)
    cpdef void clear(self)
    cpdef uint32_t checksum(self)


cdef class SequenceFile:
    # The wrapped ESL_SQFILE pointer
    cdef ESL_SQFILE* _sqfp
    # A reference to the Alphabet used if the file is in digital mode
    cdef Alphabet _alphabet

    cpdef Sequence read(self)
    cpdef Sequence read_info(self)
    cpdef Sequence read_seq(self)
    cpdef Sequence readinto(self, Sequence)
    cpdef Sequence readinto_info(self, Sequence)
    cpdef Sequence readinto_seq(self, Sequence)

    cpdef Sequence fetch(self, bytes)
    cpdef Sequence fetch_info(self, bytes)
    cpdef Sequence fetch_seq(self, bytes)
    cpdef Sequence fetchinto(self, Sequence, bytes)
    cpdef Sequence fetchinto_info(self, Sequence, bytes)
    cpdef Sequence fetchinto_seq(self, Sequence, bytes)

    cpdef void close(self)
    cpdef Alphabet guess_alphabet(self)
    cpdef void set_digital(self, Alphabet)
