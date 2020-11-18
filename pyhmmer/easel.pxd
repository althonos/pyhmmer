# coding: utf-8
# cython: language_level=3, linetrace=True

# --- C imports --------------------------------------------------------------

from libc.stdint cimport uint32_t

cimport libeasel.sq
from libeasel.alphabet cimport ESL_ALPHABET
from libeasel.bitfield cimport ESL_BITFIELD
from libeasel.keyhash cimport ESL_KEYHASH
from libeasel.sq cimport ESL_SQ
from libeasel.sqio cimport ESL_SQFILE


# --- Cython classes ---------------------------------------------------------

cdef class Alphabet:
    cdef ESL_ALPHABET* _abc

    cdef _init_default(self, int ty)


cdef class Bitfield:
    cdef ESL_BITFIELD* _b

    cdef size_t _wrap_index(self, int index)
    cpdef size_t count(self, bint value=*)
    cpdef void toggle(self, int index)


cdef class KeyHash:
    cdef ESL_KEYHASH* _kh

    cpdef void clear(self)
    cpdef KeyHash copy(self)


cdef class Sequence:
    cdef ESL_SQ* _sq

    cpdef void clear(self)
    cpdef uint32_t checksum(self)


cdef class TextSequence(Sequence):
    cpdef DigitalSequence digitize(self, Alphabet alphabet)


cdef class DigitalSequence(Sequence):
    cdef Alphabet alphabet


cdef class SequenceFile:
    cdef ESL_SQFILE* _sqfp
    cdef Alphabet _alphabet

    cpdef Sequence read(self, bint skip_info=*, bint skip_sequence=*)
    cpdef Sequence readinto(self, Sequence, bint skip_info=*, bint skip_sequence=*)

    cpdef Sequence fetch(self, bytes key, bint skip_info=*, bint skip_sequence=*)
    cpdef Sequence fetchinto(self, Sequence seq, bytes key, bint skip_info=*, bint skip_sequence=*)

    cpdef void close(self)
    cpdef Alphabet guess_alphabet(self)
    cpdef void set_digital(self, Alphabet)
