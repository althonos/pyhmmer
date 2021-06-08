# coding: utf-8
# cython: language_level=3, linetrace=True

# --- C imports --------------------------------------------------------------

from libc.stdint cimport int64_t, uint8_t, uint16_t, uint32_t
from posix.types cimport off_t

cimport libeasel.sq
from libeasel.alphabet cimport ESL_ALPHABET
from libeasel.bitfield cimport ESL_BITFIELD
from libeasel.keyhash cimport ESL_KEYHASH
from libeasel.msa cimport ESL_MSA
from libeasel.msafile cimport ESL_MSAFILE
from libeasel.random cimport ESL_RANDOMNESS
from libeasel.sq cimport ESL_SQ
from libeasel.sqio cimport ESL_SQFILE
from libeasel.ssi cimport ESL_SSI, ESL_NEWSSI


# --- Alphabet ---------------------------------------------------------------

cdef class Alphabet:
    cdef ESL_ALPHABET* _abc

    cdef void _init_default(self, int ty)
    cdef inline bint _eq(self, Alphabet other) nogil


# --- Bitfield ---------------------------------------------------------------

cdef class Bitfield:
    cdef ESL_BITFIELD* _b

    cdef size_t _wrap_index(self, int index) except -1
    cpdef size_t count(self, bint value=*)
    cpdef void toggle(self, int index)


# --- KeyHash ----------------------------------------------------------------

cdef class KeyHash:
    cdef ESL_KEYHASH* _kh

    cpdef int add(self, bytes key) except -1
    cpdef void clear(self)
    cpdef KeyHash copy(self)


# --- Matrix & Vector --------------------------------------------------------

cdef class Vector:
    cdef object _owner
    cdef int _n
    cdef readonly Py_ssize_t _shape[1]
    cdef readonly Py_ssize_t _strides[1]

cdef class VectorF(Vector):
    cdef float* _data

    cpdef int argmax(self)
    cpdef int argmin(self)
    cpdef VectorF copy(self)
    cpdef float max(self)
    cpdef float min(self)
    cpdef void normalize(self)
    cpdef void reverse(self)
    cpdef float sum(self)

cdef class VectorU8(Vector):
    cdef uint8_t* _data

    cpdef int argmax(self)
    cpdef int argmin(self)
    cpdef VectorU8 copy(self)
    cpdef uint8_t max(self)
    cpdef uint8_t min(self)
    cpdef void reverse(self)
    cpdef uint8_t sum(self)

cdef class Matrix:
    cdef object _owner
    cdef int _n
    cdef int _m
    cdef readonly Py_ssize_t _shape[2]
    cdef readonly Py_ssize_t _strides[2]

cdef class MatrixF(Matrix):
    cdef float** _data

    cpdef tuple argmax(self)
    cpdef tuple argmin(self)
    cpdef MatrixF copy(self)
    cpdef float max(self)
    cpdef float min(self)
    cpdef float sum(self)

cdef class MatrixU8(Matrix):
    cdef uint8_t** _data

    cpdef tuple argmax(self)
    cpdef tuple argmin(self)
    cpdef MatrixU8 copy(self)
    cpdef uint8_t max(self)
    cpdef uint8_t min(self)
    cpdef uint8_t sum(self)


# --- Multiple Sequences Alignment -------------------------------------------

cdef class _MSASequences:
    cdef MSA msa


cdef class MSA:
    cdef ESL_MSA* _msa

    cdef int _rehash(self) nogil except 1
    cpdef uint32_t checksum(self)
    cpdef void write(self, object fh, str format) except *


cdef class _TextMSASequences(_MSASequences):
    pass


cdef class TextMSA(MSA):
    cdef int _set_sequence(self, int idx, ESL_SQ* seq) nogil except 1
    cpdef TextMSA copy(self)
    cpdef DigitalMSA digitize(self, Alphabet alphabet)


cdef class _DigitalMSASequences(_MSASequences):
    cdef readonly Alphabet alphabet


cdef class DigitalMSA(MSA):
    cdef readonly Alphabet alphabet

    cdef int _set_sequence(self, int idx, ESL_SQ* seq) nogil except 1
    cpdef DigitalMSA copy(self)
    cpdef TextMSA textize(self)


# --- MSA File ---------------------------------------------------------------

cdef class MSAFile:
    cdef ESL_MSAFILE* _msaf
    cdef readonly Alphabet alphabet

    cpdef MSA read(self)

    cpdef void close(self)
    cpdef Alphabet guess_alphabet(self)
    cpdef Alphabet set_digital(self, Alphabet)


# --- Randomness -------------------------------------------------------------

cdef class Randomness:
    cdef ESL_RANDOMNESS* _rng
    cdef object          _owner

    cdef int _seed(self, uint32_t n) except 1

    cpdef void seed(self, object n=*)
    cpdef Randomness copy(self)
    cpdef double random(self)
    cpdef double normalvariate(self, double mu, double sigma)
    cpdef bint is_fast(self)


# --- Sequence ---------------------------------------------------------------

cdef class Sequence:
    cdef ESL_SQ* _sq

    cpdef void clear(self)
    cpdef uint32_t checksum(self)
    cpdef void write(self, object fh) except *


cdef class TextSequence(Sequence):
    cpdef TextSequence copy(self)
    cpdef DigitalSequence digitize(self, Alphabet alphabet)
    cpdef TextSequence reverse_complement(self, bint inplace=*)


cdef class DigitalSequence(Sequence):
    cdef readonly Alphabet alphabet

    cpdef DigitalSequence copy(self)
    cpdef TextSequence textize(self)
    cpdef DigitalSequence reverse_complement(self, bint inplace=*)


# --- Sequence File ----------------------------------------------------------

cdef class SequenceFile:
    cdef ESL_SQFILE* _sqfp
    cdef readonly Alphabet alphabet

    cpdef Sequence read(self, bint skip_info=*, bint skip_sequence=*)
    cpdef Sequence readinto(self, Sequence, bint skip_info=*, bint skip_sequence=*)

    # cpdef Sequence fetch(self, bytes key, bint skip_info=*, bint skip_sequence=*)
    # cpdef Sequence fetchinto(self, Sequence seq, bytes key, bint skip_info=*, bint skip_sequence=*)

    cpdef void close(self)
    cpdef Alphabet guess_alphabet(self)
    cpdef Alphabet set_digital(self, Alphabet)


# --- Sequence/Subsequence Index ---------------------------------------------

cdef class SSIReader:
    cdef ESL_SSI* _ssi

    cpdef void close(self)


cdef class SSIWriter:
    cdef ESL_NEWSSI* _newssi

    cdef void _on_write(self)

    cpdef void     add_alias(self, bytes alias, bytes key)
    cpdef uint16_t add_file(self, str filename, int format = *)
    cpdef void     add_key(
        self,
        bytes key,
        uint16_t fd,
        off_t record_offset,
        off_t data_offset = *,
        int64_t record_length = *
    )
