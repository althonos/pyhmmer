# coding: utf-8
# cython: language_level=3, linetrace=True

# --- C imports --------------------------------------------------------------

from libc.stdint cimport int64_t, uint16_t, uint32_t
from posix.types cimport off_t

cimport libeasel.sq
from libeasel.alphabet cimport ESL_ALPHABET
from libeasel.bitfield cimport ESL_BITFIELD
from libeasel.keyhash cimport ESL_KEYHASH
from libeasel.msa cimport ESL_MSA
from libeasel.msafile cimport ESL_MSAFILE
from libeasel.sq cimport ESL_SQ
from libeasel.sqio cimport ESL_SQFILE
from libeasel.ssi cimport ESL_SSI, ESL_NEWSSI


# --- Cython classes ---------------------------------------------------------

cdef class Alphabet:
    cdef ESL_ALPHABET* _abc

    cdef void _init_default(self, int ty)
    cdef inline bint _eq(self, Alphabet other) nogil


cdef class Bitfield:
    cdef ESL_BITFIELD* _b

    cdef size_t _wrap_index(self, int index) except -1
    cpdef size_t count(self, bint value=*)
    cpdef void toggle(self, int index)


cdef class KeyHash:
    cdef ESL_KEYHASH* _kh

    cpdef int add(self, bytes key) except -1
    cpdef void clear(self)
    cpdef KeyHash copy(self)


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


cdef class MSAFile:
    cdef ESL_MSAFILE* _msaf
    cdef readonly Alphabet alphabet

    cpdef MSA read(self)

    cpdef void close(self)
    cpdef Alphabet guess_alphabet(self)
    cpdef void set_digital(self, Alphabet)


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


cdef class SequenceFile:
    cdef ESL_SQFILE* _sqfp
    cdef readonly Alphabet alphabet

    cpdef Sequence read(self, bint skip_info=*, bint skip_sequence=*)
    cpdef Sequence readinto(self, Sequence, bint skip_info=*, bint skip_sequence=*)

    # cpdef Sequence fetch(self, bytes key, bint skip_info=*, bint skip_sequence=*)
    # cpdef Sequence fetchinto(self, Sequence seq, bytes key, bint skip_info=*, bint skip_sequence=*)

    cpdef void close(self)
    cpdef Alphabet guess_alphabet(self)
    cpdef void set_digital(self, Alphabet)


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
