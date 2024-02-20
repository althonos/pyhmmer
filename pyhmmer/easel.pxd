# coding: utf-8
# cython: language_level=3, linetrace=True

# --- C imports --------------------------------------------------------------

from libc.stdint cimport int64_t, uint8_t, uint16_t, uint32_t
from posix.types cimport off_t

cimport libeasel.sq
from libeasel.alphabet cimport ESL_ALPHABET
from libeasel.bitfield cimport ESL_BITFIELD
from libeasel.gencode cimport ESL_GENCODE
from libeasel.keyhash cimport ESL_KEYHASH
from libeasel.msa cimport ESL_MSA
from libeasel.msafile cimport ESL_MSAFILE
from libeasel.random cimport ESL_RANDOMNESS
from libeasel.sq cimport ESL_SQ, ESL_SQ_BLOCK
from libeasel.sqio cimport ESL_SQFILE
from libeasel.ssi cimport ESL_SSI, ESL_NEWSSI
from libeasel cimport ESL_DSQ


# --- Alphabet ---------------------------------------------------------------

cdef class Alphabet:
    cdef ESL_ALPHABET* _abc

    cdef int _init_default(self, int ty) except 1 nogil
    cdef inline bint _eq(self, Alphabet other) nogil

    cpdef bint is_dna(self)
    cpdef bint is_rna(self)
    cpdef bint is_amino(self)
    cpdef bint is_nucleotide(self)
    cpdef VectorU8 encode(self, str sequence)
    cpdef str decode(self, const ESL_DSQ[::1] sequence)


# --- GeneticCode ------------------------------------------------------------

cdef class GeneticCode:
    cdef readonly Alphabet     amino_alphabet
    cdef readonly Alphabet     nucleotide_alphabet
    cdef          ESL_GENCODE* _gcode

    cdef int _translate(
        self,
        const ESL_DSQ* seq,
        int64_t seqlen,
        ESL_DSQ* out,
        int64_t outlen
    ) except -1 nogil

    cpdef VectorU8 translate(self, const ESL_DSQ[::1] sequence)


# --- Bitfield ---------------------------------------------------------------

cdef class Bitfield:
    cdef ESL_BITFIELD* _b
    cdef readonly Py_ssize_t _shape[1]

    cdef size_t _wrap_index(self, int index) except -1
    cpdef size_t count(self, bint value=*)
    cpdef Bitfield copy(self)
    cpdef void toggle(self, int index) except *


# --- KeyHash ----------------------------------------------------------------

cdef class KeyHash:
    cdef ESL_KEYHASH* _kh

    cpdef int add(self, bytes key) except -1
    cpdef void clear(self) except *
    cpdef KeyHash copy(self)


# --- Matrix & Vector --------------------------------------------------------

cdef class Vector:
    cdef object              _owner
    cdef int                 _n
    cdef readonly Py_ssize_t _shape[1]
    cdef void*               _data

    cdef const char* _format(self) noexcept
    cdef int _allocate(self, size_t n) except -1

cdef class VectorF(Vector):
    cpdef int argmax(self) except -1
    cpdef int argmin(self) except -1
    cpdef VectorF copy(self)
    cpdef float entropy(self) except *
    cpdef float max(self) except *
    cpdef float min(self) except *
    cpdef object normalize(self)
    cpdef float relative_entropy(self, VectorF other) except *
    cpdef object reverse(self)
    cpdef float sum(self)

cdef class VectorU8(Vector):
    cpdef int argmax(self) except -1
    cpdef int argmin(self) except -1
    cpdef VectorU8 copy(self)
    cpdef uint8_t max(self) except *
    cpdef uint8_t min(self) except *
    cpdef object reverse(self)
    cpdef uint8_t sum(self)

cdef class Matrix:
    cdef object              _owner
    cdef int                 _n
    cdef int                 _m
    cdef readonly Py_ssize_t _shape[2]
    cdef void**              _data

    cdef const char* _format(self) noexcept
    cdef int _allocate(self, size_t m, size_t n) except -1

cdef class MatrixF(Matrix):
    cpdef tuple argmax(self)
    cpdef tuple argmin(self)
    cpdef MatrixF copy(self)
    cpdef float max(self)
    cpdef float min(self)
    cpdef float sum(self)

cdef class MatrixU8(Matrix):
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

    cdef int _rehash(self) except 1 nogil
    cpdef uint32_t checksum(self)
    cpdef void write(self, object fh, str format) except *


cdef class _TextMSASequences(_MSASequences):
    pass


cdef class TextMSA(MSA):
    cdef int _set_sequence(self, int idx, ESL_SQ* seq) except 1 nogil
    cpdef TextMSA copy(self)
    cpdef DigitalMSA digitize(self, Alphabet alphabet)


cdef class _DigitalMSASequences(_MSASequences):
    cdef readonly Alphabet alphabet


cdef class DigitalMSA(MSA):
    cdef readonly Alphabet alphabet

    cdef int _set_sequence(self, int idx, ESL_SQ* seq) except 1 nogil
    cpdef DigitalMSA copy(self)
    cpdef TextMSA textize(self)


# --- MSA File ---------------------------------------------------------------

cdef class MSAFile:
    cdef          ESL_MSAFILE* _msaf
    cdef readonly object       _file
    cdef readonly Alphabet     alphabet
    cdef readonly str          name

    @staticmethod
    cdef ESL_MSAFILE* _open_fileobj(object fh, int fmt) except NULL

    cpdef void close(self)
    cpdef Alphabet guess_alphabet(self)
    cpdef MSA read(self)


# --- Randomness -------------------------------------------------------------

cdef class Randomness:
    cdef ESL_RANDOMNESS* _rng
    cdef object          _owner

    cpdef void seed(self, object n=*) except *
    cpdef Randomness copy(self)
    cpdef double random(self)
    cpdef double normalvariate(self, double mu, double sigma)


# --- Sequence ---------------------------------------------------------------

cdef class Sequence:
    cdef ESL_SQ* _sq

    cpdef void clear(self) except *
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
    cpdef DigitalSequence translate(self, GeneticCode genetic_code=*)
    cpdef DigitalSequence reverse_complement(self, bint inplace=*)


# --- Sequence Block ---------------------------------------------------------

cdef class SequenceBlock:
    cdef          size_t     _length   # the number of sequences in the block
    cdef          size_t     _capacity # the total number of sequences that can be stored
    cdef          ESL_SQ**   _refs     # the array to pass the sequence references
    cdef          list       _storage  # the actual Python list where `Sequence` objects are stored
    cdef          object     _owner    # the owner, if the object is just a shallow copy
    cdef          ssize_t    _largest  # the index of the largest sequence, or -1

    cdef void _on_modification(self) except *

    cdef void _allocate(self, size_t n) except *
    cdef void _append(self, Sequence sequence) except *
    cdef Sequence _pop(self, ssize_t index=*)
    cdef void _insert(self, ssize_t index, Sequence sequence) except *
    cdef size_t _index(self, Sequence sequence, ssize_t start=*, ssize_t stop=*) except *
    cdef void _remove(self, Sequence sequence) except *

    cpdef void extend(self, object iterable) except *
    cpdef SequenceBlock copy(self)
    cpdef void clear(self) except *
    cpdef Sequence largest(self)


cdef class TextSequenceBlock(SequenceBlock):
    cpdef void append(self, TextSequence sequence) except *
    cpdef TextSequence pop(self, ssize_t index=*)
    cpdef void insert(self, ssize_t index, TextSequence sequence) except *
    cpdef size_t index(self, TextSequence sequence, ssize_t start=*, ssize_t stop=*) except *
    cpdef void remove(self, TextSequence sequence) except *

    cpdef TextSequenceBlock copy(self)
    cpdef DigitalSequenceBlock digitize(self, Alphabet alphabet)
    cpdef TextSequence largest(self)

cdef class DigitalSequenceBlock(SequenceBlock):
    cdef readonly Alphabet   alphabet

    cpdef void append(self, DigitalSequence sequence) except *
    cpdef DigitalSequence pop(self, ssize_t index=*)
    cpdef void insert(self, ssize_t index, DigitalSequence sequence) except *
    cpdef size_t index(self, DigitalSequence sequence, ssize_t start=*, ssize_t stop=*) except *
    cpdef void remove(self, DigitalSequence sequence) except *

    cpdef DigitalSequenceBlock copy(self)
    cpdef TextSequenceBlock textize(self)
    cpdef DigitalSequenceBlock translate(self, GeneticCode genetic_code=*)
    cpdef DigitalSequence largest(self)

# --- Sequence File ----------------------------------------------------------

cdef class SequenceFile:
    cdef          ESL_SQFILE* _sqfp
    cdef readonly object      _file
    cdef readonly Alphabet    alphabet
    cdef readonly str         name

    @staticmethod
    cdef ESL_SQFILE* _open_fileobj(object fh, int fmt) except NULL

    cpdef void close(self) except *
    cpdef Alphabet guess_alphabet(self)
    cpdef Sequence read(self, bint skip_info=*, bint skip_sequence=*)
    cpdef Sequence readinto(self, Sequence, bint skip_info=*, bint skip_sequence=*)
    cpdef SequenceBlock read_block(self, object sequences=*, object residues=*)
    cpdef void rewind(self) except *

# --- Sequence/Subsequence Index ---------------------------------------------

cdef class SSIReader:
    cdef ESL_SSI* _ssi

    cpdef void close(self)


cdef class SSIWriter:
    cdef ESL_NEWSSI* _newssi

    cdef void _on_write(self)

    cpdef void     add_alias(self, bytes alias, bytes key) except *
    cpdef uint16_t add_file(self, object filename, int format = *) except *
    cpdef void     add_key(
        self,
        bytes key,
        uint16_t fd,
        off_t record_offset,
        off_t data_offset = *,
        int64_t record_length = *
    ) except *
    cpdef void close(self) except *
