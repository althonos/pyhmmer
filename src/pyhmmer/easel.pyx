# coding: utf-8
# cython: language_level=3
"""High-level interface to the Easel C library.

Easel is a library developed by the `Eddy/Rivas Lab <http://eddylab.org/>`_
to facilitate the development of biological software in C. It is used by
`HMMER <http://hmmer.org/>`_ and `Infernal <http://eddylab.org/infernal/>`_.

"""

# --- C imports --------------------------------------------------------------

cimport cython
from cpython cimport Py_buffer
from cpython.buffer cimport PyBUF_FORMAT, PyBUF_READ, PyBuffer_FillInfo
from cpython.bytes cimport PyBytes_FromString, PyBytes_FromStringAndSize, PyBytes_AsString
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from cpython.memoryview cimport PyMemoryView_FromMemory
from cpython.ref cimport Py_INCREF
from cpython.tuple cimport PyTuple_New, PyTuple_SET_ITEM
from libc.stdint cimport int32_t, int64_t, uint8_t, uint16_t, uint32_t, uint64_t, SIZE_MAX
from libc.stdio cimport fclose, FILE
from libc.stdlib cimport calloc, malloc, realloc, free
from libc.string cimport memcmp, memcpy, memmove, memset, strdup, strlen, strncpy
from posix.types cimport off_t
from cpython.unicode cimport (
    PyUnicode_New,
    PyUnicode_DecodeASCII,
    PyUnicode_DATA,
    PyUnicode_KIND,
    PyUnicode_WRITE,
    PyUnicode_READ,
)

cimport libeasel
cimport libeasel.alphabet
cimport libeasel.bitfield
cimport libeasel.buffer
cimport libeasel.gencode
cimport libeasel.keyhash
cimport libeasel.matrixops
cimport libeasel.msa
cimport libeasel.msafile
cimport libeasel.msaweight
cimport libeasel.random
cimport libeasel.sq
cimport libeasel.sqio
cimport libeasel.sqio.ascii
cimport libeasel.ssi
cimport libeasel.vec
from libeasel cimport ESL_DSQ, esl_pos_t, eslERRBUFSIZE
from libeasel.buffer cimport ESL_BUFFER
from libeasel.gencode cimport ESL_GENCODE
from libeasel.msaweight cimport ESL_MSAWEIGHT_CFG
from libeasel.sq cimport ESL_SQ
from libeasel.sqio cimport ESL_SQFILE, ESL_SQASCII_DATA
from libeasel.random cimport ESL_RANDOMNESS
from capacity cimport new_capacity

from .reexports.esl_sqio_ascii cimport (
    loadbuf,
    sqascii_Position,
    sqascii_Close,
    sqascii_SetDigital,
    sqascii_GuessAlphabet,
    sqascii_IsRewindable,
    sqascii_Read,
    sqascii_ReadInfo,
    sqascii_ReadSequence,
    sqascii_ReadWindow,
    sqascii_Echo,
    sqascii_ReadBlock,
    sqascii_OpenSSI,
    sqascii_PositionByKey,
    sqascii_PositionByNumber,
    sqascii_Fetch,
    sqascii_FetchInfo,
    sqascii_FetchSubseq,
    sqascii_GetError,
    sqascii_GuessFileFormat,
    config_embl,
    config_genbank,
    config_fasta,
    config_daemon,
    inmap_embl,
    inmap_genbank,
    inmap_fasta,
    inmap_daemon,
    fileheader_hmmpgmd,
)

if TARGET_SYSTEM == "Linux":
    from .fileobj.linux cimport fileobj_linux_open as fopen_obj
elif TARGET_SYSTEM == "Darwin" or TARGET_SYSTEM.endswith("BSD"):
    from .fileobj.bsd cimport fileobj_bsd_open as fopen_obj

include "exceptions.pxi"
include "_getid.pxi"

# --- Python imports ---------------------------------------------------------

import abc
import array
import errno
import functools
import io
import itertools
import os
import operator
import collections
import pickle
import sys
import warnings

from .errors import AllocationError, UnexpectedError, AlphabetMismatch, InvalidParameter
from .utils import peekable

# --- Constants --------------------------------------------------------------

cdef dict MSA_FILE_FORMATS = {
    "stockholm": libeasel.msafile.eslMSAFILE_STOCKHOLM,
    "pfam": libeasel.msafile.eslMSAFILE_PFAM,
    "a2m": libeasel.msafile.eslMSAFILE_A2M,
    "psiblast": libeasel.msafile.eslMSAFILE_PSIBLAST,
    "selex": libeasel.msafile.eslMSAFILE_SELEX,
    "afa": libeasel.msafile.eslMSAFILE_AFA,
    "clustal": libeasel.msafile.eslMSAFILE_CLUSTAL,
    "clustallike": libeasel.msafile.eslMSAFILE_CLUSTALLIKE,
    "phylip": libeasel.msafile.eslMSAFILE_PHYLIP,
    "phylips": libeasel.msafile.eslMSAFILE_PHYLIPS,
}

cdef dict MSA_FILE_FORMATS_INDEX = {
    v:k for k,v in MSA_FILE_FORMATS.items()
}

cdef dict SEQUENCE_FILE_FORMATS = {
    "fasta": libeasel.sqio.eslSQFILE_FASTA,
    "embl": libeasel.sqio.eslSQFILE_EMBL,
    "genbank": libeasel.sqio.eslSQFILE_GENBANK,
    "ddbj": libeasel.sqio.eslSQFILE_DDBJ,
    "uniprot": libeasel.sqio.eslSQFILE_UNIPROT,
    "ncbi": libeasel.sqio.eslSQFILE_NCBI,
    "daemon": libeasel.sqio.eslSQFILE_DAEMON,
    "hmmpgmd": libeasel.sqio.eslSQFILE_DAEMON,
    "fmindex": libeasel.sqio.eslSQFILE_FMINDEX,
    **MSA_FILE_FORMATS,
}

cdef dict SEQUENCE_FILE_FORMATS_INDEX = {
    v:k for k,v in SEQUENCE_FILE_FORMATS.items()
}

cdef dict MSA_WEIGHT_PREFERENCES = {
    "conscover": libeasel.msaweight.eslMSAWEIGHT_FILT_CONSCOVER,
    "origorder": libeasel.msaweight.eslMSAWEIGHT_FILT_ORIGORDER,
    "random": libeasel.msaweight.eslMSAWEIGHT_FILT_RANDOM,
}

# --- Alphabet ---------------------------------------------------------------

cdef class Alphabet:
    """A biological alphabet, including additional marker symbols.

    This type is used to share an alphabet to several objects in the `easel`
    and `plan7` modules. Reference counting helps sharing the same instance
    everywhere, instead of reallocating memory every time an alphabet is
    needed.

    Use the factory class methods to obtain a default `Alphabet` for one of
    the three standard biological alphabets::

        >>> dna = Alphabet.dna()
        >>> rna = Alphabet.rna()
        >>> aa  = Alphabet.amino()

    """

    # --- Default constructors -----------------------------------------------

    cdef int _init_default(self, int ty) except 1 nogil:
        if self._abc != NULL:
            libeasel.alphabet.esl_alphabet_Destroy(self._abc)
        self._abc = libeasel.alphabet.esl_alphabet_Create(ty)
        if not self._abc:
            raise AllocationError("ESL_ALPHABET", sizeof(ESL_ALPHABET))
        return 0

    @classmethod
    def amino(cls):
        """Create a default amino-acid alphabet.
        """
        cdef Alphabet alphabet = Alphabet.__new__(Alphabet)
        alphabet._init_default(libeasel.alphabet.eslAMINO)
        return alphabet

    @classmethod
    def dna(cls):
        """Create a default DNA alphabet.
        """
        cdef Alphabet alphabet = Alphabet.__new__(Alphabet)
        alphabet._init_default(libeasel.alphabet.eslDNA)
        return alphabet

    @classmethod
    def rna(cls):
        """Create a default RNA alphabet.
        """
        cdef Alphabet alphabet = Alphabet.__new__(Alphabet)
        alphabet._init_default(libeasel.alphabet.eslRNA)
        return alphabet

    def __init__(self):
        raise TypeError("Cannot instantiate an alphabet directly")

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._abc = NULL

    def __dealloc__(self):
        libeasel.alphabet.esl_alphabet_Destroy(self._abc)

    def __repr__(self):
        assert self._abc != NULL

        cdef type   ty   = type(self)
        cdef str    name = ty.__name__

        if self._abc.type == libeasel.alphabet.eslRNA:
            return f"{name}.rna()"
        elif self._abc.type == libeasel.alphabet.eslDNA:
            return f"{name}.dna()"
        elif self._abc.type == libeasel.alphabet.eslAMINO:
            return f"{name}.amino()"
        else:
            return "{}({!r}, K={!r}, Kp={!r})".format(
                name,
                self._abc.sym.decode('ascii'),
                self._abc.K,
                self._abc.Kp
            )

    def __eq__(self, object other):
        assert self._abc != NULL
        # TODO: Update when we implement custom alphabet creation from Python.
        if isinstance(other, Alphabet):
            return self._eq(<Alphabet> other)
        return NotImplemented

    def __reduce__(self):
        assert self._abc != NULL
        if self._abc.type == libeasel.alphabet.eslRNA:
            return Alphabet.rna, ()
        elif self._abc.type == libeasel.alphabet.eslDNA:
            return Alphabet.dna, ()
        elif self._abc.type == libeasel.alphabet.eslAMINO:
            return Alphabet.amino, ()
        else:
            raise NotImplementedError("Alphabet.__reduce__")

    def __sizeof__(self):
        assert self._abc != NULL
        return libeasel.alphabet.esl_alphabet_Sizeof(self._abc) + sizeof(self)

    # --- Properties ---------------------------------------------------------

    @property
    def K(self):
        """`int`: The alphabet size, counting only actual alphabet symbols.

        Example:
            >>> Alphabet.dna().K
            4
            >>> Alphabet.amino().K
            20

        """
        assert self._abc != NULL
        return self._abc.K

    @property
    def Kp(self):
        """`int`: The complete alphabet size, including marker symbols.

        Example:
            >>> Alphabet.dna().Kp
            18
            >>> Alphabet.amino().Kp
            29

        """
        assert self._abc != NULL
        return self._abc.Kp

    @property
    def symbols(self):
        """`str`: The symbols composing the alphabet.

        Example:
            >>> Alphabet.dna().symbols
            'ACGT-RYMKSWHBVDN*~'
            >>> Alphabet.rna().symbols
            'ACGU-RYMKSWHBVDN*~'

        """
        assert self._abc != NULL
        return PyUnicode_DecodeASCII(self._abc.sym, self._abc.Kp, NULL)

    @property
    def type(self):
        """`str`: The alphabet type, as a short string.

        Example:
            >>> Alphabet.dna().type
            'DNA'
            >>> Alphabet.amino().type
            'amino'

        .. versionadded:: 0.8.2

        """
        assert self._abc != NULL
        return libeasel.alphabet.esl_abc_DecodeType(self._abc.type).decode('ascii')

    @property
    def gap_symbol(self):
        """`str`: The alphabet gap symbol.

        .. versionadded:: 0.11.1

        """
        assert self._abc != NULL
        return chr(libeasel.alphabet.esl_abc_CGetGap(self._abc))

    @property
    def gap_index(self):
        """`int`: The alphabet gap index.

        .. versionadded:: 0.11.1

        """
        assert self._abc != NULL
        return libeasel.alphabet.esl_abc_XGetGap(self._abc)

    # --- Utils --------------------------------------------------------------

    cdef inline bint _eq(self, Alphabet other) nogil:
        return self._abc.type == other._abc.type

    # --- Methods ------------------------------------------------------------

    cpdef bint is_dna(self):
        """Check whether the `Alphabet` object is a DNA alphabet.
        """
        assert self._abc != NULL
        return self._abc.type == libeasel.alphabet.eslDNA

    cpdef bint is_rna(self):
        """Check whether the `Alphabet` object is a RNA alphabet.
        """
        assert self._abc != NULL
        return self._abc.type == libeasel.alphabet.eslRNA

    cpdef bint is_amino(self):
        """Check whether the `Alphabet` object is a protein alphabet.
        """
        assert self._abc != NULL
        return self._abc.type == libeasel.alphabet.eslAMINO

    cpdef bint is_nucleotide(self):
        """Check whether the `Alphabet` object is a nucleotide alphabet.
        """
        assert self._abc != NULL
        return (
             self._abc.type == libeasel.alphabet.eslDNA
          or self._abc.type == libeasel.alphabet.eslRNA
        )

    cpdef VectorU8 encode(self, str sequence):
        """Encode a raw text sequence into its digital representation.

        Arguments:
            sequence (`str`): A raw sequence in text format.

        Returns:
            `~pyhmmer.easel.VectorU8`: A raw sequence in digital format.

        Example:
            >>> alphabet = easel.Alphabet.dna()
            >>> alphabet.encode("ACGT")
            VectorU8([0, 1, 2, 3])

        .. versionadded:: 0.6.3

        """
        assert self._abc != NULL

        cdef size_t   i
        cdef Py_UCS4  c
        cdef size_t   length  = len(sequence)
        cdef int      kind    = PyUnicode_KIND(sequence)
        cdef void*    data    = PyUnicode_DATA(sequence)
        cdef VectorU8 encoded = VectorU8.zeros(length)
        cdef uint8_t* buffer  = <uint8_t*> encoded._data

        for i in range(length):
            c = PyUnicode_READ(kind, data, i)
            if libeasel.alphabet.esl_abc_CIsValid(self._abc, c):
                buffer[i] = libeasel.alphabet.esl_abc_DigitizeSymbol(self._abc, c)
            else:
                raise ValueError(f"Invalid alphabet character in text sequence: {c}")

        return encoded

    cpdef str decode(self, const libeasel.ESL_DSQ[::1] sequence):
        """Decode a raw digital sequence into its textual representation.

        Arguments:
            sequence (`object`, *buffer-like*): A raw sequence in digital
                format. Any object implementing the buffer protocol (like
                `bytearray`, `~pyhmmer.easel.VectorU8`, etc.) may be given.

        Returns:
            `str`: A raw sequence in textual format.

        Example:
            >>> alphabet = easel.Alphabet.amino()
            >>> dseq = easel.VectorU8([0, 4, 2, 17, 3, 13, 0, 0, 5])
            >>> alphabet.decode(dseq)
            'AFDVEQAAG'

        .. versionadded:: 0.6.3

        """
        assert self._abc != NULL

        cdef libeasel.ESL_DSQ x
        cdef size_t           i
        cdef object           decoded
        cdef int              kind
        cdef char*            data
        cdef size_t           length  = sequence.shape[0]

        # NB(@althonos): Compatibility code for PyPy 3.6, which does
        #                not support directly writing to a string. Remove
        #                when Python 3.6 support is dropped.
        if SYS_VERSION_INFO_MAJOR <= 3 and SYS_VERSION_INFO_MINOR <= 7 and SYS_IMPLEMENTATION_NAME == "pypy":
            decoded = PyBytes_FromStringAndSize(NULL, length)
            data    = PyBytes_AsString(decoded)

            with nogil:
                for i in range(length):
                    x = sequence[i]
                    if libeasel.alphabet.esl_abc_XIsValid(self._abc, x):
                        data[i] = self._abc.sym[x]
                    else:
                        raise ValueError(f"Invalid alphabet character in digital sequence: {x}")

            return decoded.decode('ascii')

        else:
            decoded = PyUnicode_New(length, 0x7F)
            kind    = PyUnicode_KIND(decoded)
            data    = <char*> PyUnicode_DATA(decoded)

            for i in range(length):
                x = sequence[i]
                if libeasel.alphabet.esl_abc_XIsValid(self._abc, x):
                    PyUnicode_WRITE(kind, data, i, self._abc.sym[x])
                else:
                    raise ValueError(f"Invalid alphabet character in digital sequence: {x}")

            return decoded


# --- GeneticCode ------------------------------------------------------------

cdef class GeneticCode:
    """A genetic code table for translation.

    .. versionadded:: 0.7.2

    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._gcode = NULL

    def __dealloc__(self):
        libeasel.gencode.esl_gencode_Destroy(self._gcode)

    def __init__(
        self,
        int translation_table = 1,
        *,
        Alphabet nucleotide_alphabet not None = Alphabet.dna(),
        Alphabet amino_alphabet not None = Alphabet.amino(),
    ):
        """__init__(self, translation_table=1, *, nucleotide_alphabet=None, amino_alphabet=None)\n--\n

        Create a new genetic code for translating nucleotide sequences.

        Arguments:
            translation_table (`int`): The translation table to use. Check the
                `Wikipedia <https://w.wiki/47wo>`_ page listing all genetic
                codes for the available values.
            nucleotide_alphabet (`~pyhmmer.easel.Alphabet`): The nucleotide
                alphabet from which to translate the sequence.
            amino_alphabet (`~pyhmmer.easel.Alphabet`): The target alphabet
                into which to translate the sequence.

        """
        cdef int status

        if not nucleotide_alphabet.is_nucleotide():
            raise InvalidParameter("nucleotide_alphabet", nucleotide_alphabet, "nucleotide alphabet")
        if not amino_alphabet.is_amino():
            raise InvalidParameter("amino_alphabet", amino_alphabet, "amino alphabet")

        self._gcode = libeasel.gencode.esl_gencode_Create(nucleotide_alphabet._abc, amino_alphabet._abc)
        if self._gcode == NULL:
            raise AllocationError("ESL_GENCODE", sizeof(ESL_GENCODE))

        self.amino_alphabet = amino_alphabet
        self.nucleotide_alphabet = nucleotide_alphabet
        self.translation_table = translation_table

    def __repr__(self):
        cdef str ty = type(self).__name__
        if self.nucleotide_alphabet.is_dna():
            return f"{ty}({self.translation_table!r})"
        else:
            return f"{ty}({self.translation_table!r}, nucleotide_alphabet={self.nucleotide_alphabet!r})"

    def __reduce__(self):
        constructor = functools.partial(
            type(self),
            translation_table=self.translation_table,
            nucleotide_alphabet=self.nucleotide_alphabet,
            amino_alphabet=self.amino_alphabet
        )
        return constructor, ()

    # --- Properties ---------------------------------------------------------

    @property
    def translation_table(self):
        """`int`: The translation table in use.

        Can be set manually to a different number to change the
        translation table for the current `GeneticCode` object.

        """
        assert self._gcode != NULL
        return self._gcode.transl_table

    @translation_table.setter
    def translation_table(self, int translation_table):
        assert self._gcode != NULL
        status = libeasel.gencode.esl_gencode_Set(self._gcode, translation_table)
        if status == libeasel.eslENOTFOUND:
            raise InvalidParameter("translation_table", translation_table, hint="translation table code")
        elif status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_gencode_Set")

    @property
    def description(self):
        """`str`: A description of the translation table currently in use.
        """
        assert self._gcode != NULL
        return self._gcode.desc.decode('ascii')

    # --- Utils --------------------------------------------------------------

    cdef int _translate(
        self,
        const ESL_DSQ* seq,
        int64_t seqlen,
        ESL_DSQ* out,
        int64_t outlen
    ) except -1 nogil:
        cdef int     aa
        cdef int64_t i
        cdef int64_t j

        if seqlen % 3 != 0:
            raise ValueError(f"Invalid sequence of length {seqlen!r}")
        if outlen < seqlen // 3:
            raise BufferError(f"Output buffer too short for sequence of length {seqlen // 3!r}")

        for i, j in enumerate(range(0, seqlen, 3)):
            aa = libeasel.gencode.esl_gencode_GetTranslation(self._gcode, <ESL_DSQ*> &seq[j])
            if aa == -1:
                raise ValueError(f"Failed to translate codon at index {j!r}")
            out[i] = aa

        return 0

    # --- Methods ------------------------------------------------------------

    cpdef VectorU8 translate(self, const ESL_DSQ[::1] sequence):
        """Translate a raw nucleotide sequence into a protein.

        Arguments:
            sequence (`object`, *buffer-like*): A raw sequence in digital
                format. Any object implementing the buffer protocol (like
                `bytearray`, `~pyhmmer.easel.VectorU8`, etc.) may be given.

        Returns:
            `~pyhmmer.easel.VectorU8`: The translation of the input
            sequence, as a raw digital sequence.

        Raises:
            `ValueError`: When ``sequence`` could not be translated
                properly, because of a codon could not be recognized, or
                because the sequence has an invalid length.

        Note:
            The translation of a DNA/RNA codon supports ambiguous codons.
            If the amino acid is unambiguous, despite codon ambiguity,
            the correct amino acid is still determined: ``GGR`` translates
            as ``Gly``, ``UUY`` as ``Phe``, etc. If there is no single
            unambiguous amino acid translation, the codon is translated
            as ``X``. Ambiguous amino acids (such as ``J`` or ``B``) are
            never produced.

        """
        cdef int64_t  nlen = sequence.shape[0]
        cdef int64_t  alen = nlen // 3
        cdef VectorU8 prot = VectorU8.zeros(alen)

        if sequence.shape[0] > 0:
            with nogil:
                self._translate(&sequence[0], nlen, <ESL_DSQ*> prot._data, alen)

        return prot


# --- Bitfield ---------------------------------------------------------------

cdef class Bitfield:
    """A statically sized sequence of booleans stored as a packed bitfield.

    Example:
        Instantiate a bitfield from an iterable, where each object will be
        tested for truth:

            >>> bitfield = Bitfield([True, False, False])
            >>> len(bitfield)
            3
            >>> bitfield[0]
            True
            >>> bitfield[1]
            False

        Use `Bitfield.zeros` and `Bitfield.ones` to initialize a bitfield of
        a given length with all fields set to :math:`0` or :math:`1`::

            >>> Bitfield.zeros(4)
            Bitfield([False, False, False, False])
            >>> Bitfield.ones(4)
            Bitfield([True, True, True, True])

        Use indexing to access and edit individual bits::

            >>> bitfield[0] = True
            >>> bitfield[0]
            True
            >>> bitfield[0] = False
            >>> bitfield[0]
            False

    """

    # --- Class methods ------------------------------------------------------

    @classmethod
    def _from_raw_bytes(cls, object buffer, int n, str byteorder):
        f"""Create a new bitfield using the given bytes to fill its contents.
        """
        cdef const uint8_t[::1] bytes
        cdef Bitfield           bitfield = cls.zeros(n)
        cdef object             view     = memoryview(buffer)

        # fix endianness if needed
        if byteorder != SYS_BYTEORDER:
            newbuffer = array.array("Q")
            newbuffer.frombytes(view)
            newbuffer.byteswap()
            view = memoryview(newbuffer)

        # assign the items
        bytes = view.cast("B")
        assert bytes.shape[0] == bitfield._shape[0] * sizeof(uint64_t)
        if n > 0:
            with nogil:
                memcpy(bitfield._b.b, &bytes[0], bitfield._shape[0]*sizeof(uint64_t))
        return bitfield

    @classmethod
    def zeros(cls, size_t n):
        """Create a new bitfield of size ``n`` with all elements set to `False`.

        .. versionadded:: 0.7.0

        """
        if n <= 0:
            raise ValueError("Cannot create an empty `Bitfield`")
        cdef Bitfield bitfield = Bitfield.__new__(Bitfield)
        bitfield._shape[0] = (n + 63) // 64
        bitfield._b = libeasel.bitfield.esl_bitfield_Create(n)
        if not bitfield._b:
            raise AllocationError("ESL_BITFIELD", sizeof(ESL_BITFIELD))
        return bitfield

    @classmethod
    def ones(cls, size_t n):
        """Create a new bitfield of size ``n`` with all elements set to `True`.

        .. versionadded:: 0.7.0

        """
        cdef Bitfield bitfield = cls.zeros(n)
        with nogil:
            memset(bitfield._b.b, 0xFF, bitfield._shape[0]*sizeof(uint64_t))
        return bitfield

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._b = NULL
        self._shape[0] = 0

    def __dealloc__(self):
        libeasel.bitfield.esl_bitfield_Destroy(self._b)

    def __init__(self, object iterable):
        """__init__(self, iterable)\n--\n

        Create a new bitfield from an iterable of objects.

        Objects yielded by the iterable can be of any type and will be
        tested for truth before setting the corresponding field.

        Raises:
            `ValueError`: When given an empty iterable.

        """
        if not isinstance(iterable, collections.abc.Sized):
            iterable = list(iterable)

        cdef size_t n = len(iterable)
        if n <= 0:
            raise ValueError("Cannot create an empty `Bitfield`")

        self._shape[0] = (n + 63) // 64
        self._b = libeasel.bitfield.esl_bitfield_Create(n)
        if not self._b:
            raise AllocationError("ESL_BITFIELD", sizeof(ESL_BITFIELD))

        for i, item in enumerate(iterable):
            if item:
                libeasel.bitfield.esl_bitfield_Set(self._b, i)

    def __len__(self):
        assert self._b != NULL
        return self._b.nb

    def __getitem__(self, int idx):
        assert self._b != NULL
        cdef size_t index_ = self._wrap_index(idx)
        return libeasel.bitfield.esl_bitfield_IsSet(self._b, index_)

    def __setitem__(self, index, value):
        assert self._b != NULL
        cdef size_t index_ = self._wrap_index(index)
        if value:
            libeasel.bitfield.esl_bitfield_Set(self._b, index_)
        else:
            libeasel.bitfield.esl_bitfield_Clear(self._b, index_)

    def __eq__(self, object other):
        assert self._b != NULL

        cdef size_t        nu      = (self._b.nb // 64) + (self._b.nb % 64 != 0)
        cdef ESL_BITFIELD* other_b
        cdef int           cmp

        if isinstance(other, Bitfield):
            other_b = (<Bitfield> other)._b
            assert other_b != NULL
            if self._b.nb != other_b.nb:
                return False
            with nogil:
                cmp = memcmp(
                    <const void*> self._b.b,
                    <const void*> other_b.b,
                    nu
                )
            return cmp == 0

        return NotImplemented

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"{ty}({list(self)!r})"

    def __reduce_ex__(self, int protocol):
        assert self._b != NULL

        # use out-of-band pickling (only supported since protocol 5, see
        # https://docs.python.org/3/library/pickle.html#out-of-band-buffers)
        if protocol >= 5:
            buffer = pickle.PickleBuffer(self)
        else:
            buffer = array.array("Q")
            buffer.frombytes(memoryview(self).cast("B"))

        return self._from_raw_bytes, (buffer, self._b.nb, SYS_BYTEORDER)

    def __getstate__(self):
        assert self._b != NULL

        cdef size_t nu    = (self._b.nb // 64) + (self._b.nb % 64 != 0)
        cdef object mview = PyMemoryView_FromMemory(<char*> self._b.b, nu * sizeof(uint64_t), PyBUF_READ)
        cdef object b     = array.array("Q")

        b.frombytes(mview)

        return {"nb": self._b.nb, "b": b}

    def __setstate__(self, state):
        cdef size_t        nb = state["nb"]
        cdef size_t        nu = (nb // 64) + (nb % 64 != 0)
        cdef uint64_t[::1] b  = state["b"]

        if nb <= 0:
            raise ValueError("Cannot create an empty `Bitfield`")

        if self._b == NULL:
            self._b = libeasel.bitfield.esl_bitfield_Create(nb)
        else:
            self._b.nb = nb
            self._b.b  = <uint64_t*> realloc(self._b.b, nu * sizeof(uint64_t))

        with nogil:
            memcpy(self._b.b, &b[0], nu * sizeof(uint64_t))

    def __sizeof__(self):
        assert self._b != NULL
        cdef size_t nu = (self._b.nb // 64) + (self._b.nb % 64 != 0)
        return sizeof(uint64_t) * nu + sizeof(ESL_BITFIELD) + sizeof(self)

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def __getbuffer__(self, Py_buffer* buffer, int flags):
        assert self._b != NULL

        if flags & PyBUF_FORMAT:
            buffer.format = b"Q"
        else:
            buffer.format = NULL
        buffer.buf = self._b.b
        buffer.internal = NULL
        buffer.itemsize = sizeof(uint64_t)
        buffer.len = self._shape[0] * sizeof(uint64_t)
        buffer.ndim = 1
        buffer.obj = self
        buffer.readonly = 0
        buffer.shape = self._shape
        buffer.suboffsets = NULL
        buffer.strides = NULL

    # --- Utils --------------------------------------------------------------

    cdef size_t _wrap_index(self, int index) except -1:
        if index < 0:
            index += self._b.nb
        if index >= self._b.nb or index < 0:
            raise IndexError("bitfield index out of range")
        return <size_t> index

    # --- Methods ------------------------------------------------------------

    cpdef size_t count(self, bint value=1):
        """Count the number occurrences of ``value`` in the bitfield.

        If no argument is given, counts the number of `True` occurences.

        Example:
            >>> bitfield = Bitfield.zeros(8)
            >>> bitfield.count(False)
            8
            >>> bitfield[0] = bitfield[1] = True
            >>> bitfield.count()
            2

        """
        assert self._b != NULL
        cdef size_t count_
        with nogil:
            count_ = libeasel.bitfield.esl_bitfield_Count(self._b)
        return count_ if value else self._b.nb - count_

    cpdef Bitfield copy(self):
        """Return a copy of this bitfield object.

        .. versionadded:: 0.7.0

        """
        assert self._b != NULL
        cdef Bitfield copy = type(self).zeros(self._b.nb)
        with nogil:
            memcpy(copy._b.b, self._b.b, self._shape[0]*sizeof(uint64_t))
        return copy

    cpdef void toggle(self, int index) except *:
        """Switch the value of one single bit.

        Example:
            >>> bitfield = Bitfield.zeros(8)
            >>> bitfield[0]
            False
            >>> bitfield.toggle(0)
            >>> bitfield[0]
            True
            >>> bitfield.toggle(0)
            >>> bitfield[0]
            False

        """
        assert self._b != NULL
        cdef size_t index_ = self._wrap_index(index)
        with nogil:
            libeasel.bitfield.esl_bitfield_Toggle(self._b, index_)


# --- KeyHash ----------------------------------------------------------------

cdef class KeyHash:
    """A dynamically resized container to store byte keys using a hash table.

    Internally uses Bob Jenkins' *one at a time* hash, a simple and
    efficient hash function published in 1997 that exhibits
    `avalanche <https://en.wikipedia.org/wiki/Avalanche_effect>`_
    behaviour.

    Example:
        Add new keys to the key hash using the `~KeyHash.add` method
        like you would with a Python `set`::

            >>> kh = KeyHash()
            >>> kh.add(b"key")
            0

        Check if a key hash contains a given key::

            >>> b"key" in kh
            True
            >>> b"missing" in kh
            False

        Get the index associated with a key using the indexing notation::

            >>> kh[b"key"]
            0
            >>> kh[b"missing"]
            Traceback (most recent call last):
              ...
            KeyError: b'missing'

        Iterate over the keys of the key hash, in the order of insertion::

            >>> kh.add(b"key2")
            1
            >>> for k in kh:
            ...     print(k)
            b'key'
            b'key2'

    See Also:
        The Wikipedia article for Bob Jenkins' hash functions:
        https://en.wikipedia.org/wiki/Jenkins_hash_function

    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._kh = NULL

    def __dealloc__(self):
        libeasel.keyhash.esl_keyhash_Destroy(self._kh)

    def __init__(self):
        """__init__(self)\n--\n

        Create a new empty key-hash collection.

        """
        cdef int status
        with nogil:
            if self._kh == NULL:
                self._kh = libeasel.keyhash.esl_keyhash_Create()
            else:
                status = libeasel.keyhash.esl_keyhash_Reuse(self._kh)
                if status != libeasel.eslOK:
                    raise UnexpectedError(status, "esl_keyhash_Reuse")
        if self._kh == NULL:
            raise AllocationError("ESL_KEYHASH", sizeof(ESL_KEYHASH))

    def __copy__(self):
        return self.copy()

    def __len__(self):
        assert self._kh != NULL
        return libeasel.keyhash.esl_keyhash_GetNumber(self._kh)

    def __contains__(self, object value):
        assert self._kh != NULL

        if not isinstance(value, bytes):
            return False

        cdef       int    status
        cdef const char*  key    = value
        cdef       size_t length = len(value)

        with nogil:
            status = libeasel.keyhash.esl_keyhash_Lookup(self._kh, key, length, NULL)
        if status == libeasel.eslOK:
            return True
        elif status == libeasel.eslENOTFOUND:
            return False
        else:
            raise UnexpectedError(status, "esl_keyhash_Lookup")

    def __getitem__(self, bytes item):
        assert self._kh != NULL

        cdef       int    status
        cdef       int    index
        cdef const char*  key    = item
        cdef       size_t length = len(item)

        with nogil:
            status = libeasel.keyhash.esl_keyhash_Lookup(self._kh, key, length, &index)
        if status == libeasel.eslOK:
            return index
        elif status == libeasel.eslENOTFOUND:
            raise KeyError(item)
        else:
            raise UnexpectedError(status, "esl_keyhash_Lookup")

    def __iter__(self):
        assert self._kh != NULL

        cdef int   i
        cdef int   offset
        cdef char* key

        for i in range(libeasel.keyhash.esl_keyhash_GetNumber(self._kh)):
            offset = self._kh.key_offset[i]
            key = &self._kh.smem[offset]
            yield <bytes> key

    def __getstate__(self):
        assert self._kh != NULL

        cdef ssize_t i

        cdef object  smem       = array.array("b")
        cdef object  hashtable  = array.array("i")
        cdef object  key_offset = array.array("i")
        cdef object  nxt        = array.array("i")
        cdef object  sview      = PyMemoryView_FromMemory(<char*> self._kh.smem,       self._kh.salloc   * sizeof(char), PyBUF_READ)
        cdef object  hview      = PyMemoryView_FromMemory(<char*> self._kh.hashtable,  self._kh.hashsize * sizeof(int),  PyBUF_READ)
        cdef object  kview      = PyMemoryView_FromMemory(<char*> self._kh.key_offset, self._kh.nkeys    * sizeof(int),  PyBUF_READ)
        cdef object  nview      = PyMemoryView_FromMemory(<char*> self._kh.nxt,        self._kh.nkeys    * sizeof(int),  PyBUF_READ)

        smem.frombytes(sview)
        hashtable.frombytes(hview)
        key_offset.frombytes(kview)
        nxt.frombytes(nview)

        return {
            "hashtable": hashtable,
            "hashsize": self._kh.hashsize,
            "key_offset": key_offset,
            "nxt": nxt,
            "nkeys": self._kh.nkeys,
            "kalloc": self._kh.kalloc,
            "smem": smem,
            "salloc": self._kh.salloc,
            "sn": self._kh.sn,
        }

    def __setstate__(self, state):
        cdef size_t    i
        cdef char[::1] smem       = state["smem"]
        cdef int[::1]  hashtable  = state["hashtable"]
        cdef int[::1]  key_offset = state["key_offset"]
        cdef int[::1]  nxt        = state["nxt"]

        assert smem.shape[0] <= state["salloc"]
        assert hashtable.shape[0] <= state["hashsize"]
        assert key_offset.shape[0] <= state["nkeys"]
        assert nxt.shape[0] <= state["nkeys"]

        # FIXME: avoid reallocation if possible?
        if self._kh != NULL:
            libeasel.keyhash.esl_keyhash_Destroy(self._kh)
        # allocate the keyhash using the right allocation sizes
        self._kh = libeasel.keyhash.esl_keyhash_CreateCustom(
            state["hashsize"],
            state["kalloc"],
            state["salloc"]
        )
        if not self._kh:
            raise AllocationError("ESL_KEYHASH", sizeof(ESL_KEYHASH))

        # copy numeric values
        self._kh.sn = state["sn"]
        self._kh.nkeys = state["nkeys"]

        # copy data
        with nogil:
            memcpy(self._kh.smem,       &smem[0],       self._kh.salloc   * sizeof(char))
            memcpy(self._kh.hashtable,  &hashtable[0],  self._kh.hashsize * sizeof(int))
            memcpy(self._kh.key_offset, &key_offset[0], self._kh.nkeys    * sizeof(int))
            memcpy(self._kh.nxt,        &nxt[0],        self._kh.nkeys    * sizeof(int))

    def __sizeof__(self):
        assert self._kh != NULL
        return (
            sizeof(int) * self._kh.hashsize   # kh->hashtable
          + sizeof(int) * self._kh.kalloc * 2 # kh->key_offset, kh->nxt
          + sizeof(int) * self._kh.salloc     # kh->smem
          + sizeof(ESL_KEYHASH)
          + sizeof(self)
        )

    # --- Methods ------------------------------------------------------------

    cpdef int add(self, bytes key) except -1:
        """Add a new key to the hash table, and return its index.

        If ``key`` was already in the hash table, the previous index is
        returned::

            >>> kh = KeyHash()
            >>> kh.add(b"first")
            0
            >>> kh.add(b"second")
            1
            >>> kh.add(b"first")
            0

        Arguments:
            key (`bytes`): The key to add to the hash table.

        Returns:
            `int`: The index corresponding to the added ``key``.

        .. versionadded:: 0.3.0

        """
        assert self._kh != NULL

        cdef       int    status
        cdef       int    index
        cdef const char*  k      = key
        cdef       size_t length = len(key)

        with nogil:
            status = libeasel.keyhash.esl_keyhash_Store(self._kh, k, length, &index)
        if status == libeasel.eslOK or status == libeasel.eslEDUP:
            return index
        else:
            raise UnexpectedError(status, "esl_keyhash_Store")

    cpdef void clear(self) except *:
        """Remove all entries from the collection.
        """
        cdef int status
        with nogil:
            status = libeasel.keyhash.esl_keyhash_Reuse(self._kh)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_keyhash_Reuse")

    cpdef KeyHash copy(self):
        """Create and return an exact copy of this mapping.

        Example:
            >>> kh = KeyHash()
            >>> kh.add(b"key")
            0
            >>> copy = kh.copy()
            >>> b"key" in copy
            True

        """
        assert self._kh != NULL
        cdef KeyHash new = KeyHash.__new__(KeyHash)
        with nogil:
            new._kh = libeasel.keyhash.esl_keyhash_Clone(self._kh)
        if not new._kh:
            raise AllocationError("ESL_KEYHASH", sizeof(ESL_KEYHASH))
        return new


# --- Matrix & Vector --------------------------------------------------------

@cython.no_gc_clear
cdef class Vector:
    """An abstract 1D array of fixed size.

    .. versionadded:: 0.4.0

    """

    # --- Class methods ------------------------------------------------------

    @classmethod
    def zeros(cls, int n):
        """Create a vector of size ``n`` filled with zeros.
        """
        cdef Vector vec = cls()
        vec._allocate(n)
        return vec

    @classmethod
    def _from_raw_bytes(cls, object buffer, int n, str byteorder):
        f"""Create a new vector using the given bytes to fill its contents.
        """
        cdef const uint8_t[::1] bytes
        cdef Vector             vec      = cls.zeros(n)
        cdef size_t             itemsize = vec.itemsize
        cdef object             view     = memoryview(buffer)

        # fix endianness if needed
        if byteorder != SYS_BYTEORDER and vec.itemsize > 1:
            newbuffer = array.array(vec.format)
            newbuffer.frombytes(view)
            newbuffer.byteswap()
            view = memoryview(newbuffer)

        # assign the items
        bytes = view.cast("B")
        assert bytes.shape[0] == n * vec.itemsize
        if n > 0:
            with nogil:
                memcpy(vec._data, &bytes[0], n * itemsize)
        return vec

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._owner = None
        self._n = 0
        self._shape = (0,)
        self._data = NULL

    def __dealloc__(self):
        if self._owner is None and self._data != NULL:
            free(self._data)
        self._data = NULL

    def __init__(self, object iterable = ()):
        """__init__(self, iterable=())\n--\n

        Create a new vector from the given iterable of values.

        """
        raise TypeError("Can't instantiate abstract class 'Vector'")

    def __bool__(self):
        return self._n != 0

    def __len__(self):
        return self._n

    def __copy__(self):
        return self.copy()

    def __getbuffer__(self, Py_buffer* buffer, int flags):
        assert self._data != NULL

        if flags & PyBUF_FORMAT:
            buffer.format = <char*> self._format()
        else:
            buffer.format = NULL
        buffer.buf = self._data
        buffer.internal = NULL
        buffer.itemsize = self.itemsize
        buffer.len = self._n * self.itemsize
        buffer.ndim = 1
        buffer.obj = self
        buffer.readonly = 0
        buffer.shape = self._shape
        buffer.suboffsets = NULL

        if SYS_IMPLEMENTATION_NAME == "pypy":
            buffer.internal = PyMem_Malloc(sizeof(Py_ssize_t))
            if buffer.internal == NULL:
                raise AllocationError("ssize_t", sizeof(Py_ssize_t))
            buffer.strides = <Py_ssize_t*> buffer.internal
            buffer.strides[0] = self.strides[0]
        else:
            buffer.strides = NULL

    def __releasebuffer__(self, Py_buffer* buffer):
        if SYS_IMPLEMENTATION_NAME == "pypy":
            PyMem_Free(buffer.internal)
            buffer.internal = NULL

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"{ty}({list(self)!r})"

    def __sizeof__(self):
        return self._n * self.itemsize + sizeof(self)

    def __reduce_ex__(self, int protocol):
        assert self._data != NULL

        cdef object buffer

        # use out-of-band pickling (only supported since protocol 5, see
        # https://docs.python.org/3/library/pickle.html#out-of-band-buffers)
        if protocol >= 5:
            buffer = pickle.PickleBuffer(self)
        else:
            buffer = array.array(self.format)
            buffer.frombytes(memoryview(self).cast("b"))

        return self._from_raw_bytes, (buffer, self._n, SYS_BYTEORDER)

    def __add__(Vector self, object other):
        assert self._data != NULL
        cdef Vector new = self.copy()
        return new.__iadd__(other)

    def __sub__(Vector self, object other):
        assert self._data != NULL
        cdef Vector new = self.copy()
        return new.__isub__(other)

    def __mul__(Vector self, object other):
        assert self._data != NULL
        cdef Vector new = self.copy()
        return new.__imul__(other)

    # --- Properties ---------------------------------------------------------

    @property
    def shape(self):
        """`tuple`: The shape of the vector.
        """
        return tuple(self._shape)

    @property
    def strides(self):
        """`tuple`: The strides of the vector.
        """
        return (self.itemsize,)

    @property
    def itemsize(self):
        """`int`: The size of each item in the vector, in bytes.

        .. versionadded:: 0.4.6

        """

    @property
    def format(self):
        """`str`: The format of each item in the vector.

        See Also:
            The `array` module of the Python standard library for a detail
            about available type codes.

        .. versionadded:: 0.4.6

        """

    # --- Utility ------------------------------------------------------------

    cdef const char* _format(self) noexcept:
        return NULL

    cdef int _allocate(self, size_t n) except -1:
        # NB(@althonos): malloc and calloc are not guaranteed to return a
        #                pointer when called with a null allocation size,
        #                so we allocate a single item instead
        cdef int n_alloc  = 1 if n == 0 else n
        cdef int itemsize = self.itemsize

        if self._data != NULL:
            free(self._data)

        self._n = self._shape[0] = n

        with nogil:
            self._data = calloc(n_alloc, itemsize)
        if self._data == NULL:
            raise AllocationError("uint8_t", 1, self.itemsize * n_alloc)

        return 0

    # --- Methods ------------------------------------------------------------

    def argmax(self):
        """Return index of the maximum element in the vector.

        Raises:
            `ValueError`: When called on an empty vector.

        """

    def argmin(self):
        """Return index of the minimum element in the vector.

        Raises:
            `ValueError`: When called on an empty vector.

        """

    def copy(self):
        """Create a copy of the vector, allocating a new buffer.
        """

    def max(self):
        """Return value of the maximum element in the vector.

        Raises:
            `ValueError`: When called on an empty vector.

        """

    def min(self):
        """Return value of the minimum element in the vector.

        Raises:
            `ValueError`: When called on an empty vector.

        """

    def reverse(self):
        """Reverse the vector, in place.
        """

    def sum(self):
        """Returns the scalar sum of all elements in the vector.
        """


cdef class VectorF(Vector):
    """A vector storing single-precision floating point numbers.

    Individual elements of a vector can be accessed and modified with
    the usual indexing notation::

        >>> v = VectorF([1.0, 2.0, 3.0])
        >>> v[0]
        1.0
        >>> v[-1]
        3.0
        >>> v[0] = v[-1] = 4.0
        >>> v
        VectorF([4.0, 2.0, 4.0])

    Slices are also supported, and they do not copy data (use the
    `~pyhmmer.easel.VectorF.copy` method to allocate a new vector)::

        >>> v = VectorF(range(6))
        >>> v[2:5]
        VectorF([2.0, 3.0, 4.0])
        >>> v[2:-1] = 10.0
        >>> v
        VectorF([0.0, 1.0, 10.0, 10.0, 10.0, 5.0])

    Addition and multiplication is supported for scalars, in place or not::

        >>> v = VectorF([1.0, 2.0, 3.0])
        >>> v += 1
        >>> v
        VectorF([2.0, 3.0, 4.0])
        >>> v * 3
        VectorF([6.0, 9.0, 12.0])

    Pairwise operations can also be performed, but only on vectors of
    the same dimension and precision::

        >>> v = VectorF([1.0, 2.0, 3.0])
        >>> v * v
        VectorF([1.0, 4.0, 9.0])
        >>> v += VectorF([3.0, 4.0, 5.0])
        >>> v
        VectorF([4.0, 6.0, 8.0])
        >>> v *= VectorF([1.0])
        Traceback (most recent call last):
          ...
        ValueError: cannot pairwise multiply vectors of different sizes

    Objects of this type support the buffer protocol, and can be viewed
    as a `numpy.ndarray` of one dimension using the `numpy.asarray`
    function, and can be passed without copy to most `numpy` functions:

        >>> v = VectorF([1.0, 2.0, 3.0])
        >>> numpy.asarray(v)
        array([1., 2., 3.], dtype=float32)
        >>> numpy.log2(v)
        array([0.       , 1.       , 1.5849625], dtype=float32)

    .. versionadded:: 0.4.0

    """

    # --- Magic methods ------------------------------------------------------

    def __init__(self, object iterable = ()):
        """__init__(self, iterable=())\n--\n

        Create a new float vector from the given data.

        """
        cdef int        n
        cdef size_t     i
        cdef float      item
        cdef float[::1] view
        cdef int        n_alloc
        cdef float*     data

        # collect iterable if it's not a `Sized` object
        if not isinstance(iterable, collections.abc.Sized):
            iterable = array.array("f", iterable)
        n = len(iterable)

        # make sure __init__ is only called once
        if self._data != NULL:
            raise RuntimeError("Vector.__init__ must not be called more than once")
        # make sure the vector has a positive size
        if n < 0:
            raise ValueError("Cannot create a vector with negative size")

        # allocate vector storage
        self._allocate(n)

        # try to copy the memory quickly if *iterable* implements the buffer
        # protocol, otherwise fall back to cop
        data = <float*> self._data
        try:
            view = iterable
        except (TypeError, ValueError):
            for i, item in enumerate(iterable):
                data[i] = item
        else:
            with nogil:
                memcpy(data, &view[0], n * sizeof(float))

    def __eq__(self, object other):
        assert self._data != NULL

        cdef const float[:] buffer
        cdef int            i
        cdef int            status
        cdef const float*   data   = <const float*> self._data

        # check matrix type and dimensions
        try:
            buffer = other
        except ValueError:
            return NotImplemented
        if buffer.ndim != 1:
            return NotImplemented
        if buffer.shape[0] != self._n:
            return False
        elif self._n == 0:
            return True

        # check values
        with nogil:
            status = libeasel.vec.esl_vec_FCompare(&data[0], &buffer[0], self._n, 0)

        return status == libeasel.eslOK

    def __getitem__(self, object index):
        assert self._data != NULL

        cdef VectorF new
        cdef int     idx
        cdef ssize_t start
        cdef ssize_t stop
        cdef ssize_t step
        cdef float*  data  = <float*> self._data

        if isinstance(index, slice):
            start, stop, step = index.indices(self._n)
            if step != 1:
                raise ValueError(f"cannot slice a Vector with step other than 1")
            new = VectorF.__new__(VectorF)
            new._owner = self
            new._n = new._shape[0] = stop - start
            new._data = NULL if new._n == 0 else <void*> &(data[start])
            return new
        else:
            idx = index
            if idx < 0:
                idx += self._n
            if idx < 0 or idx >= self._n:
                raise IndexError("vector index out of range")
            return data[idx]

    def __setitem__(self, object index, float value):
        assert self._data != NULL

        cdef ssize_t x
        cdef float*  data = <float*> self._data

        if isinstance(index, slice):
            for x in range(*index.indices(self._n)):
                data[x] = value
        else:
            x = index
            if x < 0:
                x += self._n
            if x < 0  or x >= self._n:
                raise IndexError("vector index out of range")
            data[x] = value

    def __neg__(self):
        assert self._data != NULL

        cdef int     i
        cdef VectorF new  = self.copy()
        cdef float*  data = <float*> new._data

        with nogil:
            for i in range(self._n):
                data[i] = -data[i]

        return new

    def __iadd__(self, object other):
        assert self._data != NULL

        cdef VectorF      other_vec
        cdef const float* other_data
        cdef float        other_f
        cdef float*       data       = <float*> self._data

        if isinstance(other, VectorF):
            other_vec = other
            other_data = <float*> other_vec._data
            assert other_vec._data != NULL or other_vec._n == 0
            if self._n != other_vec._n:
                raise ValueError("cannot pairwise add vectors of different sizes")
            with nogil:
                libeasel.vec.esl_vec_FAdd(data, other_data, self._n)
        else:
            other_f = other
            with nogil:
                libeasel.vec.esl_vec_FIncrement(data, self._n, other_f)
        return self

    def __isub__(self, object other):
        assert self._data != NULL

        cdef int          i
        cdef VectorF      other_vec
        cdef const float* other_data
        cdef float        other_f
        cdef float*       data       = <float*> self._data

        if isinstance(other, VectorF):
            other_vec = other
            other_data = <float*> other_vec._data
            assert other_vec._data != NULL or other_vec._n == 0
            if self._n != other_vec._n:
                raise ValueError("cannot pairwise subtract vectors of different sizes")
            with nogil:
                for i in range(self._n):
                    data[i] -= other_data[i]
        else:
            other_f = other
            with nogil:
                for i in range(self._n):
                    data[i] -= other_f
        return self

    def __imul__(self, object other):
        assert self._data != NULL

        cdef int          i
        cdef VectorF      other_vec
        cdef const float* other_data
        cdef float        other_f
        cdef float*       data       = <float*> self._data

        if isinstance(other, VectorF):
            other_vec = other
            other_data = <float*> other_vec._data
            assert other_vec._data != NULL or other_vec._n == 0
            if self._n != other_vec._n:
                raise ValueError("cannot pairwise multiply vectors of different sizes")
            # NB(@althonos): There is no function in `vectorops.h` to do this
            # for now...
            with nogil:
                for i in range(self._n):
                    data[i] *= other_data[i]
        else:
            other_f = other
            with nogil:
                libeasel.vec.esl_vec_FScale(data, self._n, other_f)
        return self

    def __truediv__(VectorF self, object other):
        assert self._data != NULL
        cdef VectorF new = self.copy()
        return new.__itruediv__(other)

    def __itruediv__(self, object other):
        assert self._data != NULL

        cdef int          i
        cdef VectorF      other_vec
        cdef const float* other_data
        cdef float        other_f
        cdef float*       data       = <float*> self._data

        if isinstance(other, VectorF):
            other_vec = other
            other_data = <float*> other_vec._data
            assert other_vec._data != NULL or other_vec._n == 0
            if self._n != other_vec._n:
                raise ValueError("cannot pairwise divide vectors of different sizes")
            with nogil:
                for i in range(self._n):
                    data[i] /= other_data[i]
        else:
            other_f = other
            with nogil:
                for i in range(self._n):
                    data[i] /= other_f
        return self

    def __matmul__(VectorF self, object other):
        assert self._data != NULL

        cdef VectorF      other_vec
        cdef const float* other_data
        cdef const float* data       = <float*> self._data
        cdef float        res

        if not isinstance(other, VectorF):
            return NotImplemented

        other_vec = other
        other_data = <float*> other_vec._data
        assert other_data != NULL
        if self._n != other_vec._n:
            raise ValueError("cannot multiply vectors of different sizes")
        with nogil:
            res = libeasel.vec.esl_vec_FDot(data, other_data, self._n)
        return res

    # --- Properties ---------------------------------------------------------

    @property
    def itemsize(self):
        return sizeof(float)

    @property
    def format(self):
        return "f"

    # --- Utility ------------------------------------------------------------

    cdef const char* _format(self) noexcept:
        return b"f"

    # --- Methods ------------------------------------------------------------

    cpdef int argmax(self) except -1:
        assert self._data != NULL
        if self._n == 0:
            raise ValueError("argmax() called on an empty vector")
        with nogil:
            return libeasel.vec.esl_vec_FArgMax(<float*> self._data, self._n)

    cpdef int argmin(self) except -1:
        assert self._data != NULL
        if self._n == 0:
            raise ValueError("argmin() called on an empty vector")
        with nogil:
            return libeasel.vec.esl_vec_FArgMin(<float*> self._data, self._n)

    cpdef VectorF copy(self):
        assert self._data != NULL

        cdef VectorF new
        cdef int     n_alloc = 1 if self._n == 0 else self._n

        new = VectorF.__new__(VectorF)
        new._n = new._shape[0] = self._n

        new._data = calloc(n_alloc, sizeof(float))
        if new._data == NULL:
            raise AllocationError("float", sizeof(float), n_alloc)
        with nogil:
            memcpy(new._data, self._data, self._n * sizeof(float))

        return new

    cpdef float entropy(self) except *:
        """Compute the Shannon entropy of the vector.

        The Shannon entropy of a probability vector is defined as:

        .. math::

            H = \\sum_{i=0}^{N}{\\log_2 p_i}

        Example:
            >>> easel.VectorF([0.1, 0.1, 0.3, 0.5]).entropy()
            1.6854...
            >>> easel.VectorF([0.25, 0.25, 0.25, 0.25]).entropy()
            2.0

        References:
            - Cover, Thomas M., and Thomas, Joy A.
              *Entropy, Relative Entropy, and Mutual Information*. In
              Elements of Information Theory, 13–55. Wiley (2005): 2.
              :doi:`10.1002/047174882X.ch2` :isbn:`9780471241959`.

        .. versionadded:: 0.4.10

        """
        assert self._data != NULL
        with nogil:
            return libeasel.vec.esl_vec_FEntropy(<float*> self._data, self._n)

    cpdef float max(self) except *:
        assert self._data != NULL
        if self._n == 0:
            raise ValueError("max() called on an empty vector")
        with nogil:
            return libeasel.vec.esl_vec_FMax(<float*> self._data, self._n)

    cpdef float min(self) except *:
        assert self._data != NULL
        if self._n == 0:
            raise ValueError("argmin() called on an empty vector")
        with nogil:
            return libeasel.vec.esl_vec_FMin(<float*> self._data, self._n)

    cpdef object normalize(self):
        """Normalize a vector so that all elements sum to 1.

        Caution:
            If sum is zero, sets all elements to :math:`\\frac{1}{n}`,
            where :math:`n` is the size of the vector.

        """
        assert self._data != NULL
        with nogil:
            libeasel.vec.esl_vec_FNorm(<float*> self._data, self._n)
        return None

    cpdef float relative_entropy(self, VectorF other) except *:
        """Compute the relative entropy between two probability vectors.

        The Shannon relative entropy of two probability vectors :math:`p`
        and :math:`q`, also known as the Kullback-Leibler divergence, is
        defined as:

        .. math::

            D(p \\parallel q) = \\sum_i  p_i \\log_2 \\frac{p_i}{q_i}.

        with :math:`D(p \\parallel q) = \\infty` per definition if
        :math:`q_i = 0` and :math:`p_i > 0` for any :math:`i`.

        Example:
            >>> v1 = easel.VectorF([0.1, 0.1, 0.3, 0.5])
            >>> v2 = easel.VectorF([0.25, 0.25, 0.25, 0.25])
            >>> v1.relative_entropy(v2)
            0.3145...
            >>> v2.relative_entropy(v1)   # this is no symmetric relation
            0.3452...

        References:
            - Cover, Thomas M., and Thomas, Joy A.
              *Entropy, Relative Entropy, and Mutual Information*. In
              Elements of Information Theory, 13–55. Wiley (2005): 2.
              :doi:`10.1002/047174882X.ch2` :isbn:`9780471241959`.

        .. versionadded:: 0.4.10

        """
        assert self._data != NULL
        assert other._data != NULL
        if self._n != other._n:
            raise ValueError("cannot compute relative entropy of vectors of different sizes")
        with nogil:
            return libeasel.vec.esl_vec_FRelEntropy(
                <float*> self._data,
                <float*> other._data,
                self._n
            )

    cpdef object reverse(self):
        assert self._data != NULL
        with nogil:
            libeasel.vec.esl_vec_FReverse(<float*> self._data, <float*> self._data, self._n)
        return None

    cpdef float sum(self):
        """Returns the scalar sum of all elements in the vector.

        Float summations use `Kahan's algorithm <https://w.wiki/4Wa5>`_, in
        order to minimize roundoff error accumulation. Additionally, they
        are most accurate if the vector is sorted in increasing order, from
        small to large, so you may consider sorting the vector before
        summing it.

        References:
            - Kahan, W.
              *Pracniques: Further Remarks on Reducing Truncation Errors*.
              Communications of the ACM 8, no. 1 (1 January 1965): 40.
              :doi:`10.1145/363707.363723`.

        """
        assert self._data != NULL
        with nogil:
            return libeasel.vec.esl_vec_FSum(<float*> self._data, self._n)


cdef class VectorU8(Vector):
    """A vector storing byte-sized unsigned integers.

    .. versionadded:: 0.4.0

    """

    # --- Magic methods ------------------------------------------------------

    def __init__(self, object iterable = ()):
        """__init__(self, iterable=())\n--\n

        Create a new byte vector from the given data.

        """
        cdef int          n
        cdef size_t       i
        cdef uint8_t      item
        cdef uint8_t[::1] view
        cdef int          n_alloc
        cdef uint8_t*     data

        # collect iterable if it's not a `Sized` object
        if not isinstance(iterable, collections.abc.Sized):
            iterable = array.array("B", iterable)
        n = len(iterable)

        # make sure __init__ is only called once
        if self._data != NULL:
            raise RuntimeError("Vector.__init__ must not be called more than once")
        # make sure the vector has a positive size
        if n < 0:
            raise ValueError("Cannot create a vector with negative size")

        # allocate vector storage
        self._allocate(n)

        # try to copy the memory quickly if *iterable* implements the buffer
        # protocol, or fall back to copying each element with a for loop
        data = <uint8_t*> self._data
        try:
            view = iterable
        except (TypeError, ValueError):
            for i, item in enumerate(iterable):
                data[i] = item
        else:
            with nogil:
                memcpy(data, &view[0], n * sizeof(uint8_t))

    def __eq__(self, object other):
        assert self._data != NULL

        cdef unsigned char[::1] buffer
        cdef int                cmp    = 0
        cdef const uint8_t*     data   = <const uint8_t*> self._data

        try:
            buffer = other
        except ValueError:
            return NotImplemented
        if buffer.ndim != 1:
            return NotImplemented
        if buffer.shape[0] != self._n:
            return False
        with nogil:
            cmp = 0 if self._n == 0 else memcmp(&buffer[0], data, self._n)
        return cmp == 0

    def __getitem__(self, object index):
        assert self._data != NULL

        cdef VectorU8  new
        cdef int       idx
        cdef ssize_t   start
        cdef ssize_t   stop
        cdef ssize_t   step
        cdef uint8_t*  data  = <uint8_t*> self._data

        if isinstance(index, slice):
            start, stop, step = index.indices(self._n)
            if step != 1:
                raise ValueError(f"cannot slice a Vector with step other than 1")
            new = VectorU8.__new__(VectorU8)
            new._owner = self
            new._n = new._shape[0] = stop - start
            new._data = NULL if new._n == 0 else <void*> &(data[start])
            return new
        else:
            idx = index
            if idx < 0:
                idx += self._n
            if idx < 0 or idx >= self._n:
                raise IndexError("vector index out of range")
            return data[idx]

    def __setitem__(self, object index, uint8_t value):
        assert self._data != NULL

        cdef ssize_t  x
        cdef uint8_t* data = <uint8_t*> self._data

        if isinstance(index, slice):
            for x in range(*index.indices(self._n)):
                data[x] = value
        else:
            x = index
            if x < 0:
                x += self._n
            if x < 0  or x >= self._n:
                raise IndexError("vector index out of range")
            data[x] = value

    def __iadd__(self, object other):
        assert self._data != NULL

        cdef int            i
        cdef VectorU8       other_vec
        cdef const uint8_t* other_data
        cdef uint8_t        other_u
        cdef uint8_t*       data       = <uint8_t*> self._data

        if isinstance(other, VectorU8):
            other_vec = other
            other_data = <uint8_t*> other_vec._data
            assert other_vec._data != NULL or other_vec._n == 0
            if self._n != other_vec._n:
                raise ValueError("cannot pairwise add vectors of different sizes")
            with nogil:
                for i in range(self._n):
                    data[i] += other_data[i]
        else:
            other_u = other
            with nogil:
                for i in range(self._n):
                    data[i] += other_u
        return self

    def __isub__(self, object other):
        assert self._data != NULL

        cdef int            i
        cdef VectorU8       other_vec
        cdef uint8_t        other_u
        cdef const uint8_t* other_data
        cdef uint8_t*       data       = <uint8_t*> self._data

        if isinstance(other, VectorU8):
            other_vec = other
            other_data = <uint8_t*> other_vec._data
            assert other_vec._data != NULL or other_vec._n == 0
            if self._n != other_vec._n:
                raise ValueError("cannot pairwise subtract vectors of different sizes")
            with nogil:
                for i in range(self._n):
                    data[i] -= other_data[i]
        else:
            other_u = other
            with nogil:
                for i in range(self._n):
                    data[i] -= other_u
        return self

    def __imul__(self, object other):
        assert self._data != NULL

        cdef int            i
        cdef VectorU8       other_vec
        cdef uint8_t        other_u
        cdef const uint8_t* other_data
        cdef uint8_t*       data       = <uint8_t*> self._data

        if isinstance(other, VectorU8):
            other_vec = other
            other_data = <uint8_t*> other_vec._data
            assert other_vec._data != NULL or other_vec._n == 0
            if self._n != other_vec._n:
                raise ValueError("cannot pairwise multiply vectors of different sizes")
            with nogil:
                for i in range(self._n):
                    data[i] *= other_data[i]
        else:
            other_u = other
            with nogil:
                for i in range(self._n):
                    data[i] *= other_u
        return self

    def __floordiv__(VectorU8 self, object other):
        assert self._data != NULL
        cdef VectorU8 new = self.copy()
        return new.__ifloordiv__(other)

    def __ifloordiv__(self, object other):
        assert self._data != NULL

        cdef int            i
        cdef VectorU8       other_vec
        cdef uint8_t        other_u
        cdef const uint8_t* other_data
        cdef uint8_t*       data       = <uint8_t*> self._data

        if isinstance(other, VectorU8):
            other_vec = other
            other_data = <uint8_t*> other_vec._data
            assert other_vec._data != NULL or other_vec._n == 0
            if self._n != other_vec._n:
                raise ValueError("cannot pairwise multiply vectors of different sizes")
            with nogil:
                for i in range(self._n):
                    data[i] //= other_data[i]
        else:
            other_u = other
            with nogil:
                for i in range(self._n):
                    data[i] //= other_u
        return self

    def __matmul__(VectorU8 self, object other):
        assert self._data != NULL

        cdef int            i
        cdef VectorU8       other_vec
        cdef uint8_t        other_u
        cdef const uint8_t* other_data
        cdef long int       res        = 0
        cdef uint8_t*       data       = <uint8_t*> self._data

        if not isinstance(other, VectorU8):
            return NotImplemented

        other_vec = other
        other_data = <uint8_t*> other_vec._data
        assert other_vec._data != NULL or other_vec._n == 0
        if self._n != other_vec._n:
            raise ValueError("cannot get dot-product of vectors of different sizes")
        with nogil:
            for i in range(self._n):
                res += data[i] * other_data[i]
        return res

    # --- Properties ---------------------------------------------------------

    @property
    def itemsize(self):
        return sizeof(uint8_t)

    @property
    def format(self):
        return "B"

    # --- Utility ------------------------------------------------------------

    cdef const char* _format(self) noexcept:
        return b"B"

    # --- Methods ------------------------------------------------------------

    cpdef int argmax(self) except -1:
        assert self._data != NULL

        cdef int            i
        cdef int            best = 0
        cdef const uint8_t* data = <uint8_t*> self._data

        if self._n == 0:
            raise ValueError("argmax() called on an empty vector")
        with nogil:
            for i in range(1, self._n):
                if data[i] > data[best]:
                    best = i
        return best

    cpdef int argmin(self) except -1:
        assert self._data != NULL

        cdef int            i
        cdef int            best = 0
        cdef const uint8_t* data = <uint8_t*> self._data

        if self._n == 0:
            raise ValueError("argmin() called on an empty vector")
        with nogil:
            for i in range(1, self._n):
                if data[i] < data[best]:
                    best = i
        return best

    cpdef VectorU8 copy(self):
        assert self._data != NULL

        cdef VectorU8 new
        cdef int      n_alloc

        new = VectorU8.__new__(VectorU8)
        new._n = new._shape[0] = self._n
        n_alloc = 1 if self._n == 0 else self._n

        new._data = calloc(n_alloc, sizeof(uint8_t))
        if new._data == NULL:
            raise AllocationError("uint8_t", sizeof(uint8_t), n_alloc)
        with nogil:
            memcpy(new._data, self._data, self._n * sizeof(uint8_t))

        return new

    cpdef uint8_t max(self) except *:
        assert self._data != NULL

        cdef int            i
        cdef uint8_t        best
        cdef const uint8_t* data = <uint8_t*> self._data

        if self._n == 0:
            raise ValueError("max() called on an empty vector")
        with nogil:
            best = data[0]
            for i in range(1, self._n):
                if data[i] > best:
                    best = data[i]
        return best

    cpdef uint8_t min(self) except *:
        assert self._data != NULL

        cdef int            i
        cdef uint8_t        best
        cdef const uint8_t* data = <uint8_t*> self._data

        if self._n == 0:
            raise ValueError("min() called on an empty vector")
        with nogil:
            best = data[0]
            for i in range(1, self._n):
                if data[i] < best:
                    best = data[i]
        return best

    cpdef object reverse(self):
        assert self._data != NULL

        cdef int      i
        cdef uint8_t  x
        cdef uint8_t* data = <uint8_t*> self._data

        with nogil:
            for i in range(self._n // 2):
                x = data[self._n - i - 1]
                data[self._n - i - 1] = data[i]
                data[i] = x

        return None

    cpdef uint8_t sum(self):
        """Returns the scalar sum of all elements in the vector.

        Caution:
            The sum is wrapping::

                >>> vec = VectorU8([255, 2])
                >>> vec.sum()
                1

        """
        assert self._data != NULL

        cdef int            i
        cdef uint8_t        sum  = 0
        cdef const uint8_t* data = <uint8_t*> self._data

        with nogil:
            for i in range(self._n):
                sum += data[i]
        return sum


@cython.no_gc_clear
cdef class Matrix:
    """An abstract 2D array of fixed size.

    .. versionadded:: 0.4.0

    """

    # --- Class methods ------------------------------------------------------

    @classmethod
    def zeros(cls, int m, int n):
        """Create a new :math:`m \\times n` matrix filled with zeros.
        """
        cdef Matrix mat = cls()
        if m < 0 or n < 0:
            raise ValueError("Cannot create a matrix with negative dimension")
        mat._allocate(m, n)
        return mat

    @classmethod
    def _from_raw_bytes(cls, buffer, int m, int n, str byteorder):
        """Create a new matrix using the given bytes to fill its contents.
        """
        cdef const uint8_t[::1] bytes
        cdef Matrix             mat      = cls.zeros(m, n)
        cdef size_t             itemsize = mat.itemsize
        cdef object             view     = memoryview(buffer)

        # fix endianness if needed
        if byteorder != SYS_BYTEORDER and mat.itemsize > 1:
            newbuffer = array.array(mat.format)
            newbuffer.frombytes(view)
            newbuffer.byteswap()
            view = memoryview(newbuffer)

        # assign the items
        bytes = view.cast("B")
        assert bytes.shape[0] == m * n * itemsize
        if n > 0 and m > 0:
            with nogil:
                memcpy(mat._data[0], &bytes[0], m * n * itemsize)
        return mat

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._owner = None
        self._n = self._m = 0
        self._shape = (self._n, self._m)
        self._data = NULL

    def __dealloc__(self):
        if self._owner is None and self._data != NULL:
            if self._m > 0:
                free(self._data[0])
            free(self._data)
        self._data = NULL

    def __init__(self, object iterable = ()):
        raise TypeError("Can't instantiate abstract class 'Matrix'")

    def __bool__(self):
        return self._m != 0 and self._n != 0

    def __len__(self):
        return self._m

    def __copy__(self):
        return self.copy()

    def __getbuffer__(self, Py_buffer* buffer, int flags):
        assert self._data != NULL

        if flags & PyBUF_FORMAT:
            buffer.format = <char*> self._format()
        else:
            buffer.format = NULL
        buffer.buf = self._data[0]
        buffer.internal = NULL
        buffer.itemsize = self.itemsize
        buffer.len = self._m * self._n * self.itemsize
        buffer.ndim = 2
        buffer.obj = self
        buffer.readonly = 0
        buffer.shape = self._shape
        buffer.suboffsets = NULL

        if SYS_IMPLEMENTATION_NAME == "pypy":
            buffer.internal = PyMem_Malloc(2*sizeof(Py_ssize_t))
            if buffer.internal == NULL:
                raise AllocationError("Py_ssize_t", sizeof(Py_ssize_t), 2)
            buffer.strides = <Py_ssize_t*> buffer.internal
            buffer.strides[0] = self.strides[0]
            buffer.strides[1] = self.strides[1]
        else:
            buffer.strides = NULL

    def __releasebuffer__(self, Py_buffer* buffer):
        if SYS_IMPLEMENTATION_NAME == "pypy":
            PyMem_Free(buffer.internal)
            buffer.internal = NULL

    def __repr__(self):
        cdef Vector row
        cdef str    ty  = type(self).__name__
        return f"{ty}({[list(row) for row in self]!r})"

    def __sizeof__(self):
        return (
            self._n * sizeof(void*)
          + self._m * self._n * self.itemsize
          + sizeof(self)
        )

    def __reduce_ex__(self, int protocol):
        assert self._data != NULL

        cdef object buffer

        # use out-of-band pickling (only supported since protocol 5, see
        # https://docs.python.org/3/library/pickle.html#out-of-band-buffers)
        if protocol >= 5:
            buffer = pickle.PickleBuffer(self)
        else:
            buffer = array.array(self.format)
            buffer.frombytes(memoryview(self).cast("b"))

        return self._from_raw_bytes, (buffer, self._m, self._n, SYS_BYTEORDER)

    def __add__(Matrix self, object other):
        assert self._data != NULL
        cdef Matrix new = self.copy()
        return new.__iadd__(other)

    def __mul__(Matrix self, object other):
        assert self._data != NULL
        cdef Matrix new = self.copy()
        return new.__imul__(other)

    # --- Properties ---------------------------------------------------------

    @property
    def shape(self):
        """`tuple`: The shape of the matrix.

        Example:
            >>> m = MatrixF([ [1.0, 2.0], [3.0, 4.0], [5.0, 6.0] ])
            >>> m.shape
            (3, 2)

        """
        return tuple(self._shape)

    @property
    def strides(self):
        """`tuple`: The strides of the matrix.
        """
        return (self._shape[1] * self.itemsize, self.itemsize)

    @property
    def itemsize(self):
        """`int`: The size of each item in the matrix, in bytes.

        .. versionadded:: 0.4.7

        """

    @property
    def format(self):
        """`str`: The format of each item in the matrix.

        See Also:
            The `array` module of the Python standard library for a detail
            about available type codes.

        .. versionadded:: 0.4.7

        """

    # --- Utility ------------------------------------------------------------

    cdef const char* _format(self) noexcept:
        return NULL

    cdef int _allocate(self, size_t m, size_t n) except -1:
        # NB(@althonos): malloc and calloc are not guaranteed to return a
        #                pointer when called with a null allocation size,
        #                so we allocate a single item instead
        cdef int i
        cdef int m_alloc  = 1 if m == 0 else m
        cdef int mn_alloc = 1 if m == 0 or n == 0 else m * n
        cdef int itemsize = self.itemsize

        if self._data != NULL:
            free(self._data)

        self._m = self._shape[0] = m
        self._n = self._shape[1] = n

        # allocate the pointer array
        with nogil:
            self._data = <void**> calloc(m_alloc, sizeof(void*))
        if self._data == NULL:
            raise AllocationError("void*", sizeof(void*), m_alloc)

        # allocate the data array
        with nogil:
            self._data[0] = <void*> calloc(mn_alloc, itemsize)
        if self._data[0] == NULL:
            raise AllocationError("uint8_t", 1, itemsize * mn_alloc)

        # update the pointer offsets in the array of pointers
        for i in range(1, self._m):
            self._data[i] = self._data[0] + i * self._n * itemsize

        return 0


    # --- Methods ------------------------------------------------------------

    def argmax(self):
        """Return the coordinates of the maximum element in the matrix.

        Raises:
            `ValueError`: When called on an empty matrix.

        """

    def argmin(self):
        """Return the coordinates of the minimum element in the matrix.

        Raises:
            `ValueError`: When called on an empty matrix.

        """

    def copy(self):
        """Create a copy of the matrix, allocating a new buffer.
        """

    def max(self):
        """Return the value of the maximum element in the matrix.

        Raises:
            `ValueError`: When called on an empty matrix.

        """

    def min(self):
        """Return the value of the minimum element in the matrix.

        Raises:
            `ValueError`: When called on an empty matrix.

        """

    def sum(self):
        """Return the sum of all elements in the matrix.
        """


cdef class MatrixF(Matrix):
    """A matrix storing single-precision floating point numbers.

    Use indexing notation to access and edit individual elements of the
    matrix::

        >>> m = MatrixF.zeros(2, 2)
        >>> m[0, 0] = 3.0
        >>> m
        MatrixF([[3.0, 0.0], [0.0, 0.0]])

    Indexing can also be performed at the row-level to get a `VectorF`
    without copying the underlying data::

        >>> m = MatrixF([ [1.0, 2.0], [3.0, 4.0] ])
        >>> m[0]
        VectorF([1.0, 2.0])

    Objects of this type support the buffer protocol, and can be viewed
    as a `numpy.ndarray` with two dimensions using the `numpy.asarray`
    function, and can be passed without copy to most `numpy` functions::

        >>> m = MatrixF([ [1.0, 2.0], [3.0, 4.0] ])
        >>> numpy.asarray(m)
        array([[1., 2.],
               [3., 4.]], dtype=float32)
        >>> numpy.log2(m)
        array([[0.       , 1.       ],
               [1.5849625, 2.       ]], dtype=float32)

    .. versionadded:: 0.4.0

    """

    # --- Magic methods ------------------------------------------------------

    def __init__(self, object iterable = ()):
        """__init__(self, iterable=())\n--\n

        Create a new matrix from an iterable of rows.

        """
        cdef int     i
        cdef int     j
        cdef size_t  m
        cdef size_t  n
        cdef object  row
        cdef float   val
        cdef object  peeking
        cdef float** data    = NULL

        # collect iterable if it's not a `Sized` object
        if not isinstance(iterable, collections.abc.Sized):
            iterable = [array.array("f", row) for row in iterable]

        # make sure __init__ is only called once
        if self._data != NULL:
            raise RuntimeError("Matrix.__init__ must not be called more than once")

        # allow peeking the data
        peeking = peekable(iterable)
        # get the number of columns from the iterable
        m = len(iterable)
        # get the number of rows from the first element of the iterable
        n = 0 if m == 0 else len(peeking.peek())

        # allocate the buffer
        self._allocate(m, n)

        # assign the items
        data = <float**> self._data
        for i, row in enumerate(peeking):
            if len(row) != self._n:
                raise ValueError("Inconsistent number of rows in input")
            for j, val in enumerate(row):
                data[i][j] = val

    def __eq__(self, object other):
        assert self._data != NULL
        cdef MatrixF other_
        cdef int     i
        cdef int     j
        # check matrix type
        if not isinstance(other, MatrixF):
            return NotImplemented
        # check dimensions
        other_ = other
        assert other_._data != NULL
        if self._m != other_._m or self._n != other_._n:
            return False
        # check values
        for i in range(self._m):
            for j in range(self._n):
                if (<float**> self._data)[i][j] != (<float**> other_._data)[i][j]:
                    return False
        return True

    def __iadd__(self, object other):
        assert self._data != NULL

        cdef MatrixF other_mat
        cdef float   other_f

        if isinstance(other, MatrixF):
            other_mat = other
            assert other_mat._data != NULL
            if other_mat._m != self._m or other_mat._n != self._n:
                raise ValueError(f"cannot pairwise add {other_mat.shape} matrix to {self.shape} matrix")
            with nogil:
                libeasel.vec.esl_vec_FAdd(<float*> self._data[0], <float*> other_mat._data[0], self._m * self._n)
        else:
            other_f = other
            with nogil:
                libeasel.vec.esl_vec_FIncrement(<float*> self._data[0], self._m*self._n, other_f)
        return self

    def __imul__(self, object other):
        assert self._data != NULL

        cdef MatrixF other_mat
        cdef float   other_f
        cdef int     i

        if isinstance(other, MatrixF):
            other_mat = other
            assert other_mat._data != NULL
            if other_mat._m != self._m or other_mat._n != self._n:
                raise ValueError(f"cannot pairwise multiply {other_mat.shape} matrix with {self.shape} matrix")
            # NB(@althonos): There is no function in `vectorops.h` to do this
            # for now...
            for i in range(self._n * self._m):
                (<float**> self._data)[0][i] *= (<float**> other_mat._data)[0][i]
        else:
            other_f = other
            with nogil:
                libeasel.matrixops.esl_mat_FScale(<float**> self._data, self._m, self._n, other_f)
        return self

    def __getitem__(self, object index):
        assert self._data != NULL

        cdef int     x
        cdef int     y
        cdef str     ty
        cdef VectorF row
        cdef MatrixF new
        cdef float** data = <float**> self._data

        if isinstance(index, int):
            x = index
            if x < 0:
                x += self._m
            if x < 0 or x >= self._m:
                raise IndexError("vector index out of range")

            row = VectorF.__new__(VectorF)
            row._owner = self
            row._n = row._shape[0] = self._n
            row._data = <void*> &(data[x][0])
            return row

        elif isinstance(index, slice):
            start, stop, step = index.indices(self._m)
            if stop < 0 or stop > self._m or start < 0 or start >= self._m:
                raise IndexError("matrix row index out of range")

            new = MatrixF.__new__(MatrixF)
            new._owner = self
            new._m = new._shape[0] = stop - start
            new._n = new._shape[1] = self._n
            new._data = <void**> &data[start]
            return new

        elif isinstance(index, tuple):
            x, y = index
            if x < 0:
                x += self._m
            if y < 0:
                y += self._n
            if x < 0 or x >= self._m:
                raise IndexError("matrix row index out of range")
            if y < 0 or y >= self._n:
                raise IndexError("matrix column index out of range")
            return data[x][y]

        else:
            ty = type(index).__name__
            raise TypeError(f"expected integer, tuple or slice, found {ty}")

    def __setitem__(self, object index, float value):
        assert self._data != NULL

        cdef int x
        cdef int y
        cdef float** data = <float**> self._data

        if isinstance(index, tuple):
            x, y = index
            if x < 0:
                x += self._m
            if y < 0:
                y += self._n
            if x < 0 or x >= self._m:
                raise IndexError("matrix row index out of range")
            if y < 0 or y >= self._n:
                raise IndexError("matrix column index out of range")
            data[x][y] = value

        else:
            raise TypeError("Matrix.__setitem__ can only be used with a 2D index")

    # --- Properties ---------------------------------------------------------

    @property
    def itemsize(self):
        return sizeof(float)

    @property
    def format(self):
        return "f"

    # --- Utility ------------------------------------------------------------

    cdef const char* _format(self) noexcept:
        return b"f"

    # --- Methods ------------------------------------------------------------

    cpdef tuple argmax(self):
        assert self._data != NULL

        cdef int n
        cdef int x
        cdef int y

        with nogil:
            n = libeasel.vec.esl_vec_FArgMax(<float*> self._data[0], self._m*self._n)
            x = n // self._m
            y = n % self._n

        return x, y

    cpdef tuple argmin(self):
        assert self._data != NULL

        cdef int n
        cdef int x
        cdef int y

        with nogil:
            n = libeasel.vec.esl_vec_FArgMin(<float*> self._data[0], self._m*self._n)
            x = n // self._m
            y = n % self._n

        return x, y

    cpdef MatrixF copy(self):
        cdef MatrixF mat = MatrixF.__new__(MatrixF)
        mat._m = mat._shape[0] = self._m
        mat._n = mat._shape[1] = self._n
        with nogil:
            mat._data = <void**> libeasel.matrixops.esl_mat_FClone(<float**> self._data, self._m, self._n)
        if mat._data == NULL:
            raise AllocationError("float**", sizeof(float), self._m * self._n)
        return mat

    cpdef float max(self):
        assert self._data != NULL
        with nogil:
            return libeasel.matrixops.esl_mat_FMax(<float**> self._data, self._m, self._n)

    cpdef float min(self):
        assert self._data != NULL
        with nogil:
            return libeasel.vec.esl_vec_FMin(<float*> self._data[0], self._m*self._n)

    cpdef float sum(self):
        assert self._data != NULL
        with nogil:
            return libeasel.vec.esl_vec_FSum(<float*> self._data[0], self._m*self._n)


cdef class MatrixU8(Matrix):
    """A matrix storing byte-sized unsigned integers.

    .. versionadded:: 0.4.0

    """

    # --- Magic methods ------------------------------------------------------

    def __init__(self, object iterable = ()):
        """__init__(self, iterable=())\n--\n

        Create a new matrix from an iterable of rows.

        """
        cdef int       i
        cdef int       j
        cdef size_t    m
        cdef size_t    n
        cdef object    row
        cdef uint8_t   val
        cdef object    peeking
        cdef uint8_t** data    = NULL

        # collect iterable if it's not a `Sized` object
        if not isinstance(iterable, collections.abc.Sized):
            iterable = [array.array("B", row) for row in iterable]

        # make sure __init__ is only called once
        if self._data != NULL:
            raise RuntimeError("Matrix.__init__ must not be called more than once")

        # allow peeking the data
        peeking = peekable(iterable)
        # get the number of columns from the iterable
        m = len(iterable)
        # get the number of rows from the first element of the iterable
        n = 0 if m == 0 else len(peeking.peek())

        # allocate the buffer
        self._allocate(m, n)

        # assign the items
        data = <uint8_t**> self._data
        for i, row in enumerate(peeking):
            if len(row) != self._n:
                raise ValueError("Inconsistent number of rows in input")
            for j, val in enumerate(row):
                data[i][j] = val

    def __eq__(self, object other):
        assert self._data != NULL

        cdef uint8_t* data       = <uint8_t*> self._data[0]
        cdef uint8_t* other_data
        cdef MatrixU8 other_mat
        cdef int      i

        # check matrix type
        if not isinstance(other, MatrixU8):
            return NotImplemented
        # check dimensions
        other_mat = other
        assert other_mat._data != NULL
        other_data = <uint8_t*> other_mat._data[0]
        if self._m != other_mat._m or self._n != other_mat._n:
            return False
        # check values
        for i in range(self._m * self._n):
            if data[i] != other_data[i]:
                return False
        return True

    def __iadd__(self, object other):
        assert self._data != NULL

        cdef uint8_t* data       = <uint8_t*> self._data[0]
        cdef uint8_t* other_data
        cdef MatrixU8 other_mat
        cdef uint8_t  other_n
        cdef int      i

        if isinstance(other, MatrixU8):
            other_mat = other
            assert other_mat._data != NULL
            other_data = <uint8_t*> other_mat._data[0]
            if other_mat._m != self._m or other_mat._n != self._n:
                raise ValueError(f"cannot pairwise add {other_mat.shape} matrix to {self.shape} matrix")
            with nogil:
                for i in range(self._m * self._n):
                    data[i] += other_data[i]
        else:
            other_n = other
            with nogil:
                for i in range(self._m * self._n):
                    data[i] += other_n
        return self

    def __imul__(self, object other):
        assert self._data != NULL

        cdef uint8_t* data       = <uint8_t*> self._data[0]
        cdef uint8_t* other_data
        cdef MatrixU8 other_mat
        cdef uint8_t  other_n
        cdef int      i

        if isinstance(other, MatrixU8):
            other_mat = other
            assert other_mat._data != NULL
            other_data = <uint8_t*> other_mat._data[0]
            if other_mat._m != self._m or other_mat._n != self._n:
                raise ValueError(f"cannot pairwise multiply {other_mat.shape} matrix with {self.shape} matrix")
            # NB(@althonos): There is no function in `vectorops.h` to do this
            # for now...
            for i in range(self._n * self._m):
                data[i] *= other_data[i]
        else:
            other_n = other
            with nogil:
                for i in range(self._n * self._m):
                    data[i] *= other_n
        return self

    def __getitem__(self, object index):
        assert self._data != NULL

        cdef int       x
        cdef int       y
        cdef str       ty
        cdef VectorU8  row
        cdef MatrixU8  new
        cdef uint8_t** data = <uint8_t**> self._data

        if isinstance(index, int):
            x = index
            if x < 0:
                x += self._m
            if x < 0 or x >= self._m:
                raise IndexError("vector index out of range")

            row = VectorU8.__new__(VectorU8)
            row._owner = self
            row._n = row._shape[0] = self._n
            row._data = data[x]
            return row

        elif isinstance(index, slice):
            start, stop, step = index.indices(self._m)
            if stop < 0 or stop > self._m or start < 0 or start >= self._m:
                raise IndexError("matrix row index out of range")

            new = MatrixU8.__new__(MatrixU8)
            new._owner = self
            new._m = new._shape[0] = stop - start
            new._n = new._shape[1] = self._n
            new._data = <void**> &data[start]
            return new

        elif isinstance(index, tuple):
            x, y = index
            if x < 0:
                x += self._m
            if y < 0:
                y += self._n
            if x < 0 or x >= self._m:
                raise IndexError("matrix row index out of range")
            if y < 0 or y >= self._n:
                raise IndexError("matrix column index out of range")
            return data[x][y]

        else:
            ty = type(index).__name__
            raise TypeError(f"expected integer, tuple or slice, found {ty}")

    def __setitem__(self, object index, uint8_t value):
        assert self._data != NULL

        cdef int       x
        cdef int       y
        cdef uint8_t** data = <uint8_t**> self._data

        if isinstance(index, tuple):
            x, y = index
            if x < 0:
                x += self._m
            if y < 0:
                y += self._n
            if x < 0 or x >= self._m:
                raise IndexError("matrix row index out of range")
            if y < 0 or y >= self._n:
                raise IndexError("matrix column index out of range")
            data[x][y] = value

        else:
            raise TypeError("Matrix.__setitem__ can only be used with a 2D index")

    # --- Properties ---------------------------------------------------------

    @property
    def itemsize(self):
        return sizeof(uint8_t)

    @property
    def format(self):
        return "B"

    # --- Utility ------------------------------------------------------------

    cdef const char* _format(self) noexcept:
        return b"B"

    # --- Methods ------------------------------------------------------------

    cpdef tuple argmax(self):
        assert self._data != NULL

        cdef       int      i
        cdef       int      x
        cdef       int      y
        cdef       int      best = 0
        cdef const uint8_t* data = <uint8_t*> self._data[0]

        with nogil:
            for i in range(1, self._m * self._n):
                if data[i] > data[best]:
                    best = i
            x = best // self._m
            y = best % self._n
        return x, y

    cpdef tuple argmin(self):
        assert self._data != NULL

        cdef       int      i
        cdef       int      x
        cdef       int      y
        cdef       int      best = 0
        cdef const uint8_t* data = <uint8_t*> self._data[0]

        with nogil:
            for i in range(1, self._m * self._n):
                if data[i] < data[best]:
                    best = i
            x = best // self._m
            y = best % self._n
        return x, y

    cpdef MatrixU8 copy(self):

        cdef int      i
        cdef MatrixU8 mat = MatrixU8.__new__(MatrixU8)
        mat._m = mat._shape[0] = self._m
        mat._n = mat._shape[1] = self._n

        with nogil:
            # allocate array of pointers
            mat._data = <void**> malloc(sizeof(uint8_t*) * self._m)
            if mat._data == NULL:
                raise AllocationError("uint8_t*", sizeof(uint8_t*), self._m)
            # allocate memory block
            mat._data[0] = <void*> malloc(sizeof(uint8_t) * self._m * self._n)
            if mat._data == NULL:
                raise AllocationError("uint8_t", sizeof(uint8_t), self._m * self._n)
            # update array of pointers
            for i in range(self._m):
                mat._data[i] = mat._data[0] + i * self._n
            # copy data
            memcpy(mat._data[0], self._data[0], self._m * self._n * sizeof(uint8_t))

        return mat

    cpdef uint8_t max(self):
        assert self._data != NULL

        cdef       int      i
        cdef const uint8_t* data = <uint8_t*> self._data[0]
        cdef       uint8_t  best = data[0]

        with nogil:
            for i in range(1, self._m * self._n):
                if data[i] > best:
                    best = data[i]
        return best

    cpdef uint8_t min(self):
        assert self._data != NULL

        cdef       int      i
        cdef const uint8_t* data = <uint8_t*> self._data[0]
        cdef       uint8_t  best  = data[0]

        with nogil:
            for i in range(1, self._m * self._n):
                if data[i] < best:
                    best = data[i]
        return best

    cpdef uint8_t sum(self):
        assert self._data != NULL

        cdef       int      i
        cdef const uint8_t* data = <uint8_t*> self._data[0]
        cdef       uint8_t  sum  = 0

        with nogil:
            for i in range(self._m * self._n):
                sum += data[i]
        return sum


# --- Multiple Sequences Alignment -------------------------------------------

class _MSASequences(collections.abc.Sequence):
    """A read-only view over the individual sequences of an MSA.
    """

    __slots__ = ("msa",)

    def __init__(self, MSA msa):
        self.msa = msa

    def __len__(self):
        assert (<MSA> self.msa)._msa != NULL
        return (<MSA> self.msa)._msa.nseq


class _MSAAlignment(collections.abc.Sequence):
    """A read-only view over the aligned sequences (i.e. rows) of an MSA.
    """
    __slots__ = ("msa",)

    def __init__(self, MSA msa):
        self.msa = msa

    def __len__(self):
        assert (<MSA> self.msa)._msa != NULL
        return (<MSA> self.msa)._msa.nseq


class _MSAIndex(collections.abc.Mapping):
    """A read-only mapping of sequence names to sequences of an MSA.
    """

    __slots__ = ("msa",)

    def __init__(self, MSA msa):
        assert msa._msa != NULL

        cdef int          status
        cdef int          nseq   = msa._msa.nseq
        cdef ESL_KEYHASH* kh     = msa._msa.index

        if libeasel.keyhash.esl_keyhash_GetNumber(kh) != nseq:
            with nogil:
                msa._rehash()
        self.msa = msa

    def __getitem__(self, object item):
        cdef int                      status
        cdef int                      index  = -1
        cdef const unsigned char[::1] key    = item
        cdef esl_pos_t                length = key.shape[0]
        cdef MSA                      msa    = self.msa
        cdef ESL_KEYHASH*             kh     = msa._msa.index

        with nogil:
            status = libeasel.keyhash.esl_keyhash_Lookup(kh, <const char*> &key[0], length, &index)
        if status == libeasel.eslOK:
            return self.msa.sequences[index]
        elif status == libeasel.eslENOTFOUND:
            raise KeyError(item)
        else:
            raise UnexpectedError(status, "esl_keyhash_Lookup")

    def __len__(self):
        cdef ESL_KEYHASH* kh
        cdef MSA          msa = self.msa

        assert msa._msa != NULL
        assert msa._msa.index != NULL

        return libeasel.keyhash.esl_keyhash_GetNumber(msa._msa.index)

    def __iter__(self):
        cdef int          i
        cdef ESL_KEYHASH* kh
        cdef MSA          msa = self.msa

        assert msa._msa != NULL
        assert msa._msa.index != NULL

        kh = msa._msa.index
        for i in range(libeasel.keyhash.esl_keyhash_GetNumber(kh)):
            yield <bytes> libeasel.keyhash.esl_keyhash_Get(kh, i)


@cython.freelist(8)
cdef class MSA:
    """An abstract alignment of multiple sequences.

    Hint:
        Use ``len(msa)`` to get the number of columns in the alignment,
        and ``len(msa.sequences)`` to get the number of sequences (i.e.
        the number of rows).

    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._msa = NULL

    def __dealloc__(self):
        libeasel.msa.esl_msa_Destroy(self._msa)

    def __init__(self):
        raise TypeError("Can't instantiate abstract class 'MSA'")

    def __copy__(self):
        return self.copy()

    def __eq__(self, object other):
        assert self._msa != NULL

        cdef int status
        cdef MSA other_msa

        if not isinstance(other, MSA):
            return NotImplemented

        other_msa = <MSA> other
        status = libeasel.msa.esl_msa_Compare(self._msa, other_msa._msa)

        if status == libeasel.eslOK:
            return True
        elif status == libeasel.eslFAIL:
            return False
        else:
            raise UnexpectedError(status, "esl_msa_Compare")

    def __len__(self):
        assert self._msa != NULL
        if self._msa.nseq == 0:
            return 0
        return self._msa.alen

    # --- Properties ---------------------------------------------------------

    @property
    def accession(self):
        """`bytes` or `None`: The accession of the alignment, if any.
        """
        assert self._msa != NULL
        if self._msa.acc == NULL:
            return None
        return <bytes> self._msa.acc

    @accession.setter
    def accession(self, bytes accession):
        assert self._msa != NULL

        cdef       int       status
        cdef       esl_pos_t length
        cdef const char*     acc

        if accession is None:
            length = -1
            acc = NULL
        else:
            length = len(accession)
            acc = accession

        with nogil:
            status = libeasel.msa.esl_msa_SetAccession(self._msa, acc, length)
        if status == libeasel.eslEMEM:
            raise AllocationError("char", sizeof(char), length)
        elif status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_msa_SetAccession")

    @property
    def author(self):
        """`bytes` or `None`: The author of the alignment, if any.
        """
        assert self._msa != NULL
        if self._msa.au == NULL:
            return None
        return <bytes> self._msa.au

    @author.setter
    def author(self, bytes author):
        assert self._msa != NULL

        cdef       int       status
        cdef       esl_pos_t length
        cdef const char*     au

        if author is None:
            length = -1
            au = NULL
        else:
            length = len(author)
            au = author

        with nogil:
            status = libeasel.msa.esl_msa_SetAuthor(self._msa, au, length)
        if status == libeasel.eslEMEM:
            raise AllocationError("char", sizeof(char), length)
        elif status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_msa_SetAuthor")

    @property
    def name(self):
        """`bytes` or `None`: The name of the alignment, if any.
        """
        assert self._msa != NULL
        if self._msa.name == NULL:
            return None
        return <bytes> self._msa.name

    @name.setter
    def name(self, bytes name):
        assert self._msa != NULL

        cdef       int       status
        cdef       esl_pos_t length
        cdef const char*     nm

        if name is None:
            length = -1
            nm = NULL
        else:
            length = len(name)
            nm = name

        with nogil:
            status = libeasel.msa.esl_msa_SetName(self._msa, nm, length)
        if status == libeasel.eslEMEM:
            raise AllocationError("char", sizeof(char), length)
        elif status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_msa_SetName")

    @property
    def description(self):
        """`bytes` or `None`: The description of the alignment, if any.
        """
        assert self._msa != NULL
        if self._msa.desc == NULL:
            return None
        return <bytes> self._msa.desc

    @description.setter
    def description(self, bytes description):
        assert self._msa != NULL

        cdef       int       status
        cdef       esl_pos_t length
        cdef const char*     desc

        if description is None:
            length = -1
            desc = NULL
        else:
            length = len(description)
            desc = description

        with nogil:
            status = libeasel.msa.esl_msa_SetDesc(self._msa, desc, length)
        if status == libeasel.eslEMEM:
            raise AllocationError("char", sizeof(char), length)
        elif status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_msa_SetDesc")

    @property
    def names(self):
        """`tuple` of `bytes`: The name of each sequence in the alignment.

        Every sequence in the alignment is required to have a name, so
        no member of the `tuple` will ever be `None`.

        Example:
            >>> s1 = TextSequence(name=b"seq1", sequence="ATGC")
            >>> s2 = TextSequence(name=b"seq2", sequence="ATGC")
            >>> msa = TextMSA(name=b"msa", sequences=[s1, s2])
            >>> msa.names
            (b'seq1', b'seq2')

        .. versionadded:: 0.4.8

        """
        assert self._msa != NULL
        assert self._msa.sqname != NULL

        cdef int64_t i
        cdef bytes   name
        cdef tuple   names

        if self._msa.alen == 0 or self._msa.nseq == 0:
            return ()

        names = PyTuple_New(self._msa.nseq)
        for i in range(self._msa.nseq):
            # sequences must have a name
            assert self._msa.sqname[i] != NULL
            name = PyBytes_FromString(self._msa.sqname[i])
            Py_INCREF(name)
            PyTuple_SET_ITEM(names, i, name)

        return names

    @property
    def reference(self):
        """`bytes` or `None`: The reference annotation (`#=GC RF`), if any.

        .. versionadded:: 0.11.1

        """
        assert self._msa != NULL

        if self._msa.rf == NULL:
            return None
        return PyBytes_FromStringAndSize(self._msa.rf, self._msa.alen)

    @reference.setter
    def reference(self, reference: bytes):
        assert self._msa != NULL
        if reference is None:
            self._set_annotation(&self._msa.rf, NULL)
        else:
            self._set_annotation(&self._msa.rf, <char*> reference)

    @property
    def model_mask(self):
        """`bytes` or `None`: The model mask (`#=GC MM`), if any.

        .. versionadded:: 0.11.1

        """
        assert self._msa != NULL

        if self._msa.mm == NULL:
            return None
        return PyBytes_FromStringAndSize(self._msa.mm, self._msa.alen)

    @model_mask.setter
    def model_mask(self, model_mask: bytes):
        assert self._msa != NULL
        if model_mask is None:
            self._set_annotation(&self._msa.mm, NULL)
        else:
            self._set_annotation(&self._msa.mm, <char*> model_mask)

    @property
    def secondary_structure(self):
        """`bytes` or `None`: The consensus secondary structure, if any.

        .. versionadded:: 0.11.1

        """
        assert self._msa != NULL

        if self._msa.ss_cons == NULL:
            return None
        return PyBytes_FromStringAndSize(self._msa.ss_cons, self._msa.alen)

    @secondary_structure.setter
    def secondary_structure(self, secondary_structure: bytes):
        assert self._msa != NULL
        if secondary_structure is None:
            self._set_annotation(&self._msa.ss_cons, NULL)
        else:
            self._set_annotation(&self._msa.ss_cons, <char*> secondary_structure)

    @property
    def surface_accessibility(self):
        """`bytes` or `None`: The consensus surface accessibility, if any.

        .. versionadded:: 0.11.1

        """
        assert self._msa != NULL

        if self._msa.sa_cons == NULL:
            return None
        return PyBytes_FromStringAndSize(self._msa.sa_cons, self._msa.alen)

    @surface_accessibility.setter
    def surface_accessibility(self, surface_accessibility: bytes):
        assert self._msa != NULL
        if surface_accessibility is None:
            self._set_annotation(&self._msa.sa_cons, NULL)
        else:
            self._set_annotation(&self._msa.sa_cons, <char*> surface_accessibility)

    @property
    def posterior_probabilities(self):
        """`bytes` or `None`: The consensus posterior probabilities, if any.

        .. versionadded:: 0.11.1

        """
        assert self._msa != NULL

        if self._msa.pp_cons == NULL:
            return None
        return PyBytes_FromStringAndSize(self._msa.pp_cons, self._msa.alen)

    @posterior_probabilities.setter
    def posterior_probabilities(self, posterior_probabilities: bytes):
        assert self._msa != NULL
        if posterior_probabilities is None:
            self._set_annotation(&self._msa.pp_cons, NULL)
        else:
            self._set_annotation(&self._msa.pp_cons, <char*> posterior_probabilities)


    # TODO: Implement `weights` property exposing the sequence weights as
    #       a `Vector` object, needs implementation of a new `VectorD` class
    #       given that MSA.wgt is an array of `double`

    @property
    def indexed(self):
        """`~collections.abc.Mapping`: A mapping of names to sequences.

        This property can be used to access the sequence of a multiple
        sequence alignment by name. An index is created the first time this
        property is accessed.

        Raises:
            `KeyError`: When attempting to create an index for an alignment
                containing duplicate sequence names.

        Example:
            >>> s1 = TextSequence(name=b"seq1", sequence="ATGC")
            >>> s2 = TextSequence(name=b"seq2", sequence="ATTA")
            >>> msa = TextMSA(name=b"msa", sequences=[s1, s2])
            >>> msa.indexed[b'seq1'].sequence
            'ATGC'
            >>> msa.indexed[b'seq3']
            Traceback (most recent call last):
            ...
            KeyError: b'seq3'

        .. versionadded:: 0.11.1

        """
        return _MSAIndex(self)

    # --- Utils --------------------------------------------------------------

    cdef int _set_annotation(self, char** field, char* value) except 1 nogil:
        cdef size_t alen = self._msa.alen
        cdef size_t vlen
        if value == NULL:
            if field[0] == NULL:
                free(field[0])
        else:
            vlen = strlen(value)
            if vlen != alen:
                raise ValueError(f"invalid length (expected {alen}, found {vlen})")
            if field[0] == NULL:
                field[0] = <char*> calloc(alen, sizeof(char))
                if field[0] == NULL:
                    raise AllocationError("char", sizeof(char), alen)
            memcpy(field[0], value, alen * sizeof(char))
        return 0

    cdef int _rehash(self) except 1 nogil:
        """Rehash the sequence names for faster lookup.

        Raises:
            `KeyError`: When the multiple sequence alignment contains
                duplicate sequence names.

        """
        cdef int status = libeasel.msa.esl_msa_Hash(self._msa)
        if status == libeasel.eslOK:
            return 0
        elif status == libeasel.eslEDUP:
            raise KeyError("duplicate sequence names")
        else:
            raise UnexpectedError(status, "esl_msa_Hash")

    # --- Methods ------------------------------------------------------------

    cpdef uint32_t checksum(self):
        """Calculate a 32-bit checksum for the multiple sequence alignment.
        """
        cdef uint32_t checksum = 0
        cdef int status
        with nogil:
            status = libeasel.msa.esl_msa_Checksum(self._msa, &checksum)
        if status == libeasel.eslOK:
            return checksum
        else:
            raise UnexpectedError(status, "esl_msa_Checksum")

    cpdef MSA select(self, sequences = None, columns = None):
        """Select and copy a subset of the multiple sequence alignment.

        Arguments:
            sequences (iterable of `int`, or `None`): The indices of sequences
                to retain in the alignment subset. If `None` given, retain
                all sequences.
            columns (iterable of `int`, or `None`): The indices of columns to
                retain in the alignment subset. If `None` given, retain all
                columns.

        Raises:
            `IndexError`: When given indices that are out of bounds for the
                sequences or columns.

        Example:
            >>> s1 = TextSequence(name=b"seq1", sequence="ATGC")
            >>> s2 = TextSequence(name=b"seq2", sequence="ATCC")
            >>> s3 = TextSequence(name=b"seq3", sequence="ATGA")
            >>> msa = TextMSA(name=b"msa", sequences=[s1, s2, s3])
            >>> msa.select(sequences=[0, 2]).names
            (b'seq1', b'seq3')
            >>> tuple(msa.select(columns=range(1,4)).alignment)
            ('TGC', 'TCC', 'TGA')

        .. versionadded:: 0.11.1

        """
        assert self._msa != NULL

        # NB: This function does a copy first just to set the Python object
        #     properly; it would be

        cdef size_t              i
        cdef char[eslERRBUFSIZE] errbuf
        cdef int                 status
        cdef int*                mask   = NULL
        cdef MSA                 msa    = self.copy()
        cdef size_t              alen   = self._msa.alen
        cdef size_t              nseq   = self._msa.nseq

        try:
            if sequences is not None:
                # allocate mask array
                mask = <int*> realloc(mask, sizeof(int) * nseq)
                if mask == NULL:
                    raise AllocationError("int", sizeof(int), nseq)
                memset(mask, 0, sizeof(int) * nseq)
                # build array from Python arguments
                for i in sequences:
                    if i >= nseq:
                        raise IndexError(i)
                    mask[i] = True
                # clear old memory
                libeasel.msa.esl_msa_Destroy(msa._msa)
                # subset sequences
                status = libeasel.msa.esl_msa_SequenceSubset(self._msa, mask, &msa._msa)
                if status != libeasel.eslOK:
                    _reraise_error()
                    raise UnexpectedError(status, "esl_msa_SequenceSubset")
            if columns is not None:
                # allocate mask array
                mask = <int*> realloc(mask, sizeof(int) * alen)
                if mask == NULL:
                    raise AllocationError("int", sizeof(int), alen)
                memset(mask, 0, sizeof(int) * alen)
                # build array from Python arguments
                for i in columns:
                    if i >= alen:
                        raise IndexError(i)
                    mask[i] = True
                # subset sequences
                status = libeasel.msa.esl_msa_ColumnSubset(msa._msa, <char*> &errbuf, mask)
                if status != libeasel.eslOK:
                    _reraise_error()
                    raise UnexpectedError(status, "esl_msa_SequenceSubset")
        finally:
            free(mask)

        return msa


    cpdef void write(self, object fh, str format) except *:
        """Write the multiple sequence alignement to a file handle.

        Arguments:
            fh (`io.IOBase`): A Python file handle, opened in binary mode.
            format (`str`): The name of the multiple sequence alignment
                file format to use.

        .. versionadded:: 0.3.0

        """
        assert self._msa != NULL

        cdef int    fmt
        cdef int    status
        cdef FILE*  file   = NULL

        if format not in MSA_FILE_FORMATS:
            raise InvalidParameter("format", format, choices=list(MSA_FILE_FORMATS))

        fmt = MSA_FILE_FORMATS[format]

        try:
            file = fopen_obj(fh, "w")
            status = libeasel.msafile.esl_msafile_Write(file, self._msa, fmt)
        finally:
            if file is not NULL:
                fclose(file)

        if status != libeasel.eslOK:
            _reraise_error()
            raise UnexpectedError(status, "esl_sqascii_WriteFasta")


class _TextMSASequences(_MSASequences):
    """A read-only view over the sequences of an MSA in text mode.
    """

    def __init__(self, TextMSA msa):
        super().__init__(msa)

    def __getitem__(self, int idx):
        cdef int          status
        cdef TextSequence seq
        cdef TextMSA      msa    = self.msa

        assert msa._msa != NULL

        if idx < 0:
            idx += msa._msa.nseq
        if idx >= msa._msa.nseq or idx < 0:
            raise IndexError("list index out of range")

        seq = TextSequence.__new__(TextSequence)
        status = libeasel.sq.esl_sq_FetchFromMSA(msa._msa, idx, &seq._sq)
        if status == libeasel.eslOK:
            # TODO(@althonos): This is needed at the moment because of a bug
            #                  `esl_sq_FetchFromMSA`, remove when patch for
            #                  EddyRivasLab/easel#80 is released.
            seq._sq.abc = NULL
            return seq
        else:
            raise UnexpectedError(status, "esl_sq_FetchFromMSA")

    def __setitem__(self, int idx, TextSequence seq):
        cdef int status
        cdef int hash_index
        cdef MSA msa        = self.msa

        assert msa._msa != NULL
        assert seq._sq != NULL

        if idx < 0:
            idx += msa._msa.nseq
        if idx >= msa._msa.nseq or idx < 0:
            raise IndexError("list index out of range")

        # make sure the sequence has a name
        if seq.name is None:
            raise ValueError("cannot set an alignment sequence with an empty name")

        # make sure the sequence has the right length
        if len(seq) != msa._msa.alen:
            raise ValueError("sequence does not have the expected length")

        # make sure inserting the sequence will not create a name duplicate
        status = libeasel.keyhash.esl_keyhash_Lookup(
            msa._msa.index,
            seq._sq.name,
            -1,
            &hash_index
        )
        if status == libeasel.eslOK and hash_index != idx:
            raise ValueError("cannot set a sequence with a duplicate name")

        # set the new sequence
        with nogil:
            (<TextMSA> msa)._set_sequence(idx, seq._sq)
            if hash_index != idx:
                msa._rehash()

class _TextMSAAlignment(_MSAAlignment):
    """A read-only view over the alignment of an MSA in text mode.

    .. versionadded:: 0.11.1

    """

    def __init__(self, TextMSA msa):
        super().__init__(msa)

    def __getitem__(self, int idx):
        cdef int          status
        cdef MSA          msa    = self.msa

        assert msa._msa != NULL

        if idx < 0:
            idx += msa._msa.nseq
        if idx >= msa._msa.nseq or idx < 0:
            raise IndexError("list index out of range")

        return PyUnicode_DecodeASCII(msa._msa.aseq[idx], msa._msa.alen, NULL)


cdef class TextMSA(MSA):
    """A multiple sequence alignement stored in text mode.
    """

    # --- Magic methods ------------------------------------------------------

    def __init__(
        self,
        *args,
        bytes name=None,
        bytes description=None,
        bytes accession=None,
        object sequences=None,
        bytes author=None,
    ):
        """__init__(self, name=None, description=None, accession=None, sequences=None, author=None)\n--\n

        Create a new text-mode alignment with the given ``sequences``.

        Keyword Arguments:
            name (`bytes`, optional): The name of the alignment, if any.
            description (`bytes`, optional): The description of the
                alignment, if any.
            accession (`bytes`, optional): The accession of the alignment,
                if any.
            sequences (collection of `TextSequence`): The sequences to store
                in the multiple sequence alignment. All sequences must have
                the same length. They also need to have distinct names.
            author (`bytes`, optional): The author of the alignment, often
                used to record the aligner it was created with.

        Raises:
            `ValueError`: When the alignment cannot be created from the
                given sequences.
            `TypeError`: When ``sequences`` is not an iterable of
                `TextSequence` objects.

        Example:
            >>> s1 = TextSequence(name=b"seq1", sequence="ATGC")
            >>> s2 = TextSequence(name=b"seq2", sequence="ATGC")
            >>> msa = TextMSA(name=b"msa", sequences=[s1, s2])
            >>> len(msa)
            4

        .. versionchanged:: 0.3.0
           Allow creating an alignment from an iterable of `TextSequence`.

        .. deprecated:: 0.11.1
            Passing positional arguments to constructor.

        """

        # TODO: Remove in 0.12.0 (deprecation)
        if len(args) > 0:
            warnings.warn(
                "TextMSA.__init__ will not accept positional arguments after v0.12.0",
                category=DeprecationWarning
            )
            if len(args) > 0:
                name = args[0]
            if len(args) > 1:
                description = args[1]
            if len(args) > 2:
                accession = args[2]
            if len(args) > 3:
                sequences = args[3]
            if len(args) > 4:
                author = args[4]

        cdef TextSequence seq
        cdef int          i
        cdef list         seqs  = [] if sequences is None else list(sequences)
        cdef set          names = {seq.name for seq in seqs}
        cdef int64_t      alen  = len(seqs[0]) if seqs else -1
        cdef int          nseq  = len(seqs) if seqs else 1

        if len(names) != len(seqs):
            raise ValueError("duplicate names in alignment sequences")
        elif not all(len(seq) == alen for seq in seqs):
            raise ValueError("all sequences must have the same length")

        with nogil:
            self._msa = libeasel.msa.esl_msa_Create(nseq, alen)
        if self._msa == NULL:
            raise AllocationError("ESL_MSA", sizeof(ESL_MSA))

        if name is not None:
            self.name = name
        if accession is not None:
            self.accession = accession
        if description is not None:
            self.description = description
        if author is not None:
            self.author = author
        for i, seq in enumerate(seqs):
            self._set_sequence(i, seq._sq)

    # --- Properties ---------------------------------------------------------

    @property
    def alignment(self):
        """`collections.abc.Sequence`: A view of the aligned sequence data.

        This property gives access to the aligned sequences, including gap
        characters, so that they can be displayed or processed column by
        column.

        Examples:
            Use `TextMSA.alignment` to display an alignment in text
            format::

                >>> for name, aligned in zip(luxc.names, luxc.alignment):
                ...     print(name, " ", aligned[:40], "...")
                b'Q9KV99.1'   LANQPLEAILGLINEARKSWSST------------PELDP ...
                b'Q2WLE3.1'   IYSYPSEAMIEIINEYSKILCSD------------RKFLS ...
                b'Q97GS8.1'   VHDIKTEETIDLLDRCAKLWLDDNYSKK--HIETLAQITN ...
                b'Q3WCI9.1'   LLNVPLKEIIDFLVETGERIRDPRNTFMQDCIDRMAGTHV ...
                b'P08639.1'   LNDLNINNIINFLYTTGQRWKSEEYSRRRAYIRSLITYLG ...
                ...

            Use the splat operator (*) in combination with the `zip`
            builtin to iterate over the columns of an alignment:

                >>> for idx, col in enumerate(zip(*luxc.alignment)):
                ...     print(idx+1, col)
                1 ('L', 'I', 'V', 'L', 'L', ...)
                2 ('A', 'Y', 'H', 'L', 'N', ...)
                ...

        .. versionadded:: 0.4.8

        .. versionchanged:: 0.11.1
           Change the return type to a lazy `collections.abc.Sequence`.

        """
        assert self._msa != NULL
        assert not (self._msa.flags & libeasel.msa.eslMSA_DIGITAL)
        return _TextMSAAlignment(self)

    @property
    def sequences(self):
        """`collections.abc.Sequence`: A view of the alignment sequences.

        This property lets you access the individual sequences in the
        multiple sequence alignment as `TextSequence` instances.

        Example:
            Query the number of sequences in the alignment with `len`, or
            access individual members via indexing notation::

                >>> s1 = TextSequence(name=b"seq1", sequence="ATGC")
                >>> s2 = TextSequence(name=b"seq2", sequence="ATGC")
                >>> msa = TextMSA(name=b"msa", sequences=[s1, s2])
                >>> len(msa.sequences)
                2
                >>> msa.sequences[0].name
                b'seq1'

        Caution:
            Sequences in the list are copies, so editing their attributes
            will have no effect on the alignment::

                >>> msa.sequences[0].name
                b'seq1'
                >>> msa.sequences[0].name = b"seq1bis"
                >>> msa.sequences[0].name
                b'seq1'

            Support for this feature may be added in a future version, but
            can be circumvented for now by forcingly setting the updated
            version of the object::

                >>> seq = msa.sequences[0]
                >>> seq.name = b"seq1bis"
                >>> msa.sequences[0] = seq
                >>> msa.sequences[0].name
                b'seq1bis'

        .. versionadded:: 0.3.0

        """
        assert self._msa != NULL
        assert not (self._msa.flags & libeasel.msa.eslMSA_DIGITAL)
        return _TextMSASequences(self)

    # --- Utils --------------------------------------------------------------

    cdef int _set_sequence(self, int idx, ESL_SQ* seq) except 1 nogil:
        # assert seq.seq != NULL

        cdef int status

        status = libeasel.msa.esl_msa_SetSeqName(self._msa, idx, seq.name, -1)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_msa_SetSeqName")

        if seq.acc[0] != b'\0':
            status = libeasel.msa.esl_msa_SetSeqAccession(self._msa, idx, seq.acc, -1)
            if status != libeasel.eslOK:
                raise UnexpectedError(status, "esl_msa_SetSeqAccession")

        if seq.desc[0] != b'\0':
            status = libeasel.msa.esl_msa_SetSeqDescription(self._msa, idx, seq.desc, -1)
            if status != libeasel.eslOK:
                raise UnexpectedError(status, "esl_msa_SetSeqDescription")

        # assert self._msa.aseq[idx] != NULL
        strncpy(self._msa.aseq[idx], seq.seq, self._msa.alen)
        return 0

    # --- Methods ------------------------------------------------------------

    cpdef TextMSA copy(self):
        """Duplicate the text sequence alignment, and return the copy.
        """
        assert self._msa != NULL
        assert not (self._msa.flags & libeasel.msa.eslMSA_DIGITAL)

        cdef int status
        cdef TextMSA new = TextMSA.__new__(TextMSA)
        with nogil:
            new._msa = libeasel.msa.esl_msa_Clone(self._msa)
        if new._msa == NULL:
            raise AllocationError("ESL_MSA", sizeof(ESL_MSA))
        return new

    cpdef DigitalMSA digitize(self, Alphabet alphabet):
        """Convert the text alignment to a digital alignment using ``alphabet``.

        Returns:
            `DigitalMSA`: An alignment in digital mode containing the same
            sequences digitized with ``alphabet``.

        Raises:
            `ValueError`: When the text sequence contains invalid characters
                that cannot be converted according to ``alphabet.symbols``.

        """
        assert self._msa != NULL
        assert alphabet._abc != NULL

        cdef int                 status
        cdef DigitalMSA          new
        cdef char[eslERRBUFSIZE] errbuf

        new = DigitalMSA.__new__(DigitalMSA, alphabet)
        with nogil:
            new._msa = libeasel.msa.esl_msa_Clone(self._msa)
            if new._msa == NULL:
                raise AllocationError("ESL_MSA", sizeof(ESL_MSA))
            status = libeasel.msa.esl_msa_Digitize(alphabet._abc, new._msa, <char*> &errbuf)

        if status == libeasel.eslOK:
            assert new._msa.flags & libeasel.msa.eslMSA_DIGITAL
            return new
        elif status == libeasel.eslEINVAL:
            err_msg = errbuf.decode("utf-8", "replace")
            raise ValueError(f"Cannot digitize MSA with {alphabet.type} alphabet: {err_msg}")
        else:
            raise UnexpectedError(status, "esl_msa_Digitize")


class _DigitalMSASequences(_MSASequences):
    """A read-only view over the sequences of an MSA in digital mode.
    """

    __slots__ = ("msa", "alphabet")

    def __init__(self, DigitalMSA msa):
        super().__init__(msa)
        self.alphabet = msa.alphabet

    def __getitem__(self, int idx):
        cdef int             status
        cdef DigitalSequence seq
        cdef MSA             msa    = self.msa

        assert msa._msa != NULL

        if idx < 0:
            idx += msa._msa.nseq
        if idx >= msa._msa.nseq or idx < 0:
            raise IndexError("list index out of range")

        seq = DigitalSequence.__new__(DigitalSequence, self.alphabet)
        with nogil:
            status = libeasel.sq.esl_sq_FetchFromMSA(msa._msa, idx, &seq._sq)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_sq_FetchFromMSA")

        return seq

    def __setitem__(self, int idx, DigitalSequence seq):
        cdef int status
        cdef int hash_index
        cdef MSA msa        = self.msa

        assert msa._msa != NULL
        assert seq._sq != NULL

        if idx < 0:
            idx += msa._msa.nseq
        if idx >= msa._msa.nseq or idx < 0:
            raise IndexError("list index out of range")

        # make sure the sequence has a name
        if seq.name is None:
            raise ValueError("cannot set an alignment sequence with an empty name")

        # make sure the sequence has the right length
        if len(seq) != msa._msa.alen:
            raise ValueError("sequence does not have the expected length")

        # make sure the sequence has the right alphabet
        if not (<Alphabet> msa.alphabet)._eq(seq.alphabet):
            raise AlphabetMismatch(msa.alphabet, seq.alphabet)

        # make sure inserting the sequence will not create a name duplicate
        status = libeasel.keyhash.esl_keyhash_Lookup(
            msa._msa.index,
            seq._sq.name,
            -1,
            &hash_index
        )
        if status == libeasel.eslOK and hash_index != idx:
            raise ValueError("cannot set a sequence with a duplicate name")

        # set the new sequence
        with nogil:
            (<DigitalMSA> msa)._set_sequence(idx, seq._sq)
            if hash_index != idx:
                msa._rehash()


class _DigitalMSAAlignment(_MSAAlignment):
    """A read-only view over the alignment of an MSA in digital mode.

    .. versionadded:: 0.11.1

    """

    def __init__(self, DigitalMSA msa):
        super().__init__(msa)

    def __getitem__(self, int idx):
        cdef int          status
        cdef VectorU8     row
        cdef MSA          msa    = self.msa

        assert msa._msa != NULL

        if idx < 0:
            idx += msa._msa.nseq
        if idx >= msa._msa.nseq or idx < 0:
            raise IndexError("list index out of range")

        row = VectorU8.__new__(VectorU8)
        row._n = row._shape[0] = msa._msa.alen
        row._data = &msa._msa.ax[idx][1]
        row._owner = self
        return row


cdef class DigitalMSA(MSA):
    """A multiple sequence alignment stored in digital mode.

    Attributes:
        alphabet (`Alphabet`): The biological alphabet used to encode this
            sequence alignment to digits.

    """

    @classmethod
    def sample(
        cls,
        Alphabet alphabet not None,
        int max_sequences,
        int max_length,
        RandomnessOrSeed randomness = None,
    ):
        """Sample a sequence of length at most ``L`` at random.

        Arguments:
            alphabet (`~pyhmmer.easel.Alphabet`): The alphabet of the
                multiple sequence alignment.
            max_sequences (`int`): The maximum number of sequences to
                sample for the alignment (the actual sequence number is
                sampled).
            max_length (`int`): The maximum length of the alignment to
                generate (the actual sequence length is sampled).
            randomness (`~pyhmmer.easel.Randomness`, `int` or `None`): The
                random number generator to use for sampling, or a seed to
                initialize a generator. If `None` or ``0`` given, create
                a new random number generator with a random seed.

        Returns:
            `~pyhmmer.easel.DigitalMSA`: A new digital multiple sequence
            alignment generated at random.

        Hint:
            This constructor is only useful for testing and should not be
            used to generate random sequences to e.g. compute a background
            distribution for a statistical method, since this function
            samples alphabet residues at random irrespective of prior
            frequences.

        .. versionadded:: 0.11.1

        """
        cdef int        status
        cdef Randomness rng
        cdef DigitalMSA msa    = DigitalMSA.__new__(DigitalMSA, alphabet)

        if RandomnessOrSeed is Randomness:
            rng = randomness
        else:
            rng = Randomness(randomness)

        status = libeasel.msa.esl_msa_Sample(
            rng._rng,
            alphabet._abc,
            max_sequences,
            max_length,
            &msa._msa,
        )
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_msa_Sample")

        return msa

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self, Alphabet alphabet, *args, **kwargs):
        self._msa = NULL
        self.alphabet = alphabet

    def __init__(
        self,
        Alphabet alphabet,
        *args,
        bytes name=None,
        bytes description=None,
        bytes accession=None,
        object sequences=None,
        bytes author=None,
    ):
        """__init__(self, alphabet, name=None, description=None, accession=None, sequences=None, author=None)\n--\n

        Create a new digital-mode alignment with the given ``sequences``.

        Arguments:
            alphabet (`Alphabet`): The alphabet of the alignmed sequences.

        Keyword Arguments:
            name (`bytes`, optional): The name of the alignment, if any.
            description (`bytes`, optional): The description of the
                alignment, if any.
            accession (`bytes`, optional): The accession of the alignment,
                if any.
            sequences (iterable of `DigitalSequence`): The sequences to
                store in the multiple sequence alignment. All sequences must
                have the same length and alphabet. They also need to have
                distinct names set.
            author (`bytes`, optional): The author of the alignment, often
                used to record the aligner it was created with.

        .. versionchanged:: 0.3.0
           Allow creating an alignment from an iterable of `DigitalSequence`.

        .. deprecated:: 0.11.1
            Passing positional arguments other than ``alphabet``.

        """
        # TODO: Remove in 0.12.0 (deprecation)
        if len(args) > 0:
            warnings.warn(
                "DigitalMSA.__init__ will not accept positional arguments besides `alphabet` after v0.12.0",
                category=DeprecationWarning
            )
            if len(args) > 0:
                name = args[0]
            if len(args) > 1:
                description = args[1]
            if len(args) > 2:
                accession = args[2]
            if len(args) > 3:
                sequences = args[3]
            if len(args) > 4:
                author = args[4]

        cdef DigitalSequence seq
        cdef list            seqs  = [] if sequences is None else list(sequences)
        cdef set             names = { seq.name for seq in seqs }
        cdef int64_t         alen  = len(seqs[0]) if seqs else -1
        cdef int             nseq  = len(seqs) if seqs else 1

        if len(names) != len(seqs):
            raise ValueError("duplicate names in alignment sequences")

        for seq in seqs:
            if not isinstance(seq, DigitalSequence):
                ty = type(seq).__name__
                raise TypeError(f"expected DigitalSequence, found {ty}")
            elif len(seq) != alen:
                raise ValueError("all sequences must have the same length")
            elif not alphabet._eq(seq.alphabet):
                raise AlphabetMismatch(alphabet, seq.alphabet)

        self.alphabet = alphabet
        with nogil:
            self._msa = libeasel.msa.esl_msa_CreateDigital(alphabet._abc, nseq, alen)
        if self._msa == NULL:
            raise AllocationError("ESL_MSA", sizeof(ESL_MSA))

        if name is not None:
            self.name = name
        if accession is not None:
            self.accession = accession
        if description is not None:
            self.description = description
        if author is not None:
            self.author = author
        for i, seq in enumerate(seqs):
            self._set_sequence(i, seq._sq)


    # --- Properties ---------------------------------------------------------

    @property
    def alignment(self):
        """`collections.abc.Sequence`: A view of the aligned sequence data.

        This property gives access to the aligned rows of the alignment
        in their encoded form.

        .. versionadded:: 0.11.1

        """
        assert self._msa != NULL
        assert (self._msa.flags & libeasel.msa.eslMSA_DIGITAL)
        return _DigitalMSAAlignment(self)

    @property
    def sequences(self):
        """`collections.abc.Sequence`: A view of the alignment sequences.

        This property lets you access the individual sequences in the
        multiple sequence alignment as `DigitalSequence` instances.

        See Also:
            The documentation for the `TextMSA.sequences` property, which
            contains some additional information.

        .. versionadded:: 0.3.0

        """
        return _DigitalMSASequences(self)

    # --- Utils --------------------------------------------------------------

    cdef int _set_sequence(self, int idx, ESL_SQ* seq) except 1 nogil:
        # assert seq.dsq != NULL

        cdef int status

        status = libeasel.msa.esl_msa_SetSeqName(self._msa, idx, seq.name, -1)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_msa_SetSeqName")

        if seq.acc[0] != b'\0':
            status = libeasel.msa.esl_msa_SetSeqAccession(self._msa, idx, seq.acc, -1)
            if status != libeasel.eslOK:
                raise UnexpectedError(status, "esl_msa_SetSeqAccession")

        if seq.desc[0] != b'\0':
            status = libeasel.msa.esl_msa_SetSeqDescription(self._msa, idx, seq.desc, -1)
            if status != libeasel.eslOK:
                raise UnexpectedError(status, "esl_msa_SetSeqDescription")

        # assert self._msa.ax[idx] != NULL
        memcpy(self._msa.ax[idx], seq.dsq, (self._msa.alen+2) * sizeof(libeasel.ESL_DSQ))
        return 0

    # --- Methods ------------------------------------------------------------

    cpdef DigitalMSA copy(self):
        """Duplicate the digital sequence alignment, and return the copy.
        """
        assert self._msa != NULL
        assert self._msa.flags & libeasel.msa.eslMSA_DIGITAL

        cdef int           status
        cdef DigitalMSA    new    = DigitalMSA.__new__(DigitalMSA, self.alphabet)

        with nogil:
            new._msa = libeasel.msa.esl_msa_Clone(self._msa)
        if new._msa == NULL:
            raise AllocationError("ESL_MSA", sizeof(ESL_MSA))
        return new

    cpdef TextMSA textize(self):
        """Convert the digital alignment to a text alignment.

        Returns:
            `TextMSA`: A copy of the alignment in text-mode.

        .. versionadded:: 0.3.0

        """
        assert self._msa != NULL
        assert self.alphabet._abc != NULL

        cdef int     status
        cdef TextMSA new

        new = TextMSA.__new__(TextMSA)
        with nogil:
            new._msa = libeasel.msa.esl_msa_Clone(self._msa)
            if new._msa == NULL:
                raise AllocationError("ESL_MSA", sizeof(ESL_MSA))
            status = libeasel.msa.esl_msa_Textize(new._msa)

        if status == libeasel.eslOK:
            assert not (new._msa.flags & libeasel.msa.eslMSA_DIGITAL)
            return new
        else:
            raise UnexpectedError(status, "esl_msa_Textize")

    cpdef DigitalMSA identity_filter(
        self,
        float max_identity=0.8,
        float fragment_threshold=libeasel.msaweight.eslMSAWEIGHT_FRAGTHRESH,
        float consensus_fraction=libeasel.msaweight.eslMSAWEIGHT_SYMFRAC,
        bint ignore_rf=libeasel.msaweight.eslMSAWEIGHT_IGNORE_RF,
        bint sample=libeasel.msaweight.eslMSAWEIGHT_ALLOW_SAMP,
        int sample_threshold=libeasel.msaweight.eslMSAWEIGHT_SAMPTHRESH,
        int sample_count=libeasel.msaweight.eslMSAWEIGHT_NSAMP,
        int max_fragments=libeasel.msaweight.eslMSAWEIGHT_MAXFRAG,
        uint64_t seed=libeasel.msaweight.eslMSAWEIGHT_RNGSEED,
        str preference="conscover",
    ):
        r"""Filter the alignment sequences by percent identity.

        Arguments:
            max_identity (`float`): The maximum fractional identity
                allowed between two alignment sequences. Sequences with
                fractional identity :math:`\geq` this number will be
                removed, with the remaining one selected according to
                the given ``preference``.

        Keyword Arguments:
            fragment_threshold (`float`): The threshold for determining
                which sequences of the alignment are fragments. An
                sequence spanning columns :math:`i` to :math:`j` of an
                alignment of width :math:`W` will be flagged as a fragment
                if :math:`\frac{j - i}{ W } < \text{fragment_threshold}`,
            consensus_fraction (`float`): The parameter for determining
                with columns of the alignment are consensus columns.
                A column containing :math:`n` symbols and :math:`m` gaps
                will be marked a consensus column if
                :math:`\frac{n}{n + m} \ge \text{consensus_fraction}`.
            ignore_rf (`bool`): Set to `True` to ignore the *RF* line
                of the alignment (if present) and to force building the
                consensus.
            sample (`bool`): Whether or not to enable consensus determination
                by subsampling for large alignments. Set to `False` to force
                using all sequences.
            sample_threshold (`int`): The minimum number of sequences the
                alignment must contain to use subsampling for consensus
                determination (when ``sample`` is `True`).
            sample_count (`int`): The number of sequences to use when
                determining consensus by random subsampling.
            max_fragments (`int`): The maximum number of allowed fragments
                in the sample used for determining consensus. If the sample
                contains more than ``max_fragments`` fragments, the
                consensus determination is done with all sequences instead.
            seed (`int`): The seed to use for initializing the random
                number generator (used when ``preference`` is ``random``
                or when ``sample`` is `True`). If ``0`` or `None` is given,
                an arbitrary seed will be chosen using the system clock.
            preference (`str`): The strategy to use for selecting the
                representative sequence in case of duplicates. Supported
                strategies are ``conscover`` (the default), which prefers
                the sequence with an alignment span that covers more
                consensus columns; ``origorder`` to use the first sequence
                in the original alignment order; and ``random`` to select
                the sequence at random.

        Returns:
            `~pyhmmer.easel.MSA`: The multiple sequence alignments with
            duplicate sequence removed. Unparsed Sotckholm markup is not
            propagated.

        """
        assert self._msa != NULL

        cdef int status
        cdef int filterpref
        cdef ESL_MSAWEIGHT_CFG cfg
        cdef DigitalMSA msa = DigitalMSA.__new__(DigitalMSA, self.alphabet)

        # validate arguments
        if fragment_threshold < 0 or fragment_threshold > 1:
            raise InvalidParameter("fragment_threshold", fragment_threshold, hint="real number between 0 and 1")
        if consensus_fraction < 0 or consensus_fraction > 1:
            raise InvalidParameter("consensus_fraction", consensus_fraction, hint="real number between 0 and 1")
        if sample_threshold < 0:
            raise InvalidParameter("sample_threshold", sample_threshold, hint="positive integer")
        if sample_count < 0:
            raise InvalidParameter("sample_count", sample_count, hint="positive integer")
        if max_fragments < 0:
            raise InvalidParameter("max_fragments", max_fragments, hint="positive integer")
        if preference not in MSA_WEIGHT_PREFERENCES:
            raise InvalidParameter("preference", preference, choices=list(MSA_WEIGHT_PREFERENCES))

        # prepare configuration
        cfg.fragthresh = fragment_threshold
        cfg.symfrac = consensus_fraction
        cfg.sampthresh = sample_threshold
        cfg.ignore_rf = ignore_rf
        cfg.allow_samp = sample
        cfg.sampthresh = sample_threshold
        cfg.nsamp = sample_count
        cfg.maxfrag = max_fragments
        cfg.seed = seed
        cfg.filterpref = MSA_WEIGHT_PREFERENCES[preference]

        # run filtering
        with nogil:
            status = libeasel.msaweight.esl_msaweight_IDFilter_adv(&cfg, self._msa, max_identity, &msa._msa)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_msaweight_IDFilter_adv")

        return msa



# --- MSA File ---------------------------------------------------------------

cdef class MSAFile:
    """A wrapper around a multiple-alignment file.

    This class supports reading sequences stored in different formats, such
    as Stockholm, A2M, PSI-BLAST or Clustal.

    Attributes:
        name (`str`, *optional*): The name of the MSA file, if it was
            created from a filename, or `None` if it wraps a file-like
            object.

    Hint:
        Some Clustal files created by alignment tools other than Clustal
        (such as MUSCLE or MAFFT, for instance), may not contain the header
        expected by Easel for the Clustal format. If you get an error while
        trying to parse these files, use the ``"clustallike"`` format
        instead of the ``"clustal"`` format when creating the `MSAFile`.

    """

    _FORMATS = dict(MSA_FILE_FORMATS) # copy dict to prevent editing

    # --- Constructor --------------------------------------------------------

    @staticmethod
    cdef ESL_MSAFILE* _open_fileobj(object fh, int fmt) except NULL:
        """Get an ``ESL_MSAFILE*`` to read MSA from the file-like object ``fh``.

        Adapted from the ``esl_msafile_Open`` function in ``esl_msafile.c``.

        Raises:
            `~pyhmmer.errors.AllocationError`: When either the ``ESL_SQFILE``
                or one of the string allocations fails.
            `NotImplementedError`: When calling `HMMFile._open_fileobj` with
                ``eslSQFILE_NCBI`` as the sequence format.

        """
        cdef int          status
        cdef ESL_BUFFER*  buffer  = NULL
        cdef ESL_MSAFILE* msaf    = NULL
        cdef FILE*        fp      = fopen_obj(fh, "r")
        cdef bytes        fh_repr = repr(fh).encode("ascii")

        try:
            # get an ESL_BUFFER from the file object
            status = libeasel.buffer.esl_buffer_OpenStream(fp, &buffer)
            if status != libeasel.eslOK:
                raise UnexpectedError(status, "esl_buffer_OpenStream")
        except:
            fclose(fp)
            raise

        try:
            # get an ESL_MSAFILE from the ESL_BUFFER
            status = libeasel.msafile.esl_msafile_OpenBuffer(NULL, buffer, fmt, NULL, &msaf)
            if status == libeasel.eslOK:
                return msaf
            elif status == libeasel.eslENOFORMAT:
                raise ValueError("Could not determine format of file: {!r}".format(fh))
            elif status != libeasel.eslOK:
                raise UnexpectedError(status, "esl_msafile_OpenBuffer")
        except:
            libeasel.msafile.esl_msafile_Close(msaf)
            fclose(fp)
            raise

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self.alphabet = None
        self.name = None
        self._msaf = NULL

    def __init__(
        self,
        object file,
        str format = None,
        *,
        bint digital = False,
        Alphabet alphabet = None,
    ):
        """__init__(self, file, format=None, *, digital=False, alphabet=False)\n--\n

        Create a new MSA file parser wrapping the given ``file``.

        Arguments:
            file (`str` or file-like object): Either the path to a file
                containing the sequences to read, or a file-like object
                opened in **binary mode**.
            format (`str`, optional): The format of the file, or `None` to
                autodetect. Supported values are: ``stockholm``, ``pfam``,
                ``a2m``, ``psiblast``, ``selex``, ``afa`` (aligned FASTA),
                ``clustal``, ``clustallike``, ``phylip``, ``phylips``.

        Keyword Arguments:
            digital (`bool`): Whether to read the sequences in text or digital
                mode. This will affect the type of `MSA` objects returned
                later by the `read` function.
            alphabet (`~pyhmmer.easel.Alphabet`): The alphabet to use to
                digitize the sequences while reading.  If `None` given, it
                will be guessed based on the contents of the first sequence.

        Raises:
            `ValueError`: When ``format`` is not a valid MSA format.

        .. versionchanged:: 0.4.8
           Support reading from a file-like object.

        .. versionchanged:: 0.5.0
           Added the ``digital`` and ``alphabet`` keyword arguments.

        """
        cdef int   fmt
        cdef int   status
        cdef bytes fspath = None

        fmt = libeasel.msafile.eslMSAFILE_UNKNOWN
        if format is not None:
            format_ = format.lower()
            if format_ not in MSA_FILE_FORMATS:
                raise ValueError("Invalid MSA format: {!r}".format(format))
            fmt = MSA_FILE_FORMATS[format_]

        # open from either a file-like object or a path
        try:
            fspath = os.fsencode(file)
            # NOTE(@althonos): manually check the file is not a folder,
            # otherwise Easel will run into a segmentation fault after
            # failing to "slurp" the file!
            if os.path.isdir(fspath):
                raise IsADirectoryError(errno.EISDIR, f"Is a directory: {file!r}")
        except TypeError:
            self._msaf = MSAFile._open_fileobj(file, fmt)
            status = libeasel.eslOK
        else:
            self.name = os.fsdecode(fspath)
            status = libeasel.msafile.esl_msafile_Open(NULL, fspath, NULL, fmt, NULL, &self._msaf)

        # store a reference to the argument
        self._file = file

        # configure the new instance
        try:
            # check opening the file was successful
            if status == libeasel.eslENOTFOUND:
                raise FileNotFoundError(errno.ENOENT, f"No such file or directory: {file!r}")
            elif status == libeasel.eslEMEM:
                raise AllocationError("ESL_MSAFILE", sizeof(ESL_MSAFILE))
            elif status == libeasel.eslENOFORMAT:
                if format is None:
                    raise ValueError("Could not determine format of file: {!r}".format(file))
                else:
                    raise EOFError("Sequence file is empty")
            elif status != libeasel.eslOK:
                raise UnexpectedError(status, "esl_msafile_Open")
            # set digital mode if requested
            if digital:
                self.alphabet = self.guess_alphabet() if alphabet is None else alphabet
                if self.alphabet is None:
                    raise ValueError("Could not determine alphabet of file: {!r}".format(file))
                status = libeasel.msafile.esl_msafile_SetDigital(self._msaf, self.alphabet._abc)
                if status != libeasel.eslOK:
                    raise UnexpectedError(status, "esl_msafile_SetDigital")
        except Exception as err:
            self.close()
            raise err

    def __dealloc__(self):
        if self._msaf != NULL:
            warnings.warn("unclosed MSAFile", ResourceWarning)
            self.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def __iter__(self):
        return self

    def __next__(self):
        cdef MSA msa = self.read()
        if msa is None:
            raise StopIteration()
        return msa

    def __repr__(self):
        cdef str  name = type(self).__name__
        cdef list args = [
            f"format={self.format!r}",
            f"digital={self.digital}",
            f"alphabet={self.alphabet!r}"
        ]
        if self.name is None:
            return f"<{name} file={self.file!r} {' '.join(args)}>"
        else:
            return f"{name}({self.name!r}, {', '.join(args)})"

    # --- Properties ---------------------------------------------------------

    @property
    def closed(self):
        """`bool`: Whether the `MSAFile` is closed or not.
        """
        return self._msaf == NULL

    @property
    def digital(self):
        """`bool`: Whether the `MSAFile` is in digital mode or not.
        """
        return self.alphabet is not None

    @property
    def format(self):
        """`str`: The format of the `MSAFile`.
        """
        if self._msaf == NULL:
            raise ValueError("I/O operation on closed file.")
        return MSA_FILE_FORMATS_INDEX[self._msaf.format]

    # --- Methods ------------------------------------------------------------

    cpdef void close(self):
        """Close the file and free the resources used by the parser.
        """
        libeasel.msafile.esl_msafile_Close(self._msaf)
        self._msaf = NULL

    cpdef Alphabet guess_alphabet(self):
        """Guess the alphabet of an open `MSAFile`.

        This method tries to guess the alphabet of a multiple-alignment file
        by inspecting the first entry in the file. It returns the alphabet,
        or `None` if the file alphabet cannot be reliably guessed.

        Raises:
            `EOFError`: if the file is empty.
            `OSError`: if a parse error occurred.
            `ValueError`: if this methods is called after the file was closed.

        Example:
            >>> with MSAFile("tests/data/msa/LuxC.sto") as mf:
            ...     mf.guess_alphabet()
            Alphabet.amino()

        """
        cdef int      ty
        cdef int      status
        cdef Alphabet alphabet
        cdef str      msg

        if self._msaf == NULL:
            raise ValueError("I/O operation on closed file.")

        status = libeasel.msafile.esl_msafile_GuessAlphabet(self._msaf, &ty)
        if status == libeasel.eslOK:
            alphabet = Alphabet.__new__(Alphabet)
            alphabet._init_default(ty)
            return alphabet
        elif status == libeasel.eslENOALPHABET or status == libeasel.eslEOD:
            return None
        elif status == libeasel.eslENODATA:
            raise EOFError("Sequence file appears to be empty.")
        elif status == libeasel.eslEFORMAT:
            msg = self._msaf.errmsg.decode("utf-8", "replace")
            raise ValueError("Could not parse file: {}".format(msg))
        else:
            raise UnexpectedError(status, "esl_msafile_GuessAlphabet")

    cpdef MSA read(self):
        """Read the next alignment from the file.

        Returns:
            `MSA`: The next alignment in the file, or `None` if all the
            alignments were read from the file already.

        Raises:
            `ValueError`: When attempting to read an alignment from a closed
                file, or when the file could not be parsed.

        """
        cdef MSA msa

        if self.alphabet is not None:
            msa = DigitalMSA.__new__(DigitalMSA, self.alphabet)
        else:
            msa = TextMSA.__new__(TextMSA)

        if self._msaf == NULL:
            raise ValueError("I/O operation on closed file.")
        else:
            status = libeasel.msafile.esl_msafile_Read(self._msaf, &msa._msa)

        if status == libeasel.eslOK:
            return msa
        elif status == libeasel.eslEOF:
            return None
        elif status == libeasel.eslEFORMAT:
            msg = self._msaf.errmsg.decode("utf-8", "replace")
            raise ValueError("Could not parse file: {}".format(msg))
        else:
            _reraise_error()
            raise UnexpectedError(status, "esl_msafile_Read")


# --- Randomness -------------------------------------------------------------

@cython.no_gc_clear
cdef class Randomness:
    """A portable, thread-safe random number generator.

    Methods with an implementation in Easel are named after the equivalent
    methods of `random.Random`.

    .. versionadded:: 0.4.2

    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._rng = NULL
        self._owner = None

    def __dealloc__(self):
        if self._rng != NULL and self._owner is None:
            libeasel.random.esl_randomness_Destroy(self._rng)
        self._rng = NULL

    def __init__(self, object seed=None, bint fast=False):
        """__init__(self, seed=None, fast=False)\n--\n

        Create a new random number generator with the given seed.

        Arguments:
            seed (`int`): The seed to initialize the generator with. If ``0``
                or `None` is given, an arbitrary seed will be chosen using
                the system clock.
            fast (`bool`): If `True`, use a linear congruential generator
                (LCG), which is low quality and should only be used for
                integration with legacy code. With `False`, use the
                Mersenne Twister MT19937 algorithm instead.

        """
        cdef uint32_t _seed = 0 if seed is None else seed
        # Support calling __init__ multiple times
        if self._rng == NULL:
            if fast:
                self._rng = libeasel.random.esl_randomness_CreateFast(_seed)
            else:
                self._rng = libeasel.random.esl_randomness_Create(_seed)
            if self._rng == NULL:
                raise AllocationError("ESL_RANDOMNESS", sizeof(ESL_RANDOMNESS))
        else:
            self.seed(_seed)

    def __copy__(self):
        return self.copy()

    def __reduce__(self):
        state = self.getstate()
        return Randomness, (state[0], state[1]), state

    def __repr__(self):
        assert self._rng != NULL
        cdef type ty   = type(self)
        cdef str  name = ty.__name__
        cdef str  mod  = ty.__module__
        return f"{name}({self._rng.seed!r}, fast={self.fast!r})"

    def __getstate__(self):
        return self.getstate()

    def __setstate__(self, state):
        self.setstate(state)

    def __sizeof__(self):
        assert self._rng != NULL
        return sizeof(ESL_RANDOMNESS) + sizeof(self)

    # --- Properties ---------------------------------------------------------

    @property
    def fast(self):
        """`bool`: `True` when the linear congruential generator is in use.
        """
        assert self._rng != NULL
        return self._rng.type == libeasel.random.esl_randomness_type.eslRND_FAST


    # --- Methods ------------------------------------------------------------

    def getstate(self):
        """Get a tuple containing the current state.
        """
        assert self._rng != NULL
        if self.fast:
            return ( True, self._rng.seed, self._rng.x )
        else:
            return ( False, self._rng.seed, self._rng.mti, [self._rng.mt[x] for x in range(624)] )

    def setstate(self, tuple state):
        """Restores the state of the random number generator.
        """
        assert self._rng != NULL

        self._rng.seed = state[1]
        if state[0]:
            self._rng.type = libeasel.random.esl_randomness_type.eslRND_FAST
            self._rng.x = state[2]
        else:
            self._rng.type = libeasel.random.esl_randomness_type.eslRND_MERSENNE
            self._rng.mti = state[2]
            for x in range(624):
                self._rng.mt[x] = state[3][x]

    cpdef void seed(self, object n=None) except *:
        """Reinitialize the random number generator with the given seed.

        Arguments:
            n (`int`, optional): The seed to use. If ``0`` or `None`, an
                arbitrary seed will be chosen using the current time.

        """
        assert self._rng != NULL

        cdef int      status
        cdef uint32_t seed   = n if n is not None else 0

        with nogil:
            status = libeasel.random.esl_randomness_Init(self._rng, seed)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_randomness_Init")

    cpdef Randomness copy(self):
        """Return a copy of the random number generator in the same exact state.
        """
        assert self._rng != NULL

        cdef Randomness new = Randomness.__new__(Randomness)
        new._rng = <ESL_RANDOMNESS*> malloc(sizeof(ESL_RANDOMNESS))
        if new._rng == NULL:
            raise AllocationError("ESL_RANDOMNESS", sizeof(ESL_RANDOMNESS))

        memcpy(new._rng, self._rng, sizeof(ESL_RANDOMNESS))
        return new

    cpdef double random(self):
        """Generate a uniform random deviate on :math:`\\left[ 0, 1 \\right)`.
        """
        assert self._rng != NULL
        return libeasel.random.esl_random(self._rng)

    cpdef double normalvariate(self, double mu, double sigma):
        """Generate a Gaussian-distributed sample.

        Arguments:
            mu (`float`): The mean of the Gaussian being sampled.
            sigma (`float`): The standard deviation of the Gaussian being
                sampled.

        """
        assert self._rng != NULL
        return libeasel.random.esl_rnd_Gaussian(self._rng, mu, sigma)


# --- Sequence ---------------------------------------------------------------

@cython.freelist(8)
cdef class Sequence:
    """An abstract biological sequence with some associated metadata.

    Easel provides two different mode to store a sequence: text, or digital.
    In the HMMER code, changing from one mode to another mode is done in
    place, which allows recycling memory. However, doing so can be confusing
    since there is no way to know statically the representation of a sequence.

    To avoid this, ``pyhmmer`` provides two subclasses of the `Sequence`
    abstract class to maintain the mode contract: `TextSequence` and
    `DigitalSequence`. Functions expecting sequences in digital format, like
    `pyhmmer.hmmer.hmmsearch`, can then use Python type system to make sure
    they receive sequences in the right mode. This allows type checkers
    such as ``mypy`` to detect potential contract breaches at compile-time.

    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._sq = NULL

    def __init__(self):
        raise TypeError("Can't instantiate abstract class 'Sequence'")

    def __dealloc__(self):
        libeasel.sq.esl_sq_Destroy(self._sq)

    def __eq__(self, object other):
        assert self._sq != NULL

        cdef int      status
        cdef Sequence other_sq

        if not isinstance(other, Sequence):
            return NotImplemented

        other_sq = <Sequence> other
        status = libeasel.sq.esl_sq_Compare(self._sq, other_sq._sq)

        if status == libeasel.eslOK:
            return True
        elif status == libeasel.eslFAIL:
            return False
        else:
            raise UnexpectedError(status, "esl_sq_Compare")

    def __len__(self):
        assert self._sq != NULL
        if self._sq.n == -1:
            return 0
        return <int> self._sq.n

    def __copy__(self):
        return self.copy()

    def __sizeof__(self):
        assert self._sq != NULL
        cdef ssize_t i
        cdef size_t  extra_markup_size = 0
        for i in range(self._sq.nxr):
            extra_markup_size += strlen(self._sq.xr_tag[i]) * sizeof(char)
            extra_markup_size += self._sq.salloc * sizeof(char)
        return (
                sizeof(self)
            +   sizeof(ESL_SQ)
            +   self._sq.nalloc * sizeof(char)
            +   self._sq.aalloc * sizeof(char)
            +   self._sq.dalloc * sizeof(char)
            +   self._sq.salloc * sizeof(char)
            +   self._sq.srcalloc * sizeof(char)
            +   extra_markup_size
        )

    # --- Properties ---------------------------------------------------------

    @property
    def accession(self):
        """`bytes`: The accession of the sequence.
        """
        assert self._sq != NULL
        return <bytes> self._sq.acc

    @accession.setter
    def accession(self, bytes accession):
        assert self._sq != NULL

        cdef       int   status
        cdef const char* acc    = accession

        with nogil:
            status = libeasel.sq.esl_sq_SetAccession(self._sq, acc)
        if status == libeasel.eslEMEM:
            raise AllocationError("char", sizeof(char), len(accession))
        elif status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_sq_SetAccession")

    @property
    def name(self):
        """`bytes`: The name of the sequence.
        """
        assert self._sq != NULL
        return <bytes> self._sq.name

    @name.setter
    def name(self, bytes name not None):
        assert self._sq != NULL

        cdef       int   status
        cdef const char* nm     = name

        with nogil:
            status = libeasel.sq.esl_sq_SetName(self._sq, nm)
        if status == libeasel.eslEMEM:
            raise AllocationError("char", sizeof(char), len(name))
        elif status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_sq_SetName")

    @property
    def description(self):
        """`bytes`: The description of the sequence.
        """
        assert self._sq != NULL
        return <bytes> self._sq.desc

    @description.setter
    def description(self, bytes description):
        assert self._sq != NULL

        cdef       int   status
        cdef const char* desc   = description

        with nogil:
            status = libeasel.sq.esl_sq_SetDesc(self._sq, desc)
        if status == libeasel.eslEMEM:
            raise AllocationError("char", sizeof(char), len(description))
        elif status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_sq_SetDesc")

    @property
    def source(self):
        """`bytes`: The source of the sequence, if any.
        """
        assert self._sq != NULL
        return <bytes> self._sq.source

    @source.setter
    def source(self, bytes source):
        assert self._sq != NULL

        cdef       int   status
        cdef const char* src    = source

        with nogil:
            status = libeasel.sq.esl_sq_SetSource(self._sq, src)
        if status == libeasel.eslEMEM:
            raise AllocationError("char", sizeof(char), len(source))
        elif status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_sq_SetSource")

    @property
    def residue_markups(self):
        """`dict`: Extra residue markups, mapping information to each position.

        Keys and values are not decoded, since they are not necessarily
        valid UTF-8 bytestrings.

        Caution:
            The values of the dictionary must be the same size as the sequence
            itself. Trying to set a residue markup of the wrong length will
            raise a `ValueError`::

                >>> seq = TextSequence(sequence="TTAATTGGT")
                >>> seq.residue_markups = {b"quality": b"efcfffffcfee"}
                Traceback (most recent call last):
                  ...
                ValueError: Residue markup annotation has an invalid length (expected 9, got 12)

        .. versionadded:: 0.4.6

        """
        assert self._sq != NULL

        cdef int    i
        cdef bytes  tag
        cdef bytes  val
        cdef dict   xr  = {}
        cdef size_t off = 0 if libeasel.sq.esl_sq_IsText(self._sq) else 1

        for i in range(self._sq.nxr):
            tag = PyBytes_FromString(self._sq.xr_tag[i])
            val = PyBytes_FromStringAndSize(&self._sq.xr[i][off], self._sq.n)
            xr[tag] = val

        return xr

    @residue_markups.setter
    def residue_markups(self, dict xr):
        assert self._sq != NULL

        cdef const unsigned char[::1] tag
        cdef const unsigned char[::1] val

        cdef int     i
        cdef ssize_t xrlen = len(xr)
        cdef size_t  off   = 0 if libeasel.sq.esl_sq_IsText(self._sq) else 1
        cdef size_t  xrdim = self._sq.n + 2 * libeasel.sq.esl_sq_IsDigital(self._sq)

        # check the values have the right length before doing anything
        for tag, val in xr.items():
            if val.shape[0] != self._sq.n:
                raise ValueError(f"Residue markup annotation has an invalid length (expected {self._sq.n}, got {val.shape[0]})")
        # clear old values
        for i in range(self._sq.nxr):
            free(self._sq.xr_tag[i])
            self._sq.xr_tag[i] = NULL
            free(self._sq.xr[i])
            self._sq.xr[i] = NULL
        # reallocate arrays if needed
        if xrlen != self._sq.nxr:
            self._sq.nxr = xrlen
            self._sq.xr = <char**> realloc(<void*> self._sq.xr, xrlen * sizeof(char*))
            self._sq.xr_tag = <char**> realloc(<void*> self._sq.xr_tag, xrlen * sizeof(char*))
            if self._sq.xr == NULL or self._sq.xr_tag == NULL:
                raise AllocationError("char*", sizeof(char*), xrlen)
        # assign the new values
        for i, (tag, val) in enumerate(xr.items()):
            self._sq.xr_tag[i] = strdup(<const char*> &tag[0])
            self._sq.xr[i] = <char*> calloc(xrdim, sizeof(char))
            if self._sq.xr_tag[i] == NULL or self._sq.xr[i] == NULL:
                raise AllocationError("char", sizeof(char), xrdim)
            memcpy(&self._sq.xr[i][off], &val[0], val.shape[0] * sizeof(unsigned char))

    # --- Abstract methods ---------------------------------------------------

    def copy(self):
        """Duplicate the sequence, and return the copy.
        """
        raise NotImplementedError("Sequence.copy")

    # --- Methods ------------------------------------------------------------

    cpdef uint32_t checksum(self):
        """Calculate a 32-bit checksum for the sequence.
        """
        assert self._sq != NULL

        cdef int      status
        cdef uint32_t checksum = 0

        with nogil:
            status = libeasel.sq.esl_sq_Checksum(self._sq, &checksum)
        if status == libeasel.eslOK:
            return checksum
        else:
            raise UnexpectedError(status, "esl_sq_Checksum")

    cpdef void clear(self) except *:
        """Reinitialize the sequence for re-use.
        """
        assert self._sq != NULL

        cdef int status

        with nogil:
            status = libeasel.sq.esl_sq_Reuse(self._sq)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_sq_Reuse")

    cpdef void write(self, object fh) except *:
        """Write the sequence alignement to a file handle, in FASTA format.

        Arguments:
            fh (`io.IOBase`): A Python file handle, opened in binary mode.

        .. versionadded:: 0.3.0

        """
        assert self._sq != NULL

        cdef int    status
        cdef FILE*  file   = NULL

        try:
            file = fopen_obj(fh, "w")
            status = libeasel.sqio.ascii.esl_sqascii_WriteFasta(file, self._sq, False)
        finally:
            if file is not NULL:
                fclose(file)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_sqascii_WriteFasta")


cdef class TextSequence(Sequence):
    """A biological sequence stored in text mode.

    Hint:
        Use the ``sequence`` property to access the sequence letters as a
        Python string.

    .. versionadded:: 0.10.4
        `pickle` protocol support.

    """

    @classmethod
    def sample(
        cls,
        int max_length,
        RandomnessOrSeed randomness = None,
    ):
        """Sample a sequence of length at most ``L`` at random.

        Arguments:
            max_length (`int`): The maximum length of the sequence to
                generate (the actual sequence length is sampled).
            randomness (`~pyhmmer.easel.Randomness`, `int` or `None`): The
                random number generator to use for sampling, or a seed to
                initialize a generator. If `None` or ``0`` given, create
                a new random number generator with a random seed.

        Returns:
            `~pyhmmer.easel.TextSequence`: A new text sequence generated
            at random.

        Hint:
            This constructor is only useful for testing and should not be
            used to generate random sequences to e.g. compute a background
            distribution for a statistical method, since this function
            samples alphabet residues at random irrespective of prior
            frequences.

        .. versionadded:: 0.11.1

        """
        cdef int          status
        cdef Randomness   rng
        cdef TextSequence seq    = TextSequence.__new__(TextSequence)

        if RandomnessOrSeed is Randomness:
            rng = randomness
        else:
            rng = Randomness(randomness)

        status = libeasel.sq.esl_sq_Sample(
            rng._rng,
            NULL,
            max_length,
            &seq._sq
        )
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_sq_Sample")

        return seq

    # --- Magic methods ------------------------------------------------------

    def __init__(
        self,
        *args,
        bytes name=None,
        bytes description=None,
        bytes accession=None,
        str   sequence=None,
        bytes source=None,
        dict  residue_markups=None,
    ):
        """__init__(self, name=None, description=None, accession=None, sequence=None, source=None, residue_markups=None)\n--\n

        Create a new text-mode sequence with the given attributes.

        .. versionadded:: 0.10.4
            The ``residue_markups`` argument.

        .. deprecated:: 0.11.1
            Passing positional arguments to constructor.

        """
        cdef bytes sq

        # TODO: Remove in 0.12.0 (deprecation)
        if len(args) > 0:
            warnings.warn(
                "TextSequence.__init__ will not accept positional arguments after v0.12.0",
                category=DeprecationWarning
            )
            if len(args) > 0:
                name = args[0]
            if len(args) > 1:
                description = args[1]
            if len(args) > 2:
                accession = args[2]
            if len(args) > 3:
                sequence = args[3]
            if len(args) > 4:
                source = args[4]
            if len(args) > 5:
                residue_markups = args[5]

        if sequence is not None:
            sq = sequence.encode("ascii")
            self._sq = libeasel.sq.esl_sq_CreateFrom(NULL, sq, NULL, NULL, NULL)
        else:
            self._sq = libeasel.sq.esl_sq_Create()
        if self._sq == NULL:
            raise AllocationError("ESL_SQ", sizeof(ESL_SQ))
        self._sq.abc = NULL

        if name is not None:
            self.name = name
        if accession is not None:
            self.accession = accession
        if description is not None:
            self.description = description
        if source is not None:
            self.source = source
        if residue_markups is not None:
            self.residue_markups = residue_markups

        assert libeasel.sq.esl_sq_IsText(self._sq)
        assert self._sq.name != NULL
        assert self._sq.desc != NULL
        assert self._sq.acc != NULL

    def __reduce__(self):
        constructor = functools.partial(
            type(self),
            name=self.name,
            description=self.description,
            accession=self.accession,
            sequence=self.sequence,
            source=self.source,
            residue_markups=self.residue_markups,
        )
        return constructor, ()

    def __getbuffer__(self, Py_buffer* buffer, int flags):
        assert self._sq != NULL

        if flags & PyBUF_FORMAT:
            buffer.format = 'B'
        else:
            buffer.format = NULL
        buffer.internal = NULL
        buffer.itemsize = sizeof(unsigned char)
        buffer.ndim = 1
        PyBuffer_FillInfo(buffer, self, self._sq.seq, self._sq.n * sizeof(unsigned char), True, flags)

    # --- Properties ---------------------------------------------------------

    @property
    def sequence(self):
        """`str`: The raw sequence letters, as a Python string.
        """
        assert self._sq != NULL
        assert libeasel.sq.esl_sq_IsText(self._sq)
        return PyUnicode_DecodeASCII(self._sq.seq, self._sq.n, NULL)

    # --- Methods ------------------------------------------------------------

    cpdef DigitalSequence digitize(self, Alphabet alphabet):
        """Convert the text sequence to a digital sequence using ``alphabet``.

        Returns:
            `DigitalSequence`: A copy of the sequence in digital mode,
            digitized with ``alphabet``.

        """
        assert self._sq != NULL
        assert libeasel.sq.esl_sq_IsText(self._sq)

        cdef int             status
        cdef ESL_ALPHABET*   abc    = alphabet._abc
        cdef DigitalSequence new    = DigitalSequence.__new__(DigitalSequence, alphabet)

        with nogil:
            new._sq = libeasel.sq.esl_sq_CreateDigital(abc)
            if new._sq == NULL:
                raise AllocationError("ESL_SQ", sizeof(ESL_SQ))
            status = libeasel.sq.esl_sq_Copy(self._sq, new._sq)

        if status == libeasel.eslOK:
            assert libeasel.sq.esl_sq_IsDigital(new._sq)
            return new
        elif status == libeasel.eslEINVAL:
            raise ValueError(f"Cannot digitize sequence with alphabet {alphabet}: invalid chars in sequence")
        else:
            raise UnexpectedError(status, "esl_sq_Copy")

    cpdef TextSequence copy(self):
        """Duplicate the text sequence, and return the copy.
        """
        assert self._sq != NULL
        assert libeasel.sq.esl_sq_IsText(self._sq)

        cdef int          status
        cdef TextSequence new    = TextSequence.__new__(TextSequence)

        with nogil:
            new._sq = libeasel.sq.esl_sq_Create()
            if new._sq == NULL:
                raise AllocationError("ESL_SQ", sizeof(ESL_SQ))
            new._sq.abc = NULL

            status = libeasel.sq.esl_sq_Copy(self._sq, new._sq)
            if status != libeasel.eslOK:
                raise UnexpectedError(status, "esl_sq_Copy")

        assert libeasel.sq.esl_sq_IsText(new._sq)
        return new

    cpdef TextSequence reverse_complement(self, bint inplace=False):
        """Build the reverse complement of the sequence.

        This method assumes that the sequence alphabet is IUPAC/DNA. If the
        sequence contains any unknown letters, they will be replaced by
        :math:`N` in the reverse-complement.

        Arguments:
            inplace (`bool`): Whether or not to copy the sequence before
                computing its reverse complement. With `False` (the default),
                the method will return a copy of the sequence that has been
                reverse-complemented. With `True`, it will reverse-complement
                inplace and return `None`.

        Raises:
            `UserWarning`: When the sequence contains unknown characters.

        Example:
            >>> seq = TextSequence(sequence="ATGC")
            >>> seq.reverse_complement().sequence
            'GCAT'

        Caution:
            The copy made when ``inplace`` is `False` is an exact copy, so
            the `name`, `description` and `accession` of the copy will be
            the same. This could lead to duplicates if you're not careful!

        .. versionadded:: 0.3.0

        """
        assert self._sq != NULL

        cdef TextSequence rc
        cdef int          status

        if inplace:
            status = libeasel.sq.esl_sq_ReverseComplement(self._sq)
        else:
            rc = self.copy()
            status = libeasel.sq.esl_sq_ReverseComplement(rc._sq)

        if status == libeasel.eslEINVAL:
            warnings.warn(
                "reverse-complementing a text sequence with non-DNA characters",
                UserWarning,
                stacklevel=2
            )
        elif status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_sq_ReverseComplement")

        return None if inplace else rc


cdef class DigitalSequence(Sequence):
    """A biological sequence stored in digital mode.

    Attributes:
        alphabet (`Alphabet`, *readonly*): The biological alphabet used to
            encode this sequence to digits.

    Hint:
        Use the ``sequence`` property to access the sequence digits as a
        memory view, allowing to access the individual bytes. This can be
        combined with `numpy.asarray` to get the sequence as an array with
        zero-copy.

    .. versionadded:: 0.10.4
        `pickle` protocol support.

    """

    @classmethod
    def sample(
        cls,
        Alphabet alphabet not None,
        int max_length,
        RandomnessOrSeed randomness = None,
    ):
        """Sample a sequence of length at most ``L`` at random.

        Arguments:
            alphabet (`~pyhmmer.easel.Alphabet`): The alphabet of the
                sequence.
            max_length (`int`): The maximum length of the sequence to
                generate (the actual sequence length is sampled).
            randomness (`~pyhmmer.easel.Randomness`, `int` or `None`): The
                random number generator to use for sampling, or a seed to
                initialize a generator. If `None` or ``0`` given, create
                a new random number generator with a random seed.

        Returns:
            `~pyhmmer.easel.DigitalSequence`: A new digital sequence
            generated at random, including degenerate symbols.

        Hint:
            This constructor is only useful for testing and should not be
            used to generate random sequences to e.g. compute a background
            distribution for a statistical method, since this function
            samples alphabet residues at random irrespective of prior
            frequences.

        .. versionadded:: 0.11.1

        """
        cdef int             status
        cdef Randomness      rng
        cdef DigitalSequence seq    = DigitalSequence.__new__(DigitalSequence, alphabet)

        if RandomnessOrSeed is Randomness:
            rng = randomness
        else:
            rng = Randomness(randomness)

        status = libeasel.sq.esl_sq_Sample(
            rng._rng,
            alphabet._abc,
            max_length,
            &seq._sq
        )
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_sq_Sample")

        return seq

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self, Alphabet alphabet, *args, **kwargs):
        self.alphabet = alphabet

    def __init__(self,
              Alphabet              alphabet      not None,
              *args,
              bytes                 name            = None,
              bytes                 description     = None,
              bytes                 accession       = None,
        const libeasel.ESL_DSQ[::1] sequence        = None,
              bytes                 source          = None,
              dict                  residue_markups = None,
    ):
        """__init__(self, alphabet, name=None, description=None, accession=None, sequence=None, source=None, residue_markups=None)\n--\n

        Create a new digital-mode sequence with the given attributes.

        Raises:
            `ValueError`: When ``sequence`` contains digits outside the
                alphabet symbol range.

        .. versionadded:: 0.1.4

        .. versionadded:: 0.10.4
            The ``residue_markups`` argument.

        .. deprecated:: 0.11.1
            Passing positional arguments other than ``alphabet``.

        """
        cdef int     status
        cdef int64_t i
        cdef int64_t n

        # TODO: Remove in 0.12.0 (deprecation)
        if len(args) > 0:
            warnings.warn(
                "DigitalSequence.__init__ will not accept positional arguments besides `alphabet` after v0.12.0",
                category=DeprecationWarning
            )
            if len(args) > 0:
                name = args[0]
            if len(args) > 1:
                description = args[1]
            if len(args) > 2:
                accession = args[2]
            if len(args) > 3:
                sequence = args[3]
            if len(args) > 4:
                source = args[4]
            if len(args) > 5:
                residue_markups = args[5]

        # create an empty digital sequence
        self._sq = libeasel.sq.esl_sq_CreateDigital(alphabet._abc)
        if self._sq == NULL:
            raise AllocationError("ESL_SQ", sizeof(ESL_SQ))

        # NB: because the easel sequence has sentinel bytes that we hide from
        #     the user, we cannot just copy the sequence here or use the libeasel
        #     internals; instead, if a sequence is given, we need to emulate
        #     the `esl_sq_CreateDigitalFrom` but copy the sequence with different
        #     offsets.
        if sequence is not None:
            # we can release the GIL while copying memory
            with nogil:
                n = sequence.shape[0]
                # check the sequence encoding is compatible with the alphabet
                for i in range(n):
                    if not libeasel.alphabet.esl_abc_XIsValid(alphabet._abc, sequence[i]):
                        raise ValueError(f"Invalid alphabet character in digital sequence: {sequence[i]}")
                # grow the sequence buffer so it can hold `n` residues
                status = libeasel.sq.esl_sq_GrowTo(self._sq, n)
                if status != libeasel.eslOK:
                    raise UnexpectedError(status, "esl_sq_GrowTo")
                # update the digital sequence buffer
                self._sq.dsq[0] = self._sq.dsq[n+1] = libeasel.eslDSQ_SENTINEL
                memcpy(&self._sq.dsq[1], &sequence[0], n * sizeof(char))
            # set the coor bookkeeping like it would happen
            self._sq.start = 1
            self._sq.C = 0
            self._sq.end = self._sq.W = self._sq.L = self._sq.n = n

        if name is not None:
            self.name = name
        if accession is not None:
            self.accession = accession
        if description is not None:
            self.description = description
        if source is not None:
            self.source = source
        if residue_markups is not None:
            self.residue_markups = residue_markups

        assert libeasel.sq.esl_sq_IsDigital(self._sq)
        assert self._sq.name != NULL
        assert self._sq.desc != NULL
        assert self._sq.acc != NULL

    def __reduce__(self):
        constructor = functools.partial(
            type(self),
            name=self.name,
            description=self.description,
            accession=self.accession,
            sequence=self.sequence,
            source=self.source,
            residue_markups=self.residue_markups,
        )
        return constructor, (self.alphabet,)

    def __getbuffer__(self, Py_buffer* buffer, int flags):
        assert self._sq != NULL

        if flags & PyBUF_FORMAT:
            buffer.format = 'B'
        else:
            buffer.format = NULL
        buffer.internal = NULL
        buffer.itemsize = sizeof(ESL_DSQ)
        buffer.ndim = 1
        PyBuffer_FillInfo(buffer, self, &self._sq.dsq[1], self._sq.n * sizeof(ESL_DSQ), True, flags)

    # --- Properties ---------------------------------------------------------

    @property
    def sequence(self):
        """`VectorU8`: The raw sequence digits, as a byte vector.

        Note:
            The internal ``ESL_SQ`` object allocates a buffer of size
            :math:`n+2` (where :math:`n` is the number of residues in the
            sequence), with the first and the last element of the buffer
            being sentinel values. This vector does not expose the sentinel
            values, only the :math:`n` elements of the buffer in between.

        .. versionchanged:: 0.4.0
           Property is now a `VectorU8` instead of a memoryview.

        """
        assert self._sq != NULL

        cdef VectorU8 seq = VectorU8.__new__(VectorU8)
        seq._n = seq._shape[0] = self._sq.n
        seq._data = &self._sq.dsq[1]
        seq._owner = self

        return seq

    # --- Methods ------------------------------------------------------------

    cpdef DigitalSequence copy(self):
        """Duplicate the digital sequence, and return the copy.
        """
        assert self._sq != NULL
        assert libeasel.sq.esl_sq_IsDigital(self._sq)

        cdef int             status
        cdef ESL_ALPHABET*   abc    = self.alphabet._abc
        cdef DigitalSequence new    = DigitalSequence.__new__(DigitalSequence, self.alphabet)

        with nogil:
            new._sq = libeasel.sq.esl_sq_CreateDigital(abc)
            if new._sq == NULL:
                raise AllocationError("ESL_SQ", sizeof(ESL_SQ))

            status = libeasel.sq.esl_sq_Copy(self._sq, new._sq)
            if status != libeasel.eslOK:
                raise UnexpectedError(status, "esl_sq_Copy")

        assert libeasel.sq.esl_sq_IsDigital(new._sq)
        return new

    cpdef TextSequence textize(self):
        """Convert the digital sequence to a text sequence.

        Returns:
            `TextSequence`: A copy of the sequence in text-mode.

        .. versionadded:: 0.1.4

        """
        assert self._sq != NULL
        assert libeasel.sq.esl_sq_IsDigital(self._sq)

        cdef int          status
        cdef TextSequence new    = TextSequence.__new__(TextSequence)

        with nogil:
            new._sq = libeasel.sq.esl_sq_Create()

            if new._sq == NULL:
                raise AllocationError("ESL_SQ", sizeof(ESL_SQ))

            status = libeasel.sq.esl_sq_Copy(self._sq, new._sq)
            if status != libeasel.eslOK:
                raise UnexpectedError(status, "esl_sq_Copy")

        assert libeasel.sq.esl_sq_IsText(new._sq)
        return new

    cpdef DigitalSequence translate(self, GeneticCode genetic_code = GeneticCode()):
        """Translate the sequence using the given genetic code.

        Arguments:
            genetic_code (`~pyhmmer.easel.GeneticCode`): The genetic code to
                use for translating the sequence. If none provided, the
                default uses the standard translation table (1) and expects
                DNA sequences.

        Returns:
            `~pyhmmer.easel.DigitalSequence`: The translation of the
            input sequence, in digital mode.

        Raises:
            `pyhmmer.errors.AlphabetMismatch`: When the ``genetic_code``
                expects a different nucleotide alphabet than the one
                currently in use to encode the sequence.
            `ValueError`: When ``sequence`` could not be translated
                properly, because of a codon could not be recognized, or
                because the sequence has an invalid length.

        Note:
            The translation of a DNA/RNA codon supports ambiguous codons.
            If the amino acid is unambiguous, despite codon ambiguity,
            the correct amino acid is still determined: ``GGR`` translates
            as ``Gly``, ``UUY`` as ``Phe``, etc. If there is no single
            unambiguous amino acid translation, the codon is translated
            as ``X``. Ambiguous amino acids (such as ``J`` or ``B``) are
            never produced.

        .. versionadded:: 0.7.2

        """
        assert self._sq != NULL
        assert self.alphabet is not None

        cdef DigitalSequence protein
        cdef int64_t         ntlen   = len(self)
        cdef int64_t         aalen   = ntlen // 3

        # check sequence can be translated
        if not self.alphabet._eq(genetic_code.nucleotide_alphabet):
            raise AlphabetMismatch(genetic_code.nucleotide_alphabet, self.alphabet)
        if ntlen % 3 != 0:
            raise ValueError(f"Incomplete sequence of length {len(self)!r}")

        # create copy with metadata
        protein = DigitalSequence(
            genetic_code.amino_alphabet,
            name=self.name,
            description=self.description,
            accession=self.accession,
            source=self.source,
        )

        # allocate output buffer
        status = libeasel.sq.esl_sq_GrowTo(protein._sq, aalen)
        if status == libeasel.eslEMEM:
            raise AllocationError("ESL_DSQ", sizeof(ESL_DSQ), aalen + 2)
        elif status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_sq_Grow")

        # translate & rename sequence
        with nogil:
            genetic_code._translate(&self._sq.dsq[1], ntlen, &protein._sq.dsq[1], aalen)
            protein._sq.dsq[0] = protein._sq.dsq[aalen+1] = libeasel.eslDSQ_SENTINEL

        # record sequence coordinates
        protein._sq.start = 1
        protein._sq.C = 0
        protein._sq.end = protein._sq.W = protein._sq.L = protein._sq.n = aalen

        return protein

    cpdef DigitalSequence reverse_complement(self, bint inplace=False):
        """Build the reverse complement of the sequence.

        Arguments:
            inplace (`bool`): Whether or not to copy the sequence before
                computing its reverse complement. With `False` (the default),
                the method will return a copy of the sequence that has been
                reverse-complemented. With `True`, it will reverse-complement
                inplace and return `None`.

        Raises:
            `ValueError`: When the alphabet of the `DigitalSequence` does
                not have a complement mapping set (e.g., `Alphabet.amino`).

        Caution:
            The copy made when ``inplace`` is `False` is an exact copy, so
            the `name`, `description` and `accession` of the copy will be
            the same. This could lead to duplicates if you're not careful!

        .. versionadded:: 0.3.0

        """
        assert self._sq != NULL
        assert self.alphabet is not None

        cdef DigitalSequence rc
        cdef int             status

        if self.alphabet._abc.complement == NULL:
            raise ValueError(f"{self.alphabet} has no defined complement")

        if inplace:
            status = libeasel.sq.esl_sq_ReverseComplement(self._sq)
        else:
            rc = self.copy()
            status = libeasel.sq.esl_sq_ReverseComplement(rc._sq)

        if status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_sq_ReverseComplement")

        return None if inplace else rc


# --- Sequence Block ---------------------------------------------------------

class _SequenceBlockIndex(collections.abc.Mapping):
    """A read-only mapping of sequence names to sequences of a block.
    """
    __slots__ = ("block",)

    def __init__(self, SequenceBlock block):
        self.block = block
        if len(block._indexed) != len(block):
            with nogil:
                block._rehash()

    def __len__(self):
        cdef SequenceBlock block = self.block
        return libeasel.keyhash.esl_keyhash_GetNumber(block._indexed._kh)

    def __getitem__(self, object item):
        cdef int                      status
        cdef int                      index  = -1
        cdef const unsigned char[::1] key    = item
        cdef esl_pos_t                length = key.shape[0]
        cdef SequenceBlock            block  = self.block

        assert block._indexed is not None
        assert block._indexed._kh != NULL

        with nogil:
            status = libeasel.keyhash.esl_keyhash_Lookup(
                block._indexed._kh,
                <const char*> &key[0],
                length,
                &index
            )
        if status == libeasel.eslOK:
            return block[index]
        elif status == libeasel.eslENOTFOUND:
            raise KeyError(item)
        else:
            raise UnexpectedError(status, "esl_keyhash_Lookup")

    def __iter__(self):
        cdef size_t        i
        cdef size_t        length
        cdef SequenceBlock block  = self.block

        assert block._indexed is not None
        assert block._indexed._kh != NULL

        length = libeasel.keyhash.esl_keyhash_GetNumber(block._indexed._kh)
        for i in range(length):
            yield <bytes> libeasel.keyhash.esl_keyhash_Get(block._indexed._kh, i)


cdef class SequenceBlock:
    """An abstract container for storing `Sequence` objects.

    To pass the target sequences efficiently in `Pipeline.search_hmm`,
    an array is allocated so that the inner loop can iterate over the
    target sequences without having to acquire the GIL for each new
    sequence (this gave a huge performance boost in v0.4.5). However,
    there was no way to reuse this between different queries; some memory
    recycling was done, but the target sequences had to be indexed for
    every query. This class allows synchronizing a Python `list` of
    `Sequence` objects with an internal C-contiguous buffer of pointers
    to ``ESL_SQ`` structs that can be used in the HMMER search loop.

    .. versionadded:: 0.7.0

    """

    # NOTE(@althonos): This implementation of a sequence block doesn't
    #                  actually use an `ESL_SQ_BLOCK` because `ESL_SQ_BLOCK`
    #                  uses a contiguous array to store sequence data, which
    #                  is unpractical for `Sequence` objects that may be
    #                  held from elsewhere with reference counting. An array
    #                  of pointers is easier to work with in that particular
    #                  case.

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._refs = NULL
        self._length = 0
        self._capacity = 0
        self._storage = []
        self._largest = -1
        self._owner = None
        self._indexed = KeyHash()

    def __init__(self):
        raise TypeError("Can't instantiate abstract class 'SequenceBlock'")

    def __dealloc__(self):
        free(self._refs)

    def __len__(self):
        return self._length

    def __getitem__(self, object index):
        if isinstance(index, slice):
            return type(self)(self._storage[index])
        else:
            return self._storage[index]

    def __delitem__(self, object index):
        cdef size_t   i
        cdef Sequence sequence

        self._on_modification()

        if isinstance(index, slice):
            del self._storage[index]
            self._length = len(self._storage)
            self._allocate(self._length)
            for i, sequence in enumerate(self._storage):
                self._refs[i] = sequence._sq
        else:
            self.pop(index)

    def __reduce__(self):
        return type(self), (), None, iter(self)

    def __copy__(self):
        return self.copy()

    def __eq__(self, object other):
        cdef size_t         i
        cdef SequenceBlock  other_
        cdef const ESL_SQ** r1
        cdef const ESL_SQ** r2
        cdef bint           equal  = True

        if not isinstance(other, SequenceBlock):
            return NotImplemented

        other_ = other
        if self._length != other_._length:
            return False

        with nogil:
            for i in range(self._length):
                status = libeasel.sq.esl_sq_Compare(self._refs[i], other_._refs[i])
                if status == libeasel.eslOK:
                    continue
                elif status == libeasel.eslFAIL:
                    equal = False
                    break
                else:
                    raise UnexpectedError(status, "esl_sq_Compare")

        return equal

    def __sizeof__(self):
        return sizeof(self) + self._capacity * sizeof(ESL_SQ*)

    # --- Properties ---------------------------------------------------------

    @property
    def indexed(self):
        """`~collections.abc.Mapping`: A mapping of names to sequences.

        This property can be used to access the sequence of a sequence block
        by name. An index is created the first time this property is accessed.
        An error is raised if the block contains duplicate sequence names.

        Raises:
            `KeyError`: When attempting to create an index for an alignment
                containing duplicate sequence names.

        Example:
            >>> s1 = TextSequence(name=b"seq1", sequence="ATGC")
            >>> s2 = TextSequence(name=b"seq2", sequence="ATTA")
            >>> block = TextSequenceBlock([s1, s2])
            >>> block.indexed[b'seq1'].sequence
            'ATGC'
            >>> block.indexed[b'seq3']
            Traceback (most recent call last):
            ...
            KeyError: b'seq3'

        .. versionadded:: 0.11.1

        """
        return _SequenceBlockIndex(self)

    # --- C methods ----------------------------------------------------------

    cdef void _allocate(self, size_t n) except *:
        """Allocate enough storage for at least ``n`` items.
        """
        cdef size_t i
        cdef size_t capacity = new_capacity(n, self._length)
        with nogil:
            self._refs = <ESL_SQ**> realloc(self._refs, capacity * sizeof(ESL_SQ*))
        if self._refs == NULL:
            self._capacity = 0
            raise AllocationError("ESL_SQ*", sizeof(ESL_SQ*), capacity)
        else:
            self._capacity = capacity
        # for i in range(self._length, self._capacity):
        #     self._refs[i] = NULL

    cdef void _on_modification(self) noexcept:
        self._largest = -1 # invalidate cache
        self._indexed.clear()

    cdef int _rehash(self) except 1 nogil:
        assert self._indexed is not None
        assert self._indexed._kh != NULL

        cdef size_t      idx
        cdef const char* name
        cdef int         status

        status = libeasel.keyhash.esl_keyhash_Reuse(self._indexed._kh)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_keyhash_Reuse")

        for idx in range(self._length):
            status = libeasel.keyhash.esl_keyhash_Store(
                self._indexed._kh,
                self._refs[idx].name,
                -1,
                NULL,
            )
            if status == libeasel.eslEDUP:
                raise KeyError("duplicate keys in sequence block")
            elif status != libeasel.eslOK:
                raise UnexpectedError(status, "esl_keyhash_Store")

        return 0

    # --- Python methods -----------------------------------------------------

    cdef void _append(self, Sequence sequence) except *:
        if self._length == self._capacity:
            self._allocate(self._length + 1)
        self._storage.append(sequence)
        self._refs[self._length] = sequence._sq
        self._length += 1
        self._on_modification()

    cpdef void clear(self) except *:
        """Remove all sequences from the block.
        """
        cdef size_t i
        self._storage.clear()
        self._length = 0
        self._on_modification()

    cpdef void extend(self, object iterable) except *:
        """Extend block by appending sequences from the iterable.
        """
        cdef size_t hint = operator.length_hint(iterable)
        if self._length + hint > self._capacity:
            self._allocate(self._length + hint)
        self._on_modification()
        for sequence in iterable:
            self.append(sequence)

    cdef Sequence _pop(self, ssize_t index=-1):
        cdef ssize_t index_ = index
        if self._length == 0:
            raise IndexError("pop from empty block")
        if index_ < 0:
            index_ += self._length
        if index_ < 0 or <size_t> index_ >= self._length:
            raise IndexError(index)

        # remove item from storage
        item = self._storage.pop(index_)

        # update pointers in the reference array
        self._length -= 1
        if <size_t> index_ < self._length:
            memmove(&self._refs[index_], &self._refs[index_ + 1], (self._length - index_)*sizeof(ESL_SQ*))

        self._on_modification()
        return item

    cdef void _insert(self, ssize_t index, Sequence sequence) except *:
        if index < 0:
            index = 0
        elif <size_t> index > self._length:
            index = self._length

        if self._length == self._capacity - 1:
            self._allocate(self._capacity + 1)

        if <size_t> index != self._length:
            memmove(&self._refs[index + 1], &self._refs[index], (self._length - index)*sizeof(ESL_SQ*))

        self._storage.insert(index, sequence)
        self._refs[index] = sequence._sq
        self._length += 1
        self._on_modification()

    cdef size_t _index(self, Sequence sequence, ssize_t start=0, ssize_t stop=sys.maxsize) except *:
        cdef size_t i
        cdef size_t start_
        cdef size_t stop_
        cdef int    status

        # wrap once is negative indices are used
        if start < 0:
            start += <ssize_t> self._length
        if stop < 0:
            stop += <ssize_t> self._length

        # wrap a second time if indices are still negative or out of bounds
        stop_ = min(stop, <ssize_t> self._length)
        start_ = max(start, 0)

        # scan to locate the sequence
        with nogil:
            for i in range(start_, stop_):
                status = libeasel.sq.esl_sq_Compare(sequence._sq, self._refs[i])
                if status == libeasel.eslOK:
                    break
                elif status != libeasel.eslFAIL:
                    raise UnexpectedError(status, "esl_sq_Compare")
            else:
                raise ValueError(f"sequence {sequence.name!r} not in block")

        return i

    cdef void _remove(self, Sequence sequence) except *:
        self.pop(self._index(sequence))

    cpdef Sequence largest(self):
        """Return the largest sequence in the block.
        """
        cdef size_t i

        if self._length == 0:
            raise ValueError("block is empty")
        if self._largest == -1:
            self._largest = 0
            for i in range(1, self._length):
                if self._refs[i].L > self._refs[self._largest].L:
                    self._largest = i

        return self._storage[self._largest]

    cpdef SequenceBlock copy(self):
        """Return a copy of the sequence block.

        Note:
            The sequence internally refered to by this collection are not
            copied. Use `copy.deepcopy` if you also want to duplicate the
            internal storage of each sequence.

        """
        return self[:]


cdef class TextSequenceBlock(SequenceBlock):
    """A container for storing `TextSequence` objects.

    .. versionadded:: 0.7.0

    .. versionadded:: 0.10.4
        `pickle` protocol support.

    """

    # --- Magic methods ------------------------------------------------------

    def __init__(self, object iterable = ()):
        """__init__(self, iterable=())\n--\n

        Create a new block from an iterable of text sequences.

        """
        self.clear()
        self.extend(iterable)

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"{ty}({self._storage!r})"

    def __len__(self):
        return self._length

    def __getitem__(self, object index):
        if isinstance(index, slice):
            return type(self)(self._storage[index])
        else:
            return self._storage[index]

    def __setitem__(self, object index, object sequences):
        cdef size_t       i
        cdef TextSequence sequence

        if isinstance(index, slice):
            self._storage[index] = sequences
            self._length = len(self._storage)
            self._allocate(self._length)
            for i, sequence in enumerate(self._storage):
                self._refs[i] = sequence._sq
        else:
            sequence = sequences
            self._storage[index] = sequence
            self._refs[index] = sequence._sq

    def __contains__(self, object sequence):
        if isinstance(sequence, TextSequence):
            try:
                return self._index(sequence) >= 0
            except ValueError:
                return False
        return False

    def __reduce__(self):
        return type(self), (), None, iter(self)

    # --- Python methods -----------------------------------------------------

    cpdef void append(self, TextSequence sequence) except *:
        """Append ``sequence`` at the end of the block.
        """
        self._append(sequence)

    cpdef TextSequence pop(self, ssize_t index=-1):
        """Remove and return a sequence from the block (the last one by default).
        """
        return self._pop(index)

    cpdef void insert(self, ssize_t index, TextSequence sequence) except *:
        """Insert a new sequence in the block before ``index``.
        """
        self._insert(index, sequence)

    cpdef size_t index(self, TextSequence sequence, ssize_t start=0, ssize_t stop=sys.maxsize) except *:
        """Return the index of the first occurence of ``sequence``.

        Raises:
            `ValueError`: When the block does not contain ``sequence``.

        """
        return self._index(sequence, start, stop)

    cpdef void remove(self, TextSequence sequence) except *:
        """Remove the first occurence of the given sequence.
        """
        self._remove(sequence)

    cpdef DigitalSequenceBlock digitize(self, Alphabet alphabet):
        """Create a block containing sequences from this block in digital mode.
        """
        cdef size_t               i
        cdef list                 seqs  = [ DigitalSequence(alphabet) for _ in range(self._length) ]
        cdef DigitalSequenceBlock block = DigitalSequenceBlock(alphabet, seqs)

        with nogil:
            for i in range(self._length):
                libeasel.sq.esl_sq_Copy(self._refs[i], block._refs[i])

        return block

    cpdef TextSequence largest(self):
        """Return the largest sequence in the block.
        """
        return SequenceBlock.largest(self)

    cpdef TextSequenceBlock copy(self):
        """Return a copy of the text sequence block.

        Note:
            The sequence internally refered to by this collection are not
            copied. Use `copy.deepcopy` is you also want to duplicate the
            internal storage of each sequence.

        """
        cdef TextSequenceBlock new = TextSequenceBlock.__new__(TextSequenceBlock)
        new._storage = self._storage.copy()
        new._length = self._length
        new._largest = self._largest
        new._owner = self._owner
        new._allocate(self._length)
        memcpy(new._refs, self._refs, self._length * sizeof(ESL_SQ*))
        return new


cdef class DigitalSequenceBlock(SequenceBlock):
    """A container for storing `DigitalSequence` objects.

    Attributes:
        alphabet (`Alphabet`, *readonly*): The biological alphabet shared by
            all sequences in the collection.

    .. versionadded:: 0.7.0

    .. versionadded:: 0.10.4
        `pickle` protocol support.

    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self, Alphabet alphabet, *args, **kwargs):
        self.alphabet = alphabet

    def __init__(self, Alphabet alphabet not None, object iterable = ()):
        """__init__(self, alphabet, iterable=())\n--\n

        Create a new digital sequence block with the given alphabet.

        Arguments:
            alphabet (`~pyhmmer.easel.Alphabet`): The alphabet to use for all
                the sequences in the block.
            iterable (iterable of `~pyhmmer.easel.DigitalSequence`): An initial
                collection of digital sequences to add to the block.

        Raises:
            `~pyhmmer.easel.AlphabetMismatch`: When the alphabet of one of the
                sequences does not match ``alphabet``.

        """
        self.clear()
        self.extend(iterable)

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"{ty}({self.alphabet!r}, {self._storage!r})"

    def __len__(self):
        return self._length

    def __getitem__(self, object index):
        if isinstance(index, slice):
            return type(self)(self.alphabet, self._storage[index])
        else:
            return self._storage[index]

    def __reduce__(self):
        return type(self), (self.alphabet,), None, iter(self)

    def __setitem__(self, object index, object sequences):
        cdef size_t          i
        cdef DigitalSequence sequence

        if isinstance(index, slice):
            self._storage[index] = sequences
            self._length = len(self._storage)
            self._allocate(self._length)
            for i, sequence in enumerate(self._storage):
                if sequence.alphabet != self.alphabet:
                    raise AlphabetMismatch(self.alphabet, sequence.alphabet)
                self._refs[i] = sequence._sq
        else:
            sequence = sequences
            if sequence.alphabet != self.alphabet:
                raise AlphabetMismatch(self.alphabet, sequence.alphabet)
            self._storage[index] = sequence
            self._refs[index] = sequence._sq

    def __contains__(self, object sequence):
        if isinstance(sequence, DigitalSequence):
            try:
                return self._index(sequence) >= 0
            except ValueError:
                return False
        return False

    def __reduce__(self):
        return type(self), (self.alphabet,), None, iter(self)

    # --- Python methods -----------------------------------------------------

    cpdef void append(self, DigitalSequence sequence) except *:
        """Append ``sequence`` at the end of the block.
        """
        if sequence.alphabet != self.alphabet:
            raise AlphabetMismatch(self.alphabet, sequence.alphabet)
        self._append(sequence)

    cpdef DigitalSequence pop(self, ssize_t index=-1):
        """Remove and return a sequence from the block (the last one by default).
        """
        return self._pop(index)

    cpdef void insert(self, ssize_t index, DigitalSequence sequence) except *:
        """Insert a new sequence in the block before ``index``.
        """
        if sequence.alphabet != self.alphabet:
            raise AlphabetMismatch(self.alphabet, sequence.alphabet)
        self._insert(index, sequence)

    cpdef size_t index(self, DigitalSequence sequence, ssize_t start=0, ssize_t stop=sys.maxsize) except *:
        """Return the index of the first occurence of ``sequence``.

        Raises:
            `ValueError`: When the block does not contain ``sequence``.

        """
        if sequence.alphabet != self.alphabet:
            raise AlphabetMismatch(self.alphabet, sequence.alphabet)
        return self._index(sequence, start, stop)

    cpdef void remove(self, DigitalSequence sequence) except *:
        """Remove the first occurence of the given sequence.
        """
        if sequence.alphabet != self.alphabet:
            raise AlphabetMismatch(self.alphabet, sequence.alphabet)
        self._remove(sequence)

    cpdef TextSequenceBlock textize(self):
        """Create a block containing sequences from this block in text mode.
        """
        cdef size_t            i
        cdef list              seqs  = [ TextSequence() for _ in range(self._length) ]
        cdef TextSequenceBlock block = TextSequenceBlock(seqs)

        with nogil:
            for i in range(self._length):
                libeasel.sq.esl_sq_Copy(self._refs[i], block._refs[i])

        return block

    cpdef DigitalSequenceBlock translate(self, GeneticCode genetic_code = GeneticCode()):
        """Translate the sequence block using the given genetic code.

        Arguments:
            genetic_code (`~pyhmmer.easel.GeneticCode`): The genetic code to
                use for translating the sequence. If none provided, the
                default uses the standard translation table (1) and expects
                DNA sequences.

        Returns:
            `~pyhmmer.easel.DigitalSequenceBlock`: The translation of
            each sequence from the block, in digital mode.

        Raises:
            `~pyhmmer.errors.AlphabetMismatch`: When the ``genetic_code``
                expects a different nucleotide alphabet than the one
                currently for the sequences in the block.
            `ValueError`: When a sequence from the block could not be
                translated properly, because of a codon could not be
                recognized, or because the sequence has an invalid length.

        See Also:
            `DigitalSequence.translate` for more information on how
            ambiguous nucleotides are handled.

        """
        assert self.alphabet is not None

        # cdef int64_t         ntlen   = len(self)
        # cdef int64_t         aalen   = ntlen // 3
        # cdef DigitalSequence protein = DigitalSequence(genetic_code.amino_alphabet)

        cdef size_t               i
        cdef int64_t              aalen
        cdef int64_t              ntlen
        cdef int                  status
        cdef DigitalSequence      protein
        cdef DigitalSequenceBlock proteins

        # check block can be translated
        if not self.alphabet._eq(genetic_code.nucleotide_alphabet):
            raise AlphabetMismatch(genetic_code.nucleotide_alphabet, self.alphabet)

        # pre-allocate protein sequence buffers
        proteins = DigitalSequenceBlock(genetic_code.amino_alphabet)
        proteins._allocate(self._length)
        for i in range(self._length):
            assert self._refs != NULL
            assert self._refs[i] != NULL
            # get length of input and output sequences
            ntlen = self._refs[i].n
            if ntlen % 3 != 0:
                raise ValueError(f"Incomplete sequence of length {ntlen!r} at index {i!r}")
            aalen = ntlen // 3
            # create new object
            protein = DigitalSequence(
                genetic_code.amino_alphabet,
                name=self._storage[i].name,
                description=self._storage[i].description,
                accession=self._storage[i].accession,
                source=self._storage[i].source,
            )

            proteins._append(protein)
            # grow the internal sequence buffer
            status = libeasel.sq.esl_sq_GrowTo(protein._sq, aalen)
            if status == libeasel.eslEMEM:
                raise AllocationError("ESL_DSQ", sizeof(ESL_DSQ), aalen + 2)
            elif status != libeasel.eslOK:
                raise UnexpectedError(status, "esl_sq_Grow")
            # record sequence coordinates
            protein._sq.start = 1
            protein._sq.C = 0
            protein._sq.end = protein._sq.W = protein._sq.L = protein._sq.n = aalen

        # translate the sequences
        with nogil:
            for i in range(self._length):
                ntlen = self._refs[i].n
                aalen = proteins._refs[i].n
                genetic_code._translate(&self._refs[i].dsq[1], ntlen, &proteins._refs[i].dsq[1], aalen)
                proteins._refs[i].dsq[0] = proteins._refs[i].dsq[aalen+1] = libeasel.eslDSQ_SENTINEL

        proteins._largest = self._largest
        return proteins

    cpdef DigitalSequence largest(self):
        return SequenceBlock.largest(self)

    cpdef DigitalSequenceBlock copy(self):
        """Return a copy of the digital sequence block.

        Note:
            The sequence internally refered to by this collection are not
            copied. Use `copy.deepcopy` is you also want to duplicate the
            internal storage of each sequence.

        """
        cdef DigitalSequenceBlock new = DigitalSequenceBlock.__new__(DigitalSequenceBlock, self.alphabet)
        new._storage = self._storage.copy()
        new._length = self._length
        new._largest = self._largest
        new._owner = self._owner
        new._allocate(self._length)
        memcpy(new._refs, self._refs, self._length * sizeof(ESL_SQ*))
        return new


# --- Sequence File ----------------------------------------------------------

cdef class SequenceFile:
    """A wrapper around a sequence file, containing unaligned sequences.

    This class supports reading sequences stored in different formats, such
    as FASTA, GenBank or EMBL. The format of each file can be automatically
    detected, but it is also possible to pass an explicit format specifier
    when the `SequenceFile` is instantiated.

    Hint:
        `~pyhmmer.easel.SequenceFile` objects can also be used to parse
        files containing multiple sequence alignments: in that case, the
        sequences will be read sequentially, removing the gap characters::

            >>> with SequenceFile("tests/data/msa/LuxC.sto") as sf:
            ...     sequences = sf.read_block()
            >>> print(sequences[0].name[:6], sequences[0].sequence[:30])
            b'Q9KV99' LANQPLEAILGLINEARKSWSSTPELDPYR
            >>> print(sequences[1].name[:6], sequences[1].sequence[:30])
            b'Q2WLE3' IYSYPSEAMIEIINEYSKILCSDRKFLSYE

    .. versionadded:: 0.2.0
       The ``alphabet`` attribute.

    .. versionchanged:: 0.4.8
       Support reading sequences from a file-like handle.
       Support reading individual sequences from an MSA file.

    """

    _FORMATS = dict(SEQUENCE_FILE_FORMATS) # copy to prevent editing
    _EMPTY_FORMAT = "fasta" # fallback format used for empty files
    _EMPTY_ALPHABET = Alphabet.amino() # fallback alphabet used for empty files

    # --- Class methods ------------------------------------------------------

    @classmethod
    def parse(
        cls,
        const unsigned char[::1] buffer,
        str format,
        *,
        Alphabet alphabet=None
    ):
        """Parse a sequence from a binary ``buffer`` using the given ``format``.

        Argument:
            buffer (`bytes` or byte-like buffer): A buffer containing the
                sequence data to parse. Any type implementing the buffer
                protocol (such as `bytes`, `bytearray`, or `memoryview`)
                is supported.
            format (`str`): The format of the sequence data. See the
                `SequenceFile.__init__` documentation for allowed values.

        Keyword Arguments:
            alphabet (`~pyhmmer.easel.Alphabet`): The alphabet to use to
                digitize the returned sequence, if desired.

        Returns:
            `~pyhmmer.easel.Sequence`: The sequenced parsed from the buffer,
            either as a `DigitalSequence` if an alphabet was provided, or as
            a `TextSequence` if `None` was given.

        Raises:
            `ValueError`: When ``format`` is not a valid sequence format.
            `OSError`: If an internal parser error occurred while guessing
                the alphabet or the format.

        """
        cdef Sequence seq

        if alphabet is None:
            seq = TextSequence.__new__(TextSequence)
            seq._sq = libeasel.sq.esl_sq_Create()
            seq._sq.abc = NULL
        else:
            seq = DigitalSequence.__new__(DigitalSequence)
            seq.alphabet = alphabet
            seq._sq = libeasel.sq.esl_sq_CreateDigital(alphabet._abc)
        if not seq._sq:
            raise AllocationError("ESL_SQ", sizeof(ESL_SQ))

        return cls.parseinto(seq, buffer, format)

    @classmethod
    def parseinto(
        cls,
        Sequence seq,
        const unsigned char[::1] buffer,
        str format
    ):
        """Parse a sequence from a binary ``buffer`` into ``seq``.

        Argument:
            seq (`~pyhmmer.easel.Sequence`): The sequence object into which
                the deseriazlied sequence data will be written.
            buffer (`bytes` or byte-like buffer): A buffer containing the
                sequence data to parse. Any type implementing the buffer
                protocol (such as `bytes`, `bytearray`, or `memoryview`)
                is supported.
            format (`str`): The format of the sequence data. See the
                `SequenceFile.__init__` documentation for allowed values.

        Raises:
            `ValueError`: When ``format`` is not a valid sequence format.
            `OSError`: If an internal parser error occurred while guessing
                the alphabet or the format.

        Returns:
            `~pyhmmer.easel.Sequence`: The sequence given as argument, or
            `None` if the end of the file was reached.

        """
        assert seq._sq != NULL

        cdef int fmt = libeasel.sqio.eslSQFILE_UNKNOWN
        if format is not None:
            format_ = format.lower()
            if format_ not in SEQUENCE_FILE_FORMATS:
                raise InvalidParameter("format", format, choices=list(SEQUENCE_FILE_FORMATS))
            fmt = SEQUENCE_FILE_FORMATS[format_]

        cdef int status = libeasel.sqio.esl_sqio_Parse(
            <char*> &buffer[0],
            buffer.shape[0],
            seq._sq, fmt
        )
        if status == libeasel.eslEFORMAT:
            raise AllocationError("ESL_SQFILE", sizeof(ESL_SQFILE))
        elif status == libeasel.eslOK:
            return seq
        else:
            raise UnexpectedError(status, "esl_sqio_Parse")

    # --- Constructor --------------------------------------------------------

    @staticmethod
    cdef ESL_SQFILE* _open_fileobj(object fh, int fmt) except NULL:
        """Get an ``ESL_SQFILE*`` to read sequences from the file-like object ``fh``.

        Adapted from the ``esl_sqfile_Open`` function in ``esl_sqio.c`` and
        ``esl_sqascii_Open`` in ``esl_sqio_ascii.c`` to support using a
        file-like object instead of a filename only.

        Raises:
            `~pyhmmer.errors.AllocationError`: When either the ``ESL_SQFILE``
                or one of the string allocations fails.
            `NotImplementedError`: When calling `HMMFile._open_fileobj` with
                ``eslSQFILE_NCBI`` as the sequence format.

        """
        cdef int               n
        cdef int               status
        cdef ESL_SQFILE*       sqfp    = NULL
        cdef ESL_SQASCII_DATA* ascii   = NULL
        cdef FILE*             fp      = fopen_obj(fh, "r")
        cdef bytes             fh_repr = repr(fh).encode("ascii")

        # bail out early if format is not supported
        if fmt == libeasel.sqio.eslSQFILE_NCBI:
            fclose(fp)
            raise NotImplementedError("Cannot use a file-like object to read sequences from an NCBI database")

        # attempt to allocate space for the ESL_SQFILE
        sqfp = <ESL_SQFILE*> malloc(sizeof(ESL_SQFILE))
        if sqfp == NULL:
            fclose(fp)
            raise AllocationError("ESL_SQFILE", sizeof(ESL_SQFILE))

        # reset options
        sqfp.filename   = NULL
        sqfp.do_digital = False
        sqfp.abc        = NULL
        sqfp.format     = fmt

        # initialize function pointers
        sqfp.position       = &sqascii_Position
        sqfp.close          = &sqascii_Close
        sqfp.set_digital    = &sqascii_SetDigital
        sqfp.guess_alphabet = &sqascii_GuessAlphabet
        sqfp.is_rewindable  = &sqascii_IsRewindable
        sqfp.read           = &sqascii_Read
        sqfp.read_info      = &sqascii_ReadInfo
        sqfp.read_seq       = &sqascii_ReadSequence
        sqfp.read_window    = &sqascii_ReadWindow
        sqfp.echo           = &sqascii_Echo
        sqfp.read_block     = &sqascii_ReadBlock
        sqfp.open_ssi       = &sqascii_OpenSSI
        sqfp.pos_by_key     = &sqascii_PositionByKey
        sqfp.pos_by_number  = &sqascii_PositionByNumber
        sqfp.fetch          = &sqascii_Fetch
        sqfp.fetch_info     = &sqascii_FetchInfo
        sqfp.fetch_subseq   = &sqascii_FetchSubseq
        sqfp.get_error      = &sqascii_GetError

        # we don't need to look through the environment folder here, and
        # we know we are going to read from an ASCII file, so we can just
        # inline the `esl_sqascii_Open` function.
        ascii                  = &sqfp.data.ascii
        ascii.fp               = fp
        ascii.do_gzip          = False
        ascii.do_stdin         = False
        ascii.do_buffer        = False
        ascii.mem              = NULL
        ascii.allocm           = 0
        ascii.mn               = 0
        ascii.mpos             = 0
        ascii.moff             = -1
        ascii.is_recording     = False
        ascii.buf              = NULL
        ascii.boff             = 0
        ascii.balloc           = 0
        ascii.nc               = 0
        ascii.bpos             = 0
        ascii.L                = 0
        ascii.linenumber       = 0
        ascii.bookmark_offset  = 0
        ascii.bookmark_linenum = 0
        ascii.is_linebased     = False
        ascii.eof_is_ok        = False
        ascii.parse_header     = NULL
        ascii.skip_header      = NULL
        ascii.parse_end        = NULL
        ascii.afp              = NULL
        ascii.msa              = NULL
        ascii.idx              = -1
        ascii.ssifile          = NULL
        ascii.rpl              = -1
        ascii.bpl              = -1
        ascii.prvrpl           = -1
        ascii.prvbpl           = -1
        ascii.currpl           = -1
        ascii.curbpl           = -1
        ascii.ssi              = NULL

        try:
            # use the repr string of the file-like object as a filename
            sqfp.filename = strdup(<const char*> fh_repr)
            if sqfp.filename == NULL:
                raise AllocationError("char", sizeof(char), len(fh_repr))

            # if we don't know the format yet, try to autodetect
            if fmt == libeasel.sqio.eslSQFILE_UNKNOWN:
                status = sqascii_GuessFileFormat(sqfp, &fmt)
                if status == libeasel.eslOK:
                    sqfp.format = fmt
                elif status != libeasel.eslEFORMAT:
                    raise UnexpectedError(status, "sqascii_GuessFileFormat")

            # if we still couldn't guess, it may be an MSA, try opening it as such
            if fmt == libeasel.sqio.eslSQFILE_UNKNOWN or libeasel.sqio.esl_sqio_IsAlignment(fmt):
                ascii.afp = MSAFile._open_fileobj(fh, fmt)
                sqfp.format = fmt = ascii.afp.format

            # NOTE: at this point, we either successfully determined the format
            #       as a sequential format, or called MSAFile._open_fileobj,
            #       which raises an exception when it cannot determine the format
            #       either: in the following code, the format is necessarily
            #       known.

            # configure the parser and inmaps for this format
            if not libeasel.sqio.esl_sqio_IsAlignment(fmt):
                if fmt == libeasel.sqio.eslSQFILE_EMBL or fmt == libeasel.sqio.eslSQFILE_UNIPROT:
                    config_embl(sqfp)
                    inmap_embl(sqfp, NULL)
                elif fmt == libeasel.sqio.eslSQFILE_GENBANK or fmt == libeasel.sqio.eslSQFILE_DDBJ:
                    config_genbank(sqfp)
                    inmap_genbank(sqfp, NULL)
                elif fmt == libeasel.sqio.eslSQFILE_FASTA or fmt == libeasel.sqio.eslSQFILE_HMMPGMD:
                    config_fasta(sqfp)
                    inmap_fasta(sqfp, NULL)
                elif fmt == libeasel.sqio.eslSQFILE_DAEMON:
                    config_daemon(sqfp)
                    inmap_daemon(sqfp,  NULL)
                else:
                    raise ValueError("Unknown format code for `SequenceFile._open_fileobj`: {}".format(fmt))

                # preload the first line
                status = loadbuf(sqfp)
                if status == libeasel.eslEOF:
                    raise EOFError("Sequence file is empty")
                elif status != libeasel.eslOK:
                    raise UnexpectedError(status, "loadbuf")

                # skip the first line of HMMPGMD files, which is a header
                if fmt == libeasel.sqio.eslSQFILE_HMMPGMD:
                    status = fileheader_hmmpgmd(sqfp)
                    if status != libeasel.eslOK:
                        raise UnexpectedError(status, "fileheader_hmmpgmd")
            else:
                ascii.is_linebased = True
                ascii.eof_is_ok    = False
                ascii.parse_header = NULL
                ascii.skip_header  = NULL
                ascii.parse_end    = NULL

            # return the newly created sequence reader
            return sqfp
        except:
            # on error, make sure to clean up resources
            libeasel.sqio.esl_sqfile_Close(sqfp)
            raise

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self.alphabet = None
        self.name = None
        self._sqfp = NULL

    def __init__(
        self,
        object file,
        str format = None,
        *,
        bint digital = False,
        Alphabet alphabet = None,
    ):
        """__init__(self, file, format=None, *, digital=False, alphabet=None)\n--\n

        Create a new sequence file parser wrapping the given ``file``.

        Arguments:
            file (`str` or file-like object): Either the path to a file
                containing the sequences to read, or a file-like object
                opened in **binary mode**.
            format (`str`, optional): The format of the file, or `None` to
                autodetect. Supported values are: ``fasta``, ``embl``,
                ``genbank``, ``ddbj``, ``uniprot``, ``ncbi``, ``daemon``,
                ``hmmpgmd``, ``fmindex``, plus any format also supported
                by `~pyhmmer.easel.MSAFile`.
            digital (`bool`): Whether to read the sequences in text or digital
                mode. This will affect the type of `Sequence` objects returned
                later by the `read` function.
            alphabet (`~pyhmmer.easel.Alphabet`): The alphabet to use to
                digitize the sequences while reading.  If `None` given, it
                will be guessed based on the contents of the first sequence.

        Raises:
            `ValueError`: When ``format`` is not a valid sequence format.
            `OSError`: If an internal parser error occurred while guessing
                the alphabet or the format.

        Caution:
            `SequenceFile` can generally read sequences from binary-mode
            file-like objects, except for sequences in an NCBI BLAST
            database, since it is composed of multiple files. Reading from
            an NCBI BLAST database passed from a filename is however
            supported.

        .. versionchanged:: 0.4.4
           Added the ``ignore_gaps`` parameter.

        .. versionchanged:: 0.4.8
           Support reading from a file-like object (except NCBI format).

        .. versionchanged:: 0.5.0
           Added the ``digital`` and ``alphabet`` keyword arguments.

        .. deprecated:: 0.6.0
           The ``ignore_gaps`` keyword argument, use ``afa`` format instead.

        .. versionchanged:: 0.8.0
           Removed the ``ignore_gaps`` keyword argument.

        """
        cdef int   fmt
        cdef int   status
        cdef bytes fspath

        # get format from string passed as input
        fmt = libeasel.sqio.eslSQFILE_UNKNOWN
        if format is not None:
            format_ = format.lower()
            if format_ not in SEQUENCE_FILE_FORMATS:
                raise InvalidParameter("format", format, choices=list(SEQUENCE_FILE_FORMATS))
            fmt = SEQUENCE_FILE_FORMATS[format_]

        # open the given filename
        try:
            fspath = os.fsencode(file)
            # NOTE(@althonos): manually check the file is not a folder,
            # otherwise Easel will run into a segmentation fault after
            # failing to "slurp" the file!
            if os.path.isdir(fspath):
                raise IsADirectoryError(errno.EISDIR, f"Is a directory: {file!r}")
        except TypeError:
            self._sqfp = SequenceFile._open_fileobj(file, fmt)
            status = libeasel.eslOK
        else:
            status = libeasel.sqio.esl_sqfile_Open(fspath, fmt, NULL, &self._sqfp)
            self.name = os.fsdecode(fspath)

        # store a reference to the argument
        self._file = file

        # configure the file
        try:
            # check opening the file was successful
            if status == libeasel.eslENOTFOUND:
                raise FileNotFoundError(errno.ENOENT, "No such file or directory: {!r}".format(file))
            elif status == libeasel.eslEMEM:
                raise AllocationError("ESL_SQFILE", sizeof(ESL_SQFILE))
            elif status == libeasel.eslEFORMAT:
                if format is None:
                    raise ValueError("Could not determine format of file: {!r}".format(file))
                else:
                    raise EOFError("Sequence file is empty")
            elif status != libeasel.eslOK:
                raise UnexpectedError(status, "esl_sqfile_Open")
            # set digital mode if requested
            if digital:
                self.alphabet = self.guess_alphabet() if alphabet is None else alphabet
                if self.alphabet is None:
                    raise ValueError("Could not determine alphabet of file: {!r}".format(file))
                status = libeasel.sqio.esl_sqfile_SetDigital(self._sqfp, self.alphabet._abc)
                if status != libeasel.eslOK:
                    raise UnexpectedError(status, "esl_sqfile_SetDigital")
        except Exception as err:
            self.close()
            raise err

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def __dealloc__(self):
        if self._sqfp:
            warnings.warn("unclosed sequence file", ResourceWarning)
            self.close()

    def __iter__(self):
        return self

    def __next__(self):
        cdef Sequence seq = self.read()
        if seq is None:
            raise StopIteration()
        return seq

    def __repr__(self):
        cdef str  name = type(self).__name__
        cdef list args = [
            f"format={self.format!r}",
            f"digital={self.digital}",
            f"alphabet={self.alphabet!r}"
        ]
        if self.name is None:
            return f"<{name} file={self.file!r} {' '.join(args)}>"
        else:
            return f"{name}({self.name!r}, {', '.join(args)})"

    # --- Properties ---------------------------------------------------------

    @property
    def closed(self):
        """`bool`: Whether the `SequenceFile` is closed or not.
        """
        return self._sqfp == NULL

    @property
    def digital(self):
        """`bool`: Whether the `SequenceFile` is in digital mode or not.

        .. versionadded:: 0.5.0

        """
        return self.alphabet is not None

    @property
    def format(self):
        """`str`: The format of the `SequenceFile`.
        """
        if self._sqfp == NULL:
            raise ValueError("I/O operation on closed file.")
        return SEQUENCE_FILE_FORMATS_INDEX[self._sqfp.format]

    # --- Methods ------------------------------------------------------------

    cpdef void close(self) except *:
        """Close the file and free the resources used by the parser.
        """
        libeasel.sqio.esl_sqfile_Close(self._sqfp)
        self._sqfp = NULL

    cpdef Alphabet guess_alphabet(self):
        """Guess the alphabet of an open `SequenceFile`.

        This method tries to guess the alphabet of a sequence file by
        inspecting the first sequence in the file. It returns the alphabet,
        or `None` if the file alphabet cannot be reliably guessed.

        Raises:
            `EOFError`: if the file is empty.
            `OSError`: if a parse error occurred.
            `ValueError`: if this methods is called on a closed file.

        Example:
            >>> with SequenceFile("tests/data/seqs/bmyD.fna") as sf:
            ...     sf.guess_alphabet()
            Alphabet.dna()
            >>> with SequenceFile("tests/data/seqs/LuxC.faa") as sf:
            ...     sf.guess_alphabet()
            Alphabet.amino()

        .. versionadded:: 0.6.3

        """
        cdef int         ty
        cdef int         status
        cdef Alphabet    alphabet
        cdef const char* errbuf
        cdef str         msg

        if self._sqfp == NULL:
            raise ValueError("I/O operation on closed file.")

        status = libeasel.sqio.esl_sqfile_GuessAlphabet(self._sqfp, &ty)
        if status == libeasel.eslOK:
            alphabet = Alphabet.__new__(Alphabet)
            alphabet._init_default(ty)
            return alphabet
        elif status == libeasel.eslENOALPHABET or status == libeasel.eslEOD:
            return None
        elif status == libeasel.eslENODATA:
            raise EOFError("Sequence file appears to be empty.")
        elif status == libeasel.eslEFORMAT:
            errbuf = libeasel.sqio.esl_sqfile_GetErrorBuf(self._sqfp)
            msg = errbuf.decode("utf-8", "replace")
            raise ValueError("Could not parse file: {}".format(msg))
        else:
            raise UnexpectedError(status, "esl_sqfile_GuessAlphabet")

    cpdef Sequence read(self, bint skip_info=False, bint skip_sequence=False):
        """Read the next sequence from the file.

        Arguments:
            skip_info (`bool`): Pass `True` to disable reading the sequence
                *metadata*, and only read the sequence *letters*. Defaults to
                `False`.
            skip_sequence (`bool`): Pass `True` to disable reading the
                sequence *letters*, and only read the sequence *metadata*.
                Defaults to `False`.

        Returns:
            `Sequence`: The next sequence in the file, or `None` if all
            sequences were read from the file.

        Raises:
            `ValueError`: When attempting to read a sequence from a closed
                file, or when the file could not be parsed.

        Hint:
            This method allocates a new sequence, which is not efficient in
            case the sequences are being read within a tight loop. Use
            `SequenceFile.readinto` with an already initialized `Sequence`
            if you can to recycle the internal buffers.

        """
        cdef Sequence seq
        if self._sqfp == NULL:
            raise ValueError("I/O operation on closed file.")
        if self.alphabet is None:
            seq = TextSequence()
        else:
            seq = DigitalSequence(self.alphabet)
        return self.readinto(seq, skip_info=skip_info, skip_sequence=skip_sequence)

    cpdef Sequence readinto(self, Sequence seq, bint skip_info=False, bint skip_sequence=False):
        """Read the next sequence from the file, using ``seq`` to store data.

        Arguments:
            seq (`~pyhmmer.easel.Sequence`): A sequence object to use to store
                the next entry in the file. If this sequence was used before,
                it must be properly reset (using the `Sequence.clear` method)
                before using it again with `readinto`.
            skip_info (`bool`): Pass `True` to disable reading the sequence
                *metadata*, and only read the sequence *letters*. Defaults to
                `False`.
            skip_sequence (`bool`): Pass `True` to disable reading the
                sequence *letters*, and only read the sequence *metadata*.
                Defaults to `False`.

        Returns:
            `~pyhmmer.easel.Sequence`: A reference to ``seq`` that was passed
            as an input, or `None` if no sequences are left in the file.

        Raises:
            `ValueError`: When attempting to read a sequence from a closed
                file, or when the file could not be parsed.

        Example:
            Use `SequenceFile.readinto` to loop over the sequences in a file
            while recycling the same `Sequence` buffer:

            >>> with SequenceFile("tests/data/seqs/LuxC.faa") as sf:
            ...     seq = TextSequence()
            ...     while sf.readinto(seq) is not None:
            ...         # ... process seq here ... #
            ...         seq.clear()

        """
        assert seq._sq != NULL

        cdef str         funcname
        cdef int         status
        cdef const char* errbuf
        cdef str         msg

        if self._sqfp == NULL:
            raise ValueError("I/O operation on closed file.")

        if not skip_info and not skip_sequence:
            funcname = "esl_sqio_Read"
            status = libeasel.sqio.esl_sqio_Read(self._sqfp, seq._sq)
        elif not skip_info:
            funcname = "esl_sqio_ReadInfo"
            status = libeasel.sqio.esl_sqio_ReadInfo(self._sqfp, seq._sq)
        elif not skip_sequence:
            funcname = "esl_sqio_ReadSequence"
            status = libeasel.sqio.esl_sqio_ReadSequence(self._sqfp, seq._sq)
        else:
            raise ValueError("Cannot skip reading both sequence and metadata.")

        if status == libeasel.eslOK:
            return seq
        elif status == libeasel.eslEOF:
            return None
        elif status == libeasel.eslEFORMAT:
            errbuf = libeasel.sqio.esl_sqfile_GetErrorBuf(self._sqfp)
            msg = errbuf.decode("utf-8", "replace")
            raise ValueError("Could not parse file: {}".format(msg))
        else:
            raise UnexpectedError(status, funcname)

    cpdef SequenceBlock read_block(self, object sequences=None, object residues=None):
        """Read several sequences into a sequence block.

        Arguments:
            sequences (`int`, *optional*): The maximum number of sequences
                to read before returning a block. Leave as `None` to read all
                remaining sequences from the file.
            residues (`int`, *optional*): The number of residues to read
                before returning the block. Leave as `None` to keep reading
                sequences without a residue limit.

        Returns:
            `~pyhmmer.easel.SequenceBlock`: A sequence block object, which may
            be empty if there are no sequences to read anymore. The concrete
            type depends on whether the `SequenceFile` was opened in text or
            digital mode.

        Raises:
            `ValueError`: When attempting to read a sequence from a closed
                file, or when the file could not be parsed.

        Example:
            Read a block of at most 4 sequences from a sequence file::

                >>> with SequenceFile("tests/data/seqs/LuxC.faa") as sf:
                ...     block = sf.read_block(sequences=4)
                >>> len(block)
                4

            Read sequences until the block contains at least 1000 residues:

                >>> with SequenceFile("tests/data/seqs/LuxC.faa") as sf:
                ...     block = sf.read_block(residues=1000)
                >>> len(block)
                3
                >>> len(block[0]) + len(block[1]) + len(block[2])
                1444

            Note that the last sequence will not be truncated, so the block
            will always contain more than ``max_residues`` unless the end of
            the file was reached.

        .. versionadded:: 0.7.0

        """
        cdef SequenceBlock block
        cdef object        iterator
        cdef size_t        n_residues    = 0
        cdef size_t        n_sequences   = 0
        cdef size_t        max_sequences = SIZE_MAX if sequences is None else sequences
        cdef size_t        max_residues  = SIZE_MAX if residues is None else residues

        if self._sqfp == NULL:
            raise ValueError("I/O operation on closed file.")

        if self.alphabet is None:
            block = TextSequenceBlock()
        else:
            block = DigitalSequenceBlock(self.alphabet)
        if sequences is not None:
            block._allocate(sequences)

        while n_sequences < max_sequences and n_residues < max_residues:
            sequence = self.read()
            if sequence is None:
                break
            else:
                n_residues += len(sequence)
                n_sequences += 1
                block.append(sequence)

        return block

    cpdef void rewind(self) except *:
        """Rewind the file back to the beginning.

        For sequential formats, this method is supported for both *path*-based
        and *file object*-based sequence files. For multiple-sequence
        alignment formats, the underlying `MSAFile` needs to be reopened,
        so this is only supported for *path*-based files.

        Raises:
            `io.UnsupportedOperation`: When attempting to rewind a sequence
                file where the underlying stream is a file-like object that
                does not support the `~io.IOBase.seek` method.

        """
        if self._sqfp == NULL:
            raise ValueError("I/O operation on closed file")
        status = libeasel.sqio.esl_sqfile_Position(self._sqfp, 0)
        if status == libeasel.eslEMEM:
            raise AllocationError("ESL_SQFILE", sizeof(ESL_SQFILE))
        elif status == libeasel.eslESYS or status == libeasel.eslEINVAL:
            raise io.UnsupportedOperation("underlying stream is not seekable")
        elif status == libeasel.eslENOTFOUND:
            raise io.UnsupportedOperation("failed to re-open file")
        elif status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_sqfile_Position")

# --- Sequence/Subsequence Index ---------------------------------------------

cdef class SSIReader:
    """A read-only handler for sequence/subsequence index file.
    """

    Entry = collections.namedtuple(
        "Entry",
        ["fd", "record_offset", "data_offset", "record_length"]
    )

    FileInfo = collections.namedtuple(
        "FileInfo",
        ["name", "format"]
    )

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._ssi = NULL

    def __init__(self, object file):
        """__init__(self, file)\n--\n

        Create a new SSI file reader for the file at the given location.

        Arguments:
            file (`str`, `bytes` or `os.PathLike`): The path to a
                sequence/subsequence index file to read.

        """
        cdef int      status
        cdef bytes    fspath = os.fsencode(file)

        status = libeasel.ssi.esl_ssi_Open(fspath, &self._ssi)
        if status == libeasel.eslENOTFOUND:
            raise FileNotFoundError(2, "No such file or directory: {!r}".format(file))
        elif status == libeasel.eslEFORMAT:
            raise ValueError("File is not in correct SSI format")
        elif status == libeasel.eslERANGE:
            raise RuntimeError("File has 64-bit file offsets, which are unsupported on this system")
        elif status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_ssi_Open")

    def __dealloc__(self):
        libeasel.ssi.esl_ssi_Close(self._ssi)

    def __enter__(self):
        return self

    def __exit__(self, exc_value, exc_type, traceback):
        self.close()
        return False

    # --- Methods ------------------------------------------------------------

    def file_info(self, uint16_t fd):
        """Retrieve the `~pyhmmer.easel.SSIReader.FileInfo` of the descriptor.
        """
        cdef int   status
        cdef char* filename
        cdef int   format

        if self._ssi == NULL:
            raise ValueError("I/O operation on closed file.")
        if fd >= self._ssi.nfiles:
            raise IndexError(fd)

        status = libeasel.ssi.esl_ssi_FileInfo(self._ssi, fd, &filename, &format)
        if status == libeasel.eslOK:
            return self.FileInfo(os.fsdecode(filename), format)
        else:
            raise UnexpectedError(status, "esl_ssi_FileInfo")

    def find_name(self, bytes key):
        """Retrieve the `~pyhmmer.easel.SSIReader.Entry` for the given name.
        """
        cdef uint16_t ret_fh
        cdef off_t    ret_roff
        cdef off_t    opt_doff
        cdef int64_t  opt_L
        cdef int      status

        if self._ssi == NULL:
            raise ValueError("I/O operation on closed file.")

        status = libeasel.ssi.esl_ssi_FindName(
            self._ssi, key, &ret_fh, &ret_roff, &opt_doff, &opt_L
        )

        if status == libeasel.eslOK:
            return self.Entry(ret_fh, ret_roff, opt_doff or None, opt_L or None)
        elif status == libeasel.eslENOTFOUND:
            raise KeyError(key)
        elif status == libeasel.eslEFORMAT:
            raise ValueError("malformed index")
        else:
            raise UnexpectedError(status, "esl_ssi_FindName")

    cpdef void close(self):
        """Close the SSI file reader.
        """
        libeasel.ssi.esl_ssi_Close(self._ssi)
        self._ssi = NULL


cdef class SSIWriter:
    """A writer for sequence/subsequence index files.
    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._newssi = NULL

    def __init__(self, object file, bint exclusive = False):
        """__init__(self, file, exclusive=False)\n--\n

        Create a new SSI file write for the file at the given location.

        Arguments:
            file (`str`, `bytes` or `os.PathLike`): The path to a
                sequence/subsequence index file to write.
            exclusive (`bool`): Whether or not to create a file if one does
                not exist.

        Raises:
            `FileNotFoundError`: When the path to the file cannot be resolved.
            `FileExistsError`: When the file exists and ``exclusive`` is `True`.

        """
        cdef int   status
        cdef bytes fspath = os.fsencode(file)

        status = libeasel.ssi.esl_newssi_Open(fspath, not exclusive, &self._newssi)
        if status == libeasel.eslENOTFOUND:
            raise FileNotFoundError(2, "No such file or directory: {!r}".format(file))
        elif status == libeasel.eslEOVERWRITE:
            raise FileExistsError(17, "File exists: {!r}".format(file))
        elif status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_newssi_Open")

    def __dealloc__(self):
        if self._newssi != NULL:
            warnings.warn("unclosed SSI file", ResourceWarning)
            self.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    # --- Properties -----------------------------------------------------------

    @property
    def closed(self):
        return self._newssi == NULL

    # --- Utils --------------------------------------------------------------

    cdef void _on_write(self):
        if self.closed:
            raise ValueError("I/O operation on closed file.")

    # --- Methods ------------------------------------------------------------

    cpdef uint16_t add_file(self, object filename, int format = 0) except *:
        """Add a new file to the index.

        Arguments:
            filename (`str`, `bytes` or `os.PathLike`): The name of the
                file to register.
            format (`int`): A format code to associate with the file,
                or *0* by default.

        Returns:
            `int`: The filehandle associated with the new indexed file.

        """
        cdef int        status
        cdef uint16_t   fd

        self._on_write()
        name = os.fsencode(filename)
        status = libeasel.ssi.esl_newssi_AddFile(self._newssi, name, format, &fd)
        if status == libeasel.eslOK:
            return fd
        elif status == libeasel.eslERANGE:
            raise ValueError("Too many files registered in index.")
        else:
            _reraise_error()
            raise UnexpectedError(status, "esl_newssi_AddFile")

    cpdef void add_key(
        self,
        bytes key,
        uint16_t fd,
        off_t record_offset,
        off_t data_offset = 0,
        int64_t record_length = 0
    ) except *:
        """Add a new entry to the index with the given ``key``.
        """
        cdef int status

        self._on_write()
        status = libeasel.ssi.esl_newssi_AddKey(
            self._newssi, key, fd, record_offset, data_offset, record_length
        )

        if status == libeasel.eslERANGE:
            raise ValueError("Too many primary keys registered in index.")
        elif status == libeasel.eslENOTFOUND:
            raise OSError("Could not open external temporary files.")
        elif status != libeasel.eslOK:
            _reraise_error()
            raise UnexpectedError(status, "esl_newssi_AddKey")

    cpdef void add_alias(self, bytes alias, bytes key) except *:
        """Make ``alias`` an alias of ``key`` in the index.
        """
        cdef int status

        self._on_write()
        status = libeasel.ssi.esl_newssi_AddAlias(self._newssi, alias, key)
        if status == libeasel.eslOK:
            return
        elif status == libeasel.eslERANGE:
            raise ValueError("Too many secondary keys registed in index")
        elif status == libeasel.eslENOTFOUND:
            raise OSError("Could not open external temporary files.")
        else:
            _reraise_error()
            raise UnexpectedError(status, "esl_newssi_AddAlias")

    cpdef void close(self) except *:
        """Close the SSI file writer.
        """
        cdef int status

        if not self.closed:

            status = libeasel.ssi.esl_newssi_Write(self._newssi)
            if status == libeasel.eslERANGE:
                raise ValueError("Too many keys registered in index.")
            elif status == libeasel.eslESYS:
                raise RuntimeError("Extern sorting of keys failed.")
            elif status == libeasel.eslEDUP:
                raise ValueError("Index contains duplicate keys.")

            libeasel.ssi.esl_newssi_Close(self._newssi)
            self._newssi = NULL
