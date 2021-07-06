# coding: utf-8
# cython: language_level=3, linetrace=True
"""High-level interface to the Easel C library.

Easel is a library developed by the `Eddy/Rivas Lab <http://eddylab.org/>`_
to facilitate the development of biological software in C. It is used by
`HMMER <http://hmmer.org/>`_ and `Infernal <http://eddylab.org/infernal/>`_.

"""

# --- C imports --------------------------------------------------------------

cimport cython
from cpython cimport Py_buffer
from cpython.buffer cimport PyBUF_READ, PyBUF_WRITE, PyBUF_FORMAT, PyBUF_F_CONTIGUOUS
from cpython.unicode cimport PyUnicode_FromUnicode
from libc.stdint cimport int64_t, uint8_t, uint16_t, uint32_t, uint64_t
from libc.stdio cimport fclose
from libc.stdlib cimport calloc, malloc, free
from libc.string cimport memcpy, memmove, strdup, strlen, strncpy
from posix.types cimport off_t

cimport libeasel
cimport libeasel.alphabet
cimport libeasel.bitfield
cimport libeasel.keyhash
cimport libeasel.matrixops
cimport libeasel.msa
cimport libeasel.msafile
cimport libeasel.random
cimport libeasel.sq
cimport libeasel.sqio
cimport libeasel.sqio.ascii
cimport libeasel.ssi
cimport libeasel.vec
from libeasel cimport ESL_DSQ, esl_pos_t
from libeasel.sq cimport ESL_SQ
from libeasel.random cimport ESL_RANDOMNESS

DEF eslERRBUFSIZE = 128

IF UNAME_SYSNAME == "Linux":
    include "fileobj/linux.pxi"
ELIF UNAME_SYSNAME == "Darwin" or UNAME_SYSNAME.endswith("BSD"):
    include "fileobj/bsd.pxi"

# --- Python imports ---------------------------------------------------------

import abc
import os
import collections
import sys
import warnings

from .errors import AllocationError, UnexpectedError, AlphabetMismatch


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

    cdef void _init_default(self, int ty):
        with nogil:
            self._abc = libeasel.alphabet.esl_alphabet_Create(ty)
        if not self._abc:
            raise AllocationError("ESL_ALPHABET")

    @classmethod
    def amino(cls):
        """amino(cls)\n--

        Create a default amino-acid alphabet.

        """
        cdef Alphabet alphabet = Alphabet.__new__(Alphabet)
        alphabet._init_default(libeasel.alphabet.eslAMINO)
        return alphabet

    @classmethod
    def dna(cls):
        """dna(cls)\n--

        Create a default DNA alphabet.

        """
        cdef Alphabet alphabet = Alphabet.__new__(Alphabet)
        alphabet._init_default(libeasel.alphabet.eslDNA)
        return alphabet

    @classmethod
    def rna(cls):
        """rna(cls)\n--

        Create a default RNA alphabet.

        """
        cdef Alphabet alphabet = Alphabet.__new__(Alphabet)
        alphabet._init_default(libeasel.alphabet.eslRNA)
        return alphabet

    # def __init__(self, str alphabet, int K, int Kp):
    #     buffer = alphabet.encode('ascii')
    #     self._alphabet = libeasel.alphabet.esl_alphabet_CreateCustom(<char*> buffer, K, Kp)
    #     if not self._alphabet:
    #         raise AllocationError("ESL_ALPHABET")

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._abc = NULL

    def __dealloc__(self):
        libeasel.alphabet.esl_alphabet_Destroy(self._abc)

    def __repr__(self):
        assert self._abc != NULL
        if self._abc.type == libeasel.alphabet.eslRNA:
            return "Alphabet.rna()"
        elif self._abc.type == libeasel.alphabet.eslDNA:
            return "Alphabet.dna()"
        elif self._abc.type == libeasel.alphabet.eslAMINO:
            return "Alphabet.amino()"
        else:
            return "Alphabet({!r}, K={!r}, Kp={!r})".format(
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

    def __getstate__(self):
        assert self._abc != NULL
        return {"type": self._abc.type}

    def __setstate__(self, state):
        self._init_default(state["type"])

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
        cdef bytes symbols = self._abc.sym[:self._abc.Kp]
        return symbols.decode("ascii")

    # --- Utils --------------------------------------------------------------

    cdef inline bint _eq(self, Alphabet other) nogil:
        return self._abc.type == other._abc.type


# --- Bitfield ---------------------------------------------------------------

cdef class Bitfield:
    """A statically sized sequence of booleans stored as a packed bitfield.

    A bitfield is instantiated with a fixed length, and all booleans are set
    to `False` by default::

        >>> bitfield = Bitfield(8)
        >>> len(bitfield)
        8
        >>> bitfield[0]
        False

    Use indexing to access and edit individual bits::

        >>> bitfield[0] = True
        >>> bitfield[0]
        True
        >>> bitfield[0] = False
        >>> bitfield[0]
        False

    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._b = NULL

    def __dealloc__(self):
        libeasel.bitfield.esl_bitfield_Destroy(self._b)

    def __init__(self, size_t length):
        """__init__(self, length)\n--

        Create a new bitfield with the given ``length``.

        """
        with nogil:
            self._b = libeasel.bitfield.esl_bitfield_Create(length)
        if not self._b:
            raise AllocationError("ESL_BITFIELD")

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

    def __sizeof__(self):
        assert self._b != NULL
        cdef size_t nu = (self._b.nb // 64) + (self._b.nb % 64 != 0)
        return sizeof(uint64_t) * nu + sizeof(ESL_BITFIELD) + sizeof(self)

    # --- Utils --------------------------------------------------------------

    cdef size_t _wrap_index(self, int index) except -1:
        if index < 0:
            index += self._b.nb
        if index >= self._b.nb or index < 0:
            raise IndexError("bitfield index out of range")
        return <size_t> index

    # --- Methods ------------------------------------------------------------

    cpdef size_t count(self, bint value=1):
        """count(self, value=True)\n--

        Count the number occurrences of ``value`` in the bitfield.

        If no argument is given, counts the number of `True` occurences.

        Example:
            >>> bitfield = Bitfield(8)
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

    cpdef void toggle(self, int index):
        """toggle(self, index)\n--

        Switch the value of one single bit.

        Example:
            >>> bitfield = Bitfield(8)
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
        """__init__(self)\n--

        Create a new empty key-hash collection.

        """
        with nogil:
            self._kh = libeasel.keyhash.esl_keyhash_Create()
        if not self._kh:
            raise AllocationError("ESL_KEYHASH")

    def __copy__(self):
        return self.copy()

    def __len__(self):
        assert self._kh != NULL
        return libeasel.keyhash.esl_keyhash_GetNumber(self._kh)

    def __contains__(self, object value):
        assert self._kh != NULL

        if not isinstance(value, bytes):
            return False

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

        cdef const char*  key    = item
        cdef       size_t length = len(item)
        cdef       int    index

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
        """add(self, item)\n--

        Add a new key to the hash table, and return its index.

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

    cpdef void clear(self):
        """clear(self)\n--

        Remove all entries from the collection.

        """
        cdef int status
        with nogil:
            status = libeasel.keyhash.esl_keyhash_Reuse(self._kh)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_keyhash_Reuse")

    cpdef KeyHash copy(self):
        """copy(self)\n--

        Create and return an exact copy of this mapping.

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
            raise AllocationError("ESL_KEYHASH")
        return new


# --- Matrix & Vector --------------------------------------------------------

cdef class Vector:
    """An abstract 1D array of fixed size.

    .. versionadded:: 0.4.0

    """

    # --- Class methods ------------------------------------------------------

    @classmethod
    def zeros(cls, int n):
        """Create a vector of size ``n`` filled with zeros.
        """
        raise TypeError("Can't instantiate abstract class 'Vector'")

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._owner = None
        self._n = 0
        self._shape = (0,)

    def __init__(self, object iterable = None):
        raise TypeError("Can't instantiate abstract class 'Vector'")

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
        return tuple(self._strides)

    # --- Methods ------------------------------------------------------------

    def argmax(self):
        """Return index of the maximum element in the vector.
        """

    def argmin(self):
        """Return index of the minimum element in the vector.
        """

    def copy(self):
        """Create a copy of the vector, allocating a new buffer.
        """

    def max(self):
        """Return value of the maximum element in the vector.
        """

    def min(self):
        """Return value of the minimum element in the vector.
        """

    def reverse(self):
        """Reverse the vector, in place.
        """

    def sum(self):
        """Returns the scalar sum of all elements in the vector.

        Floating point summations use Kahan compensated summation, in order
        to minimize roundoff error accumulation.  Additionally, they are most
        accurate if the vector is sorted in increasing order, from small to
        large, so you may consider sorting the vector before summing it.

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

        >>> v = VectorF(range(10))
        >>> v[2:5]
        VectorF([2.0, 3.0, 4.0])
        >>> v[5:-1] = 10.0
        >>> v
        VectorF([0.0, 1.0, 2.0, 3.0, 4.0, 10.0, 10.0, 10.0, 10.0, 9.0])

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
        ValueError: cannot pairwise multiply vectors of different size

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

    # --- Class methods ------------------------------------------------------

    @classmethod
    def zeros(cls, n):
        if n <= 0:
            raise ValueError("Cannot create a vector with negative or null size")
        cdef VectorF vec = VectorF.__new__(VectorF)
        vec._n = vec._shape[0] = n
        vec._data = <float*> calloc(n, sizeof(float))

        if vec._data == NULL:
            raise AllocationError("float*")
        return vec

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._data = NULL
        self._strides = (sizeof(float),)

    def __dealloc__(self):
        if self._owner is None and self._data != NULL:
            free(self._data)
        self._data = NULL

    def __init__(self, object iterable):

        cdef size_t i
        cdef int    n = len(iterable)

        if n <= 0:
            raise ValueError("Cannot create a vector with negative or null size")

        self._n = self._shape[0] = n
        if self._data == NULL: # avoid realloc if __init__ called more than once
            self._data = <float*> malloc(n * sizeof(float))
            if self._data == NULL:
                raise AllocationError("float*")
        for i, item in enumerate(iterable):
            self._data[i] = item

    def __copy__(self):
        return self.copy()

    def __len__(self):
        return self._n

    def __getitem__(self, object index):
        assert self._data != NULL

        cdef VectorF new
        cdef int     idx

        if isinstance(index, slice):
            start, stop, step = index.indices(self._n)
            if step != 1:
                raise ValueError(f"cannot slice a Vector with step other than 1")
            new = VectorF.__new__(VectorF)
            new._owner = self
            new._n = new._shape[0] = stop - start
            new._data = &(self._data[start])
            return new
        else:
            idx = index
            if idx < 0:
                idx += self._n
            if idx < 0 or idx >= self._n:
                raise IndexError("vector index out of range")
            return self._data[idx]

    def __setitem__(self, object index, float value):
        assert self._data != NULL

        cdef int x

        if isinstance(index, slice):
            for x in range(*index.indices(self._n)):
                self._data[x] = value
        else:
            x = index
            if x < 0:
                x += self._n
            if x < 0  or x >= self._n:
                raise IndexError("vector index out of range")
            self._data[x] = value

    def __getbuffer__(self, Py_buffer* buffer, int flags):
        assert self._data != NULL

        if flags & PyBUF_FORMAT:
            buffer.format = "f"
        else:
            buffer.format = NULL
        buffer.buf = <void*> &(self._data[0])
        buffer.internal = NULL
        buffer.itemsize = sizeof(float)
        buffer.len = self._n * sizeof(float)
        buffer.ndim = 1
        buffer.obj = self
        buffer.readonly = 0
        buffer.shape = self._shape
        buffer.strides = self._strides
        buffer.suboffsets = NULL

    def __releasebuffer__(self, Py_buffer* buffer):
        pass

    def __add__(VectorF self, object other):
        assert self._data != NULL
        cdef VectorF new = self.copy()
        return new.__iadd__(other)

    def __iadd__(self, object other):
        assert self._data != NULL

        cdef VectorF other_vec
        cdef float   other_f

        if isinstance(other, VectorF):
            other_vec = other
            assert other_vec._data != NULL
            if self._n != other_vec._n:
                raise ValueError("cannot add vectors of different size")
            with nogil:
                libeasel.vec.esl_vec_FAdd(self._data, other_vec._data, self._n)
        else:
            other_f = other
            with nogil:
                libeasel.vec.esl_vec_FIncrement(self._data, self._n, other_f)
        return self

    def __mul__(VectorF self, object other):
        assert self._data != NULL
        cdef VectorF new = self.copy()
        return new.__imul__(other)

    def __imul__(self, object other):
        assert self._data != NULL

        cdef VectorF other_vec
        cdef float   other_f
        cdef int     i

        if isinstance(other, VectorF):
            other_vec = other
            if self._n != other_vec._n:
                raise ValueError("cannot pairwise multiply vectors of different size")
            assert other_vec._data != NULL
            # NB(@althonos): There is no function in `vectorops.h` to do this
            # for now...
            for i in range(self._n):
                self._data[i] *= other_vec._data[i]
        else:
            other_f = other
            with nogil:
                libeasel.vec.esl_vec_FScale(self._data, self._n, other_f)
        return self

    def __matmul__(VectorF self, object other):
        assert self._data != NULL

        cdef float* other_data
        cdef float* self_data  = self._data
        cdef int    n          = self._n
        cdef float  res

        if isinstance(other, VectorF):
            other_data = (<VectorF> other)._data
            assert other_data != NULL
            if len(self) != len(other):
                raise ValueError("cannot multiply vectors of different size")
            with nogil:
                res = libeasel.vec.esl_vec_FDot(self_data, other_data, n)
            return res

        return NotImplemented

    def __repr__(self):
        return f"VectorF({list(self)!r})"

    def __sizeof__(self):
        return self._n * sizeof(float) + sizeof(self)

    # --- Methods ------------------------------------------------------------

    cpdef int argmax(self):
        assert self._data != NULL
        with nogil:
            return libeasel.vec.esl_vec_FArgMax(self._data, self._n)

    cpdef int argmin(self):
        assert self._data != NULL
        with nogil:
            return libeasel.vec.esl_vec_FArgMin(self._data, self._n)

    cpdef VectorF copy(self):
        assert self._data != NULL

        cdef VectorF new = VectorF.__new__(VectorF)
        new._n = new._shape[0] = self._n
        new._data = <float*> malloc(self._n * sizeof(float))
        if new._data == NULL:
            raise AllocationError("float*")
        with nogil:
            libeasel.vec.esl_vec_FCopy(self._data, self._n, new._data)
        return new

    cpdef float max(self):
        assert self._data != NULL
        with nogil:
            return libeasel.vec.esl_vec_FMax(self._data, self._n)

    cpdef float min(self):
        assert self._data != NULL
        with nogil:
            return libeasel.vec.esl_vec_FMin(self._data, self._n)

    cpdef void normalize(self):
        r"""Normalize a vector so that all elements sum to 1.

        Caution:
            If sum is zero, sets all elements to :math:`\frac{1}{n}`,
            where :math:`n` is the size of the vector.

        """
        assert self._data != NULL
        with nogil:
            libeasel.vec.esl_vec_FNorm(self._data, self._n)

    cpdef void reverse(self):
        assert self._data != NULL
        with nogil:
            libeasel.vec.esl_vec_FReverse(self._data, self._data, self._n)

    cpdef float sum(self):
        assert self._data != NULL
        with nogil:
            return libeasel.vec.esl_vec_FSum(self._data, self._n)


cdef class VectorU8(Vector):
    """A vector storing byte-sized unsigned integers.

    .. versionadded:: v0.4.0

    """

    # --- Class methods ------------------------------------------------------

    @classmethod
    def zeros(cls, n):
        if n <= 0:
            raise ValueError("Cannot create a vector with negative or null size")
        cdef VectorU8 vec = VectorU8.__new__(VectorU8)
        vec._n = vec._shape[0] = n
        vec._data = <uint8_t*> calloc(n, sizeof(uint8_t))

        if vec._data == NULL:
            raise AllocationError("uint8_t*")
        return vec

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._data = NULL
        self._strides = (sizeof(uint8_t),)

    def __dealloc__(self):
        if self._owner is None and self._data != NULL:
            free(self._data)
        self._data = NULL

    def __init__(self, object iterable):

        cdef size_t i
        cdef int    n = len(iterable)

        if n <= 0:
            raise ValueError("Cannot create a vector with negative or null size")

        self._n = self._shape[0] = n
        if self._data == NULL: # avoid realloc if __init__ called more than once
            self._data = <uint8_t*> malloc(n * sizeof(uint8_t))
            if self._data == NULL:
                raise AllocationError("uint8_t*")
        for i, item in enumerate(iterable):
            self._data[i] = item

    def __copy__(self):
        return self.copy()

    def __len__(self):
        return self._n

    def __getitem__(self, object index):
        assert self._data != NULL

        cdef VectorU8 new
        cdef int      idx

        if isinstance(index, slice):
            start, stop, step = index.indices(self._n)
            if step != 1:
                raise ValueError(f"cannot slice a Vector with step other than 1")
            new = VectorU8.__new__(VectorU8)
            new._owner = self
            new._n = new._shape[0] = stop - start
            new._data = &(self._data[start])
            return new
        else:
            idx = index
            if idx < 0:
                idx += self._n
            if idx < 0 or idx >= self._n:
                raise IndexError("vector index out of range")
            return self._data[idx]

    def __setitem__(self, object index, uint8_t value):
        assert self._data != NULL

        cdef int x

        if isinstance(index, slice):
            for x in range(*index.indices(self._n)):
                self._data[x] = value
        else:
            x = index
            if x < 0:
                x += self._n
            if x < 0  or x >= self._n:
                raise IndexError("vector index out of range")
            self._data[x] = value

    def __getbuffer__(self, Py_buffer* buffer, int flags):
        assert self._data != NULL

        if flags & PyBUF_FORMAT:
            buffer.format = "B"
        else:
            buffer.format = NULL
        buffer.buf = <void*> &(self._data[0])
        buffer.internal = NULL
        buffer.itemsize = sizeof(uint8_t)
        buffer.len = self._n * sizeof(uint8_t)
        buffer.ndim = 1
        buffer.obj = self
        buffer.readonly = 0
        buffer.shape = self._shape
        buffer.strides = self._strides
        buffer.suboffsets = NULL

    def __releasebuffer__(self, Py_buffer* buffer):
        pass

    def __add__(VectorU8 self, object other):
        assert self._data != NULL
        cdef VectorU8 new = self.copy()
        return new.__iadd__(other)

    def __iadd__(self, object other):
        assert self._data != NULL

        cdef VectorU8 other_vec
        cdef uint8_t  other_f
        cdef int      i

        if isinstance(other, VectorU8):
            other_vec = other
            assert other_vec._data != NULL
            if self._n != other_vec._n:
                raise ValueError("cannot add vectors of different size")
            for i in range(self._n):
                self._data[i] += other_vec._data[i]
        else:
            other_f = other
            for i in range(self._n):
                self._data[i] += other_f
        return self

    def __mul__(VectorU8 self, object other):
        assert self._data != NULL
        cdef VectorU8 new = self.copy()
        return new.__imul__(other)

    def __imul__(self, object other):
        assert self._data != NULL

        cdef VectorU8 other_vec
        cdef uint8_t  other_f
        cdef int      i

        if isinstance(other, VectorU8):
            other_vec = other
            if self._n != other_vec._n:
                raise ValueError("cannot pairwise multiply vectors of different size")
            assert other_vec._data != NULL
            # NB(@althonos): There is no function in `vectorops.h` to do this
            # for now...
            for i in range(self._n):
                self._data[i] *= other_vec._data[i]
        else:
            other_f = other
            for i in range(self._n):
                self._data[i] *= other_f
        return self

    def __matmul__(VectorU8 self, object other):
        assert self._data != NULL

        cdef long int  res = 0
        cdef int       i

        if isinstance(other, VectorU8):
            other_data = (<VectorU8> other)._data
            assert other_data != NULL
            if len(self) != len(other):
                raise ValueError("cannot multiply vectors of different size")
            for i in range(self._n):
                res += self._data[i] * other_data[i]
            return res

        return NotImplemented

    def __repr__(self):
        return f"VectorU8({list(self)!r})"

    def __sizeof__(self):
        return self._n * sizeof(uint8_t) + sizeof(self)

    # --- Methods ------------------------------------------------------------

    cpdef int argmax(self):
        assert self._data != NULL

        cdef int i
        cdef int best = 0;

        with nogil:
            for i in range(1, self._n):
                if self._data[i] > self._data[best]:
                    best = i
        return best

    cpdef int argmin(self):
        assert self._data != NULL

        cdef int i
        cdef int best = 0;

        with nogil:
            for i in range(1, self._n):
                if self._data[i] < self._data[best]:
                    best = i
        return best

    cpdef VectorU8 copy(self):
        assert self._data != NULL

        cdef VectorU8 new = VectorU8.__new__(VectorU8)
        new._n = new._shape[0] = self._n
        new._data = <uint8_t*> malloc(self._n * sizeof(uint8_t))

        if new._data == NULL:
          raise AllocationError("uint8_t*")
        with nogil:
            memcpy(<void*> new._data, <void*> self._data, self._n * sizeof(uint8_t))

        return new

    cpdef uint8_t max(self):
        assert self._data != NULL

        cdef int     i
        cdef uint8_t best

        with nogil:
            best = self._data[0]
            for i in range(1, self._n):
                if self._data[i] > best:
                    best = self._data[i]
        return best

    cpdef uint8_t min(self):
        assert self._data != NULL

        cdef int     i
        cdef uint8_t best

        with nogil:
            best = self._data[0]
            for i in range(1, self._n):
                if self._data[i] < best:
                    best = self._data[i]
        return best

    cpdef void reverse(self):
        assert self._data != NULL

        cdef int     i
        cdef uint8_t x

        with nogil:
            for i in range(self._n // 2):
                x = self._data[self._n - i - 1]
                self._data[self._n - i - 1] = self._data[i]
                self._data[i] = x

    cpdef uint8_t sum(self):
        """Returns the scalar sum of all elements in the vector.

        Caution:
            The sum is wrapping::

                >>> vec = VectorU8([255, 2])
                >>> vec.sum()
                1

        """

        assert self._data != NULL

        cdef int     i
        cdef uint8_t sum = 0

        with nogil:
            for i in range(self._n):
                sum += self._data[i]
        return sum


cdef class Matrix:
    """An abstract 2D array of fixed size.

    .. versionadded:: 0.4.0

    """

    # --- Class methods ------------------------------------------------------

    @classmethod
    def zeros(cls, int m, int n):
        r"""Create a new :math:`m \times n` matrix filled with zeros.
        """
        raise TypeError("Can't instantiate abstract class 'Matrix'")

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._owner = None
        self._n = self._m = 0
        self._shape = (self._n, self._m)

    def __init__(self, object iterable = None):
      raise TypeError("Can't instantiate abstract class 'Matrix'")

    def __len__(self):
        return self._m

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
        return tuple(self._strides)

    # --- Methods ------------------------------------------------------------

    def argmax(self):
        """Return the coordinates of the maximum element in the matrix.
        """

    def argmin(self):
        """Return the coordinates of the minimum element in the matrix.
        """

    def copy(self):
        """Create a copy of the matrix, allocating a new buffer.
        """

    def max(self):
        """Return the value of the maximum element in the matrix.
        """

    def min(self):
        """Return the value of the minimum element in the matrix.
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

    # --- Class methods ------------------------------------------------------

    @classmethod
    def zeros(cls, int m, int n):
        if n <= 0:
            raise ValueError("Cannot create a matrix with negative or null dimension")
        cdef MatrixF mat = MatrixF.__new__(MatrixF)
        mat._m = mat._shape[0] = m
        mat._n = mat._shape[1] = n
        mat._strides = (sizeof(float), mat._m * sizeof(float))

        # allocate the array of pointers
        mat._data    = <float**> calloc(m, sizeof(float*))
        if mat._data == NULL:
            raise AllocationError("float**")

        # allocate the data array
        mat._data[0] = <float*> calloc(m*n, sizeof(float))
        if mat._data[0] == NULL:
            raise AllocationError("float*")

        # update the pointer offsets in the array of pointers
        cdef int i
        for i in range(1, m):
            mat._data[i] = mat._data[0] + i*n

        return mat

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._data = NULL
        self._strides = (sizeof(float), self._n * sizeof(float))

    def __dealloc__(self):
        if self._owner is None and self._data != NULL:
            free(self._data[0])
            free(self._data)
        self._data = NULL

    def __init__(self, object iterable):

        self._m = self._shape[0] = len(iterable)
        if self._m <= 0:
            raise ValueError("Cannot create a matrix with null number of columns")

        self._n = self._shape[1] = len(iterable[0])
        if self._n <= 0:
            raise ValueError("Cannot create a matrix with null number of rows")

        self._strides = (self._m * sizeof(float), sizeof(float))

        if self._data == NULL:
            self._data = libeasel.matrixops.esl_mat_FCreate(self._m, self._n)
            if self._data == NULL:
                raise AllocationError("float**")

        cdef size_t i
        cdef size_t j
        for i, row in enumerate(iterable):
            if len(row) != self._n:
                raise ValueError("Inconsistent number of rows in input")
            for j, val in enumerate(row):
                self._data[i][j] = val

    def __copy__(self):
        return self.copy()

    def __add__(MatrixF self, object other):
        assert self._data != NULL
        cdef MatrixF new = self.copy()
        return new.__iadd__(other)

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
                libeasel.vec.esl_vec_FAdd( self._data[0], other_mat._data[0], self._m * self._n )
        else:
            other_f = other
            with nogil:
                libeasel.vec.esl_vec_FIncrement(self._data[0], self._m*self._n, other_f)
        return self

    def __mul__(MatrixF self, object other):
        assert self._data != NULL
        cdef MatrixF new = self.copy()
        return new.__imul__(other)

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
                self._data[0][i] *= other_mat._data[0][i]
        else:
            other_f = other
            with nogil:
                libeasel.matrixops.esl_mat_FScale(self._data, self._m, self._n, other_f)
        return self

    def __getitem__(self, object index):
        assert self._data != NULL

        cdef int x
        cdef int y

        if isinstance(index, slice):
            return NotImplemented

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
            return self._data[x][y]

        x = index
        if x < 0:
            x += self._m
        if x < 0 or x >= self._m:
            raise IndexError("vector index out of range")

        cdef VectorF row = VectorF.__new__(VectorF)
        row._owner = self
        row._n = row._shape[0] = self._n
        row._data = self._data[x]
        return row

    def __setitem__(self, object index, float value):
        assert self._data != NULL

        cdef int x
        cdef int y

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
            self._data[x][y] = value

        else:
            raise TypeError("Matrix.__setitem__ can only be used with a 2D index")

    def __getbuffer__(self, Py_buffer* buffer, int flags):
        assert self._data != NULL

        if flags & PyBUF_FORMAT:
            buffer.format = "f"
        else:
            buffer.format = NULL
        buffer.buf = <void*> &(self._data[0][0])
        buffer.internal = NULL
        buffer.itemsize = sizeof(float)
        buffer.len = self._m * self._n * sizeof(float)
        buffer.ndim = 2
        buffer.obj = self
        buffer.readonly = 0
        buffer.shape = self._shape
        buffer.strides = self._strides
        buffer.suboffsets = NULL

    def __releasebuffer__(self, Py_buffer* buffer):
        pass

    def __repr__(self):
        return f"MatrixF({[list(row) for row in self]!r})"

    def __sizeof__(self):
        return (
            self._n * sizeof(float*)
          + self._m * self._n * sizeof(float)
          + sizeof(self)
        )

    # --- Methods ------------------------------------------------------------

    cpdef tuple argmax(self):
        assert self._data != NULL

        cdef int n
        cdef int x
        cdef int y

        with nogil:
            n = libeasel.vec.esl_vec_FArgMax(self._data[0], self._m*self._n)
            x = n // self._m
            y = n % self._n

        return x, y

    cpdef tuple argmin(self):
        assert self._data != NULL

        cdef int n
        cdef int x
        cdef int y

        with nogil:
            n = libeasel.vec.esl_vec_FArgMin(self._data[0], self._m*self._n)
            x = n // self._m
            y = n % self._n

        return x, y

    cpdef MatrixF copy(self):
        cdef MatrixF mat = MatrixF.__new__(MatrixF)
        mat._m = mat._shape[0] = self._m
        mat._n = mat._shape[1] = self._n
        mat._strides = self._strides
        with nogil:
            mat._data = libeasel.matrixops.esl_mat_FClone(self._data, self._m, self._n)
        if mat._data == NULL:
            raise AllocationError("float**")
        return mat

    cpdef float max(self):
        assert self._data != NULL
        with nogil:
            return libeasel.matrixops.esl_mat_FMax(self._data, self._m, self._n)

    cpdef float min(self):
        assert self._data != NULL
        with nogil:
            return libeasel.vec.esl_vec_FMin(self._data[0], self._m*self._n)

    cpdef float sum(self):
        assert self._data != NULL
        with nogil:
            return libeasel.vec.esl_vec_FSum(self._data[0], self._m*self._n)


cdef class MatrixU8(Matrix):
    """A matrix storing byte-sized unsigned integers.

    .. versionadded:: v0.4.0

    """

    # --- Class methods ------------------------------------------------------

    @classmethod
    def zeros(cls, int m, int n):
        if n <= 0:
            raise ValueError("Cannot create a matrix with negative or null dimension")
        cdef MatrixU8 mat = MatrixU8.__new__(MatrixU8)
        mat._m = mat._shape[0] = m
        mat._n = mat._shape[1] = n
        mat._strides = (sizeof(uint8_t), mat._m * sizeof(uint8_t))

        # allocate the array of pointers
        mat._data    = <uint8_t**> calloc(m, sizeof(uint8_t*))
        if mat._data == NULL:
            raise AllocationError("uint8_t**")

        # allocate the data array
        mat._data[0] = <uint8_t*> calloc(m*n, sizeof(uint8_t))
        if mat._data[0] == NULL:
            raise AllocationError("uint8_t*")

        # update the pointer offsets in the array of pointers
        cdef int i
        for i in range(1, m):
            mat._data[i] = mat._data[0] + i*n

        return mat

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._data = NULL
        self._strides = (sizeof(uint8_t), self._n * sizeof(uint8_t))

    def __dealloc__(self):
        if self._owner is None and self._data != NULL:
            free(self._data[0])
            free(self._data)
        self._data = NULL

    def __init__(self, object iterable):

        cdef int i
        cdef int j

        self._m = self._shape[0] = len(iterable)
        if self._m <= 0:
            raise ValueError("Cannot create a matrix with null number of columns")

        self._n = self._shape[1] = len(iterable[0])
        if self._n <= 0:
            raise ValueError("Cannot create a matrix with null number of rows")

        if self._data == NULL:
            # allocate the array of pointers
            self._data    = <uint8_t**> calloc(self._m, sizeof(uint8_t*))
            if self._data == NULL:
                raise AllocationError("uint8_t**")
            # allocate the data array
            self._data[0] = <uint8_t*> calloc(self._m*self._n, sizeof(uint8_t))
            if self._data[0] == NULL:
                raise AllocationError("uint8_t*")
            # update the pointer offsets in the array of pointers
            for i in range(1, self._m):
                self._data[i] = self._data[0] + i * self._n

        self._strides = (self._m * sizeof(uint8_t), sizeof(uint8_t))
        for i, row in enumerate(iterable):
            if len(row) != self._n:
                raise ValueError("Inconsistent number of rows in input")
            for j, val in enumerate(row):
                self._data[i][j] = val

    def __copy__(self):
        return self.copy()

    def __add__(MatrixU8 self, object other):
        assert self._data != NULL
        cdef MatrixU8 new = self.copy()
        return new.__iadd__(other)

    def __iadd__(self, object other):
        assert self._data != NULL

        cdef int      i
        cdef MatrixU8 other_mat
        cdef uint8_t  other_n

        if isinstance(other, MatrixU8):
            other_mat = other
            assert other_mat._data != NULL
            if other_mat._m != self._m or other_mat._n != self._n:
                raise ValueError(f"cannot pairwise add {other_mat.shape} matrix to {self.shape} matrix")
            with nogil:
                for i in range(self._m * self._n):
                    self._data[0][i] += other_mat._data[0][i]
        else:
            other_n = other
            with nogil:
                for i in range(self._m * self._n):
                    self._data[0][i] += other_n
        return self

    def __mul__(MatrixU8 self, object other):
        assert self._data != NULL
        cdef MatrixU8 new = self.copy()
        return new.__imul__(other)

    def __imul__(self, object other):
        assert self._data != NULL

        cdef MatrixU8 other_mat
        cdef uint8_t  other_n
        cdef int      i

        if isinstance(other, MatrixU8):
            other_mat = other
            assert other_mat._data != NULL
            if other_mat._m != self._m or other_mat._n != self._n:
                raise ValueError(f"cannot pairwise multiply {other_mat.shape} matrix with {self.shape} matrix")
            # NB(@althonos): There is no function in `vectorops.h` to do this
            # for now...
            for i in range(self._n * self._m):
                self._data[0][i] *= other_mat._data[0][i]
        else:
            other_n = other
            with nogil:
                for i in range(self._n * self._m):
                    self._data[0][i] *= other_n
        return self

    def __getitem__(self, object index):
        assert self._data != NULL

        cdef int x
        cdef int y

        if isinstance(index, slice):
            return NotImplemented

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
            return self._data[x][y]

        x = index
        if x < 0:
            x += self._m
        if x < 0 or x >= self._m:
            raise IndexError("vector index out of range")

        cdef VectorU8 row = VectorU8.__new__(VectorU8)
        row._owner = self
        row._n = row._shape[0] = self._n
        row._data = self._data[x]
        return row

    def __setitem__(self, object index, uint8_t value):
        assert self._data != NULL

        cdef int x
        cdef int y

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
            self._data[x][y] = value

        else:
            raise TypeError("Matrix.__setitem__ can only be used with a 2D index")

    def __getbuffer__(self, Py_buffer* buffer, int flags):
        assert self._data != NULL

        if flags & PyBUF_FORMAT:
            buffer.format = "B"
        else:
            buffer.format = NULL
        buffer.buf = <void*> &(self._data[0][0])
        buffer.internal = NULL
        buffer.itemsize = sizeof(uint8_t)
        buffer.len = self._m * self._n * sizeof(uint8_t)
        buffer.ndim = 2
        buffer.obj = self
        buffer.readonly = 0
        buffer.shape = self._shape
        buffer.strides = self._strides
        buffer.suboffsets = NULL

    def __releasebuffer__(self, Py_buffer* buffer):
        pass

    def __repr__(self):
        return f"MatrixU8({[list(row) for row in self]!r})"

    def __sizeof__(self):
        return (
            self._n * sizeof(uint8_t*)
          + self._m * self._n * sizeof(uint8_t)
          + sizeof(self)
        )

    # --- Methods ------------------------------------------------------------

    cpdef tuple argmax(self):
        assert self._data != NULL

        cdef int i
        cdef int best = 0
        cdef int x
        cdef int y

        with nogil:
            for i in range(1, self._m * self._n):
                if self._data[0][i] > self._data[0][best]:
                    best = i
            x = best // self._m
            y = best % self._n
        return x, y

    cpdef tuple argmin(self):
        assert self._data != NULL

        cdef int i
        cdef int best = 0
        cdef int x
        cdef int y

        with nogil:
            for i in range(1, self._m * self._n):
                if self._data[0][i] < self._data[0][best]:
                    best = i
            x = best // self._m
            y = best % self._n
        return x, y

    cpdef MatrixU8 copy(self):

        cdef int      i
        cdef MatrixU8 mat = MatrixU8.__new__(MatrixU8)
        mat._m = mat._shape[0] = self._m
        mat._n = mat._shape[1] = self._n
        mat._strides = self._strides

        with nogil:
            # allocate array of pointers
            mat._data = <uint8_t**> malloc(sizeof(uint8_t*) * self._m)
            if mat._data == NULL:
                raise AllocationError("uint8_t**")
            # allocate memory block
            mat._data[0] = <uint8_t*> malloc(sizeof(uint8_t) * self._m * self._n)
            if mat._data == NULL:
                raise AllocationError("uint8_t*")
            # update array of pointers
            for i in range(self._m):
                mat._data[i] = mat._data[0] + i * self._n
            # copy data
            memcpy(
                <void*> mat._data[0],
                <void*> self._data[0],
                self._m * self._n * sizeof(uint8_t)
            )

        return mat

    cpdef uint8_t max(self):
        assert self._data != NULL

        cdef int     i
        cdef uint8_t best = self._data[0][0]

        with nogil:
            for i in range(1, self._m * self._n):
                if self._data[0][i] > best:
                    best = self._data[0][i]
        return best

    cpdef uint8_t min(self):
        assert self._data != NULL

        cdef int     i
        cdef uint8_t best = self._data[0][0]

        with nogil:
            for i in range(1, self._m * self._n):
                if self._data[0][i] < best:
                    best = self._data[0][i]
        return best

    cpdef uint8_t sum(self):
        assert self._data != NULL

        cdef int     i
        cdef uint8_t sum = 0

        with nogil:
            for i in range(self._m * self._n):
                sum += self._data[0][i]
        return sum


# --- Multiple Sequences Alignment -------------------------------------------

cdef class _MSASequences:
    """A read-only view over the individual sequences of an MSA.
    """

    def __cinit__(self):
        self.msa = None

    def __init__(self, MSA msa):
        self.msa = msa

    def __len__(self):
        assert self.msa._msa != NULL
        return self.msa._msa.nseq


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
        cdef       esl_pos_t length = len(accession)
        cdef const char*     acc    = accession

        with nogil:
            status = libeasel.msa.esl_msa_SetAccession(self._msa, acc, length)
        if status == libeasel.eslEMEM:
            raise AllocationError("char*")
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
        cdef       esl_pos_t length = len(author)
        cdef const char*     au   = author

        with nogil:
            status = libeasel.msa.esl_msa_SetAuthor(self._msa, au, length)
        if status == libeasel.eslEMEM:
            raise AllocationError("char*")
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
        cdef       esl_pos_t length = len(name)
        cdef const char*     nm     = name

        with nogil:
            status = libeasel.msa.esl_msa_SetName(self._msa, nm, length)
        if status == libeasel.eslEMEM:
            raise AllocationError("char*")
        elif status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_msa_SetName")

    @property
    def description(self):
        """`bytes` or `None`: The description of the sequence, if any.
        """
        assert self._msa != NULL
        if self._msa.desc == NULL:
            return None
        return <bytes> self._msa.desc

    @description.setter
    def description(self, bytes description):
        assert self._msa != NULL

        cdef       int       status
        cdef       esl_pos_t length = len(description)
        cdef const char*     desc   = description

        with nogil:
            status = libeasel.msa.esl_msa_SetDesc(self._msa, desc, length)
        if status == libeasel.eslEMEM:
            raise AllocationError("char*")
        elif status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_msa_SetDesc")

    # --- Utils --------------------------------------------------------------

    cdef int _rehash(self) nogil except 1:
        """Rehash the sequence names for faster lookup.
        """
        cdef int status = libeasel.msa.esl_msa_Hash(self._msa)
        if status == libeasel.eslOK:
            return 0
        else:
            raise UnexpectedError(status, "esl_msa_Hash")

    # --- Methods ------------------------------------------------------------

    cpdef uint32_t checksum(self):
        """checksum(self)\n--

        Calculate a 32-bit checksum for the multiple sequence alignment.

        """
        cdef uint32_t checksum = 0
        cdef int status
        with nogil:
            status = libeasel.msa.esl_msa_Checksum(self._msa, &checksum)
        if status == libeasel.eslOK:
            return checksum
        else:
            raise UnexpectedError(status, "esl_msa_Checksum")

    cpdef void write(self, object fh, str format) except *:
        """write(self, fh, format)\n--

        Write the multiple sequence alignement to a file handle.

        Arguments:
            fh (`io.IOBase`): A Python file handle, opened in binary mode.
            format (`str`): The name of the multiple sequence alignment
                file format to use.

        .. versionadded:: 0.3.0

        """
        assert self._msa != NULL

        cdef int    status
        cdef int    fmt
        cdef FILE*  file

        if format not in MSAFile._formats:
            raise ValueError("Invalid MSA format: {!r}".format(format))

        fmt = MSAFile._formats[format]
        file = fopen_obj(fh, mode="w")
        status = libeasel.msafile.esl_msafile_Write(file, self._msa, fmt)
        fclose(file)

        if status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_sqascii_WriteFasta")


cdef class _TextMSASequences(_MSASequences):
    """A read-only view over the sequences of an MSA in text mode.
    """

    def __init__(self, TextMSA msa):
        self.msa = msa

    def __getitem__(self, int idx):
        assert self.msa._msa != NULL

        cdef int          status
        cdef TextSequence seq

        if idx < 0:
            idx += self.msa._msa.nseq
        if idx >= self.msa._msa.nseq or idx < 0:
            raise IndexError("list index out of range")

        seq = TextSequence.__new__(TextSequence)
        status = libeasel.sq.esl_sq_FetchFromMSA(self.msa._msa, idx, &seq._sq)

        if status == libeasel.eslOK:
            return seq
        else:
            raise UnexpectedError(status, "esl_sq_FetchFromMSA")

    def __setitem__(self, int idx, TextSequence seq):
        assert self.msa is not None
        assert seq._sq != NULL

        cdef int status
        cdef int hash_index

        if idx < 0:
            idx += self.msa._msa.nseq
        if idx >= self.msa._msa.nseq or idx < 0:
            raise IndexError("list index out of range")

        # make sure the sequence has a name
        if seq.name is None:
            raise ValueError("cannot set an alignment sequence with an empty name")

        # make sure the sequence has the right length
        if len(seq) != len(self.msa):
            raise ValueError("sequence does not have the expected length")

        # make sure inserting the sequence will not create a name duplicate
        status = libeasel.keyhash.esl_keyhash_Lookup(
            self.msa._msa.index,
            seq._sq.name,
            -1,
            &hash_index
        )
        if status == libeasel.eslOK and hash_index != idx:
            raise ValueError("cannot set a sequence with a duplicate name")

        # set the new sequence
        with nogil:
            (<TextMSA> self.msa)._set_sequence(idx, (<TextSequence> seq)._sq)
            if hash_index != idx:
                self.msa._rehash()


cdef class TextMSA(MSA):
    """A multiple sequence alignement stored in text mode.
    """

    # --- Magic methods ------------------------------------------------------

    def __init__(
        self,
        bytes name=None,
        bytes description=None,
        bytes accession=None,
        object sequences=None,
        bytes author=None,
    ):
        """__init__(self, name=None, description=None, accession=None, sequences=None, author=None)\n--

        Create a new text-mode alignment with the given ``sequences``.

        Arguments:
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

        """
        cdef list    seqs  = [] if sequences is None else list(sequences)

        for seq in seqs:
            if not isinstance(seq, TextSequence):
                ty = type(seq).__name__
                raise TypeError(f"expected TextSequence, found {ty}")

        cdef set     names = {seq.name for seq in seqs}
        cdef int64_t alen  = len(seqs[0]) if seqs else -1
        cdef int     nseq  = len(seqs) if seqs else 1

        if len(names) != len(seqs):
            raise ValueError("duplicate names in alignment sequences")
        elif not all(len(seq) == alen for seq in seqs):
            raise ValueError("all sequences must have the same length")

        with nogil:
            self._msa = libeasel.msa.esl_msa_Create(nseq, alen)
        if self._msa == NULL:
            raise AllocationError("ESL_MSA")

        if name is not None:
            self.name = name
        if accession is not None:
            self.accession = accession
        if description is not None:
            self.description = description
        if author is not None:
            self.author = author
        for i, seq in enumerate(seqs):
            self._set_sequence(i, (<TextSequence> seq)._sq)

    # --- Properties ---------------------------------------------------------

    @property
    def sequences(self):
        """`_TextMSASequences`: A view of the sequences in the alignment.

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

            Support for this feature will be added in a future version, but
            can be circumvented for now by forcingly setting the updated
            version of the object::

                >>> seq = msa.sequences[0]
                >>> seq.name = b"seq1bis"
                >>> msa.sequences[0] = seq
                >>> msa.sequences[0].name
                b'seq1bis'

        .. versionadded:: 0.3.0

        """
        return _TextMSASequences(self)

    # --- Utils --------------------------------------------------------------

    cdef int _set_sequence(self, int idx, ESL_SQ* seq) nogil except 1:
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
        """copy(self)\n--

        Duplicate the text sequence alignment, and return the copy.

        """
        assert self._msa != NULL
        assert not (self._msa.flags & libeasel.msa.eslMSA_DIGITAL)

        cdef int status
        cdef TextMSA new = TextMSA.__new__(TextMSA)
        with nogil:
            new._msa = libeasel.msa.esl_msa_Clone(self._msa)
        if new._msa == NULL:
            raise AllocationError("ESL_MSA")
        return new

    cpdef DigitalMSA digitize(self, Alphabet alphabet):
        """digitize(self, alphabet)\n--

        Convert the text alignment to a digital alignment using ``alphabet``.

        Returns:
            `DigitalMSA`: An alignment in digital mode containing the same
            sequences digitized with ``alphabet``.

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
                raise AllocationError("ESL_MSA")
            status = libeasel.msa.esl_msa_Digitize(alphabet._abc, new._msa, <char*> &errbuf)

        if status == libeasel.eslOK:
            assert new._msa.flags & libeasel.msa.eslMSA_DIGITAL
            return new
        elif status == libeasel.eslEINVAL:
            err_msg = errbuf.decode("utf-8")
            raise ValueError(f"Cannot digitize MSA with alphabet {alphabet}: {err_msg}")
        else:
            raise UnexpectedError(status, "esl_msa_Digitize")


cdef class _DigitalMSASequences(_MSASequences):
    """A read-only view over the sequences of an MSA in digital mode.
    """

    def __init__(self, DigitalMSA msa):
        self.msa = msa
        self.alphabet = msa.alphabet

    def __getitem__(self, int idx):
        assert self.msa._msa != NULL

        cdef int             status
        cdef DigitalSequence seq

        if idx < 0:
            idx += self.msa._msa.nseq
        if idx >= self.msa._msa.nseq or idx < 0:
            raise IndexError("list index out of range")

        seq = DigitalSequence.__new__(DigitalSequence, self.alphabet)
        with nogil:
            status = libeasel.sq.esl_sq_FetchFromMSA(self.msa._msa, idx, &seq._sq)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_sq_FetchFromMSA")

        return seq

    def __setitem__(self, int idx, DigitalSequence seq):
        assert self.msa is not None
        assert seq._sq != NULL

        cdef int status
        cdef int hash_index

        if idx < 0:
            idx += self.msa._msa.nseq
        if idx >= self.msa._msa.nseq or idx < 0:
            raise IndexError("list index out of range")

        # make sure the sequence has a name
        if seq.name is None:
            raise ValueError("cannot set an alignment sequence with an empty name")

        # make sure the sequence has the right length
        if len(seq) != len(self.msa):
            raise ValueError("sequence does not have the expected length")

        # make sure the sequence has the right alphabet
        if not (<Alphabet> self.msa.alphabet)._eq(seq.alphabet):
            raise AlphabetMismatch(self.msa.alphabet, seq.alphabet)

        # make sure inserting the sequence will not create a name duplicate
        status = libeasel.keyhash.esl_keyhash_Lookup(
            self.msa._msa.index,
            seq._sq.name,
            -1,
            &hash_index
        )
        if status == libeasel.eslOK and hash_index != idx:
            raise ValueError("cannot set a sequence with a duplicate name")

        # set the new sequence
        with nogil:
            (<TextMSA> self.msa)._set_sequence(idx, (<TextSequence> seq)._sq)
            if hash_index != idx:
                self.msa._rehash()


cdef class DigitalMSA(MSA):
    """A multiple sequence alignment stored in digital mode.

    Attributes:
        alphabet (`Alphabet`): The biological alphabet used to encode this
            sequence alignment to digits.

    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self, Alphabet alphabet, *args, **kwargs):
        self._msa = NULL
        self.alphabet = alphabet

    def __init__(
        self,
        Alphabet alphabet,
        bytes name=None,
        bytes description=None,
        bytes accession=None,
        object sequences=None,
        bytes author=None,
    ):
        """__init__(self, alphabet, name=None, description=None, accession=None, sequences=None, author=None)\n--

        Create a new digital-mode alignment with the given ``sequences``.

        Arguments:
            alphabet (`Alphabet`): The alphabet of the alignmed sequences.
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

        """
        cdef list    seqs  = [] if sequences is None else list(sequences)
        cdef set     names = { seq.name for seq in seqs }
        cdef int64_t alen  = len(seqs[0]) if seqs else -1
        cdef int     nseq  = len(seqs) if seqs else 1

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
            raise AllocationError("ESL_MSA")

        if name is not None:
            self.name = name
        if accession is not None:
            self.accession = accession
        if description is not None:
            self.description = description
        if author is not None:
            self.author = author
        for i, seq in enumerate(seqs):
            self._set_sequence(i, (<DigitalSequence> seq)._sq)


    # --- Properties ---------------------------------------------------------

    @property
    def sequences(self):
        """`_DigitalMSASequences`: A view of the sequences in the alignment.

        This property lets you access the individual sequences in the
        multiple sequence alignment as `DigitalSequence` instances.

        See Also:
            The documentation for the `TextMSA.sequences` property, which
            contains some additional information.

        .. versionadded:: 0.3.0

        """
        return _DigitalMSASequences(self)

    # --- Utils --------------------------------------------------------------

    cdef int _set_sequence(self, int idx, ESL_SQ* seq) nogil except 1:
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
        memcpy(self._msa.ax[idx], seq.dsq, self._msa.alen+2)
        return 0

    # --- Methods ------------------------------------------------------------

    cpdef DigitalMSA copy(self):
        """Duplicate the digital sequence alignment, and return the copy.
        """
        assert self._msa != NULL
        assert self._msa.flags & libeasel.msa.eslMSA_DIGITAL

        cdef int           status
        cdef ESL_ALPHABET* abc    = self.alphabet._abc
        cdef DigitalMSA    new    = DigitalMSA.__new__(DigitalMSA, self.alphabet)

        with nogil:
            new._msa = libeasel.msa.esl_msa_Clone(self._msa)
        if new._msa == NULL:
            raise AllocationError("ESL_MSA")
        return new

    cpdef TextMSA textize(self):
        """textize(self)\n--

        Convert the digital alignment to a text alignment.

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
                raise AllocationError("ESL_MSA")
            status = libeasel.msa.esl_msa_Textize(new._msa)

        if status == libeasel.eslOK:
            assert not (new._msa.flags & libeasel.msa.eslMSA_DIGITAL)
            return new
        else:
            raise UnexpectedError(status, "esl_msa_Textize")


# --- MSA File ---------------------------------------------------------------

cdef class MSAFile:
    """A wrapper around a multiple-alignment file.

    This class supports reading sequences stored in different formats, such
    as Stockholm, A2M, PSI-BLAST or Clustal.

    """

    _formats = {
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

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self.alphabet = None
        self._msaf = NULL

    def __init__(self, str file, str format = None):
        cdef int fmt = libeasel.msafile.eslMSAFILE_UNKNOWN
        if format is not None:
            format_ = format.lower()
            if format_ not in self._formats:
                raise ValueError("Invalid MSA format: {!r}".format(format))
            fmt = self._formats[format_]

        cdef bytes fspath = os.fsencode(file)
        cdef int status = libeasel.msafile.esl_msafile_Open(NULL, fspath, NULL, fmt, NULL, &self._msaf)
        if status == libeasel.eslENOTFOUND:
            raise FileNotFoundError(2, "No such file or directory: {!r}".format(file))
        elif status == libeasel.eslEMEM:
            raise AllocationError("ESL_MSAFILE")
        elif status == libeasel.eslENOFORMAT:
            if format is None:
                raise ValueError("Could not determine format of file: {!r}".format(file))
            else:
                raise EOFError("Sequence file is empty")
        elif status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_msafile_Open")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def __iter__(self):
        return self

    def __next__(self):
        seq = self.read()
        if seq is None:
            raise StopIteration()
        return seq

    # --- Read methods -------------------------------------------------------

    cpdef MSA read(self):
        """read(self)\n--

        Read the next alignment from the file.

        Returns:
            `MSA`: The next alignment in the file, or `None` if all the
            alignments were read from the file already.

        Raises:
            `ValueError`: When attempting to read an alignment from a closed
                file, or when the file could not be parsed.

        Hint:
            This method allocates a new alignment, which is not efficient in
            case the sequences are being read within a tight loop. Use
            `SequenceFile.readinto` with an already initialized `Sequence`
            if you can to recycle the internal buffers.

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
            msg = <bytes> self._msaf.errmsg
            raise ValueError("Could not parse file: {}".format(msg.decode()))
        else:
            raise UnexpectedError(status, "esl_msafile_Read")

    # --- Utils --------------------------------------------------------------

    cpdef void close(self):
        """close(self)\n--

        Close the file and free the resources used by the parser.

        """
        libeasel.msafile.esl_msafile_Close(self._msaf)
        self._msaf = NULL

    cpdef Alphabet guess_alphabet(self):
        """guess_alphabet(self)\n--

        Guess the alphabet of an open `MSAFile`.

        This method tries to guess the alphabet of a multiple-alignment file
        by inspecting the first entry in the file. It returns the alphabet,
        or `None` if the file alphabet cannot be reliably guessed.

        Raises:
            `EOFError`: if the file is empty.
            `OSError`: if a parse error occurred.
            `ValueError`: if this methods is called after the file was closed.

        """
        cdef int ty
        cdef int status
        cdef Alphabet alphabet

        if self._msaf == NULL:
            raise ValueError("I/O operation on closed file.")

        status = libeasel.msafile.esl_msafile_GuessAlphabet(self._msaf, &ty)
        if status == libeasel.eslOK:
            alphabet = Alphabet.__new__(Alphabet)
            alphabet._init_default(ty)
            return alphabet
        elif status == libeasel.eslENOALPHABET:
            return None
        elif status == libeasel.eslENODATA:
            raise EOFError("Sequence file appears to be empty.")
        elif status == libeasel.eslEFORMAT:
            msg = <bytes> self._msaf.errmsg
            raise ValueError("Could not parse file: {}".format(msg.decode()))
        else:
            raise UnexpectedError(status, "esl_msafile_GuessAlphabet")

    cpdef Alphabet set_digital(self, Alphabet alphabet):
        """set_digital(self, alphabet)\n--

        Set the `MSAFile` to read in digital mode with ``alphabet``.

        This method can be called even after the first alignment have been
        read; it only affects subsequent sequences in the file.

        Returns:
            `~pyhmmer.easel.Alphabet`: The alphabet it was given. Useful to
            wrap `~MSAFile.guess_alphabet` calls and get the resulting
            alphabet.

        .. versionchanged:: 0.4.0
            Returns the `Alphabet` given as argument instead of `None`.

        """
        if self._msaf == NULL:
            raise ValueError("I/O operation on closed file.")

        cdef int status = libeasel.msafile.esl_msafile_SetDigital(self._msaf, alphabet._abc)
        if status == libeasel.eslOK:
            self.alphabet = alphabet
            return alphabet
        else:
            raise UnexpectedError(status, "esl_msafile_SetDigital")


# --- Randomness -------------------------------------------------------------

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
        """__init__(self, seed=None, fast=False)\n--

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
                raise AllocationError("ESL_RANDOMNESS")
        else:
            self._seed(_seed)

    def __copy__(self):
        return self.copy()

    def __getstate__(self):
        return self.getstate()

    def __setstate__(self, state):
        if self._rng == NULL:
            self._rng = <ESL_RANDOMNESS*> malloc(sizeof(ESL_RANDOMNESS))
            if self._rng == NULL:
                raise AllocationError("ESL_RANDOMNESS")
        self.setstate(state)

    def __sizeof__(self):
        assert self._rng != NULL
        return sizeof(ESL_RANDOMNESS) + sizeof(self)

    # --- Methods ------------------------------------------------------------

    def getstate(self):
        """getstate(self)\n--

        Get a tuple containing the current state.

        """
        assert self._rng != NULL
        if self.is_fast():
            return ( True, self._rng.seed, self._rng.x )
        else:
            return ( False, self._rng.seed, self._rng.mti, [self._rng.mt[x] for x in range(624)] )

    def setstate(self, tuple state):
        """setstate(self, state)\n--

        Restores the state of the random number generator.

        """
        assert self._rng != NULL
        if self._rng == NULL:
            self._rng = <ESL_RANDOMNESS*> malloc(sizeof(ESL_RANDOMNESS))
            if self._rng == NULL:
                raise AllocationError("ESL_RANDOMNESS")

        self._rng.seed = state[1]
        if state[0]:
            self._rng.type = libeasel.random.esl_randomness_type.eslRND_FAST
            self._rng.x = state[2]
        else:
            self._rng.type = libeasel.random.esl_randomness_type.eslRND_MERSENNE
            self._rng.mti = state[2]
            for x in range(624):
                self._rng.mt[x] = state[3][x]

    cdef int _seed(self, uint32_t n) except 1:
        status = libeasel.random.esl_randomness_Init(self._rng, n)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_randomness_Init")

    cpdef void seed(self, object n=None):
        """seed([n])\n--

        Reinitialize the random number generator with the given seed.

        Arguments:
            n (`int`, optional): The seed to use. If ``0`` or `None`, an
                arbitrary seed will be chosen using the current time.

        """
        assert self._rng != NULL
        cdef uint32_t seed = n if n is not None else 0
        self._seed(seed)

    cpdef Randomness copy(self):
        """copy(self)\n--

        Return a copy of the random number generator in the same exact state.

        """
        assert self._rng != NULL

        cdef Randomness new = Randomness.__new__(Randomness)
        new._rng = <ESL_RANDOMNESS*> malloc(sizeof(ESL_RANDOMNESS))
        if new._rng == NULL:
            raise AllocationError("ESL_RANDOMNESS")

        memcpy(<void*> new._rng, <void*> self._rng, sizeof(ESL_RANDOMNESS))
        return new

    cpdef double random(self):
        """random(self)\n--

        Generate a uniform random deviate on :math:`\\left[ 0, 1 \\right)`.

        """
        assert self._rng != NULL
        return libeasel.random.esl_random(self._rng)

    cpdef double normalvariate(self, double mu, double sigma):
        """normalvariate(self, mu, sigma)\n--

        Generate a Gaussian-distributed sample.

        Arguments:
            mu (`float`): The mean of the Gaussian being sampled.
            sigma (`float`): The standard deviation of the Gaussian being
                sampled.

        """
        assert self._rng != NULL
        return libeasel.random.esl_rnd_Gaussian(self._rng, mu, sigma)

    cpdef bint is_fast(self):
        """is_fast(self)\n--

        Returns whether or not the linear congruential generator is in use.

        """
        assert self._rng != NULL
        return self._rng.type == libeasel.random.esl_randomness_type.eslRND_FAST


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
    `pyhmmer.hmmsearch`, can then use Python type system to make sure they
    receive sequences in the right mode. This allows type checkers such as
    ``mypy`` to detect potential contract breaches at compile-time.

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
            raise AllocationError("char*")
        elif status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_sq_SetAccession")

    @property
    def name(self):
        """`bytes`: The name of the sequence.
        """
        assert self._sq != NULL
        return <bytes> self._sq.name

    @name.setter
    def name(self, bytes name):
        assert self._sq != NULL

        cdef       int   status
        cdef const char* nm     = name

        with nogil:
            status = libeasel.sq.esl_sq_SetName(self._sq, nm)
        if status == libeasel.eslEMEM:
            raise AllocationError("char*")
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
            raise AllocationError("char*")
        elif status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_sq_SetDesc")

    @property
    def source(self):
        """`bytes`: The source of the sequence, if any.
        """
        return <bytes> self._sq.source

    @source.setter
    def source(self, bytes source):
        assert self._sq != NULL

        cdef       int   status
        cdef const char* src    = source

        with nogil:
            status = libeasel.sq.esl_sq_SetSource(self._sq, src)
        if status == libeasel.eslEMEM:
            raise AllocationError("char*")
        elif status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_sq_SetSource")


    # --- Abstract methods ---------------------------------------------------

    def copy(self):
        """copy(self)\n--

        Duplicate the sequence, and return the copy.

        """
        raise NotImplementedError("Sequence.copy")


    # --- Methods ------------------------------------------------------------

    cpdef uint32_t checksum(self):
        """checksum(self)\n--

        Calculate a 32-bit checksum for the sequence.

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

    cpdef void clear(self):
        """clear(self)\n--

        Reinitialize the sequence for re-use.

        """
        assert self._sq != NULL

        cdef int status

        with nogil:
            status = libeasel.sq.esl_sq_Reuse(self._sq)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_sq_Reuse")

    cpdef void write(self, object fh) except *:
        """write(self, fh)\n--

        Write the sequence alignement to a file handle, in FASTA format.

        Arguments:
            fh (`io.IOBase`): A Python file handle, opened in binary mode.

        .. versionadded:: 0.3.0

        """
        assert self._sq != NULL

        cdef int    status
        cdef FILE*  file   = fopen_obj(fh, mode="w")

        status = libeasel.sqio.ascii.esl_sqascii_WriteFasta(file, self._sq, False)
        fclose(file)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_sqascii_WriteFasta")


cdef class TextSequence(Sequence):
    """A biological sequence stored in text mode.

    Hint:
        Use the ``sequence`` property to access the sequence letters as a
        Python string.

    """

    def __init__(
        self,
        bytes name=None,
        bytes description=None,
        bytes accession=None,
        str   sequence=None,
        bytes source=None,
    ):
        """__init__(self, name=None, description=None, accession=None, sequence=None, source=None)\n--

        Create a new text-mode sequence with the given attributes.

        """
        cdef bytes sq

        if sequence is not None:
            sq = sequence.encode("ascii")
            self._sq = libeasel.sq.esl_sq_CreateFrom(NULL, sq, NULL, NULL, NULL)
        else:
            self._sq = libeasel.sq.esl_sq_Create()
        if self._sq == NULL:
            raise AllocationError("ESL_SQ")
        self._sq.abc = NULL

        if name is not None:
            self.name = name
        if accession is not None:
            self.accession = accession
        if description is not None:
            self.description = description
        if source is not None:
            self.source = source

        assert libeasel.sq.esl_sq_IsText(self._sq)
        assert self._sq.name != NULL
        assert self._sq.desc != NULL
        assert self._sq.acc != NULL

    @property
    def sequence(self):
        """`str`: The raw sequence letters, as a Python string.
        """
        assert self._sq != NULL
        assert libeasel.sq.esl_sq_IsText(self._sq)
        return self._sq.seq.decode("ascii")

    cpdef DigitalSequence digitize(self, Alphabet alphabet):
        """digitize(self, alphabet)\n--

        Convert the text sequence to a digital sequence using ``alphabet``.

        Returns:
            `DigitalSequence`: A copy of the sequence in digital-model,
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
                raise AllocationError("ESL_SQ")

            status = libeasel.sq.esl_sq_Copy(self._sq, new._sq)
            if status != libeasel.eslOK:
                raise UnexpectedError(status, "esl_sq_Copy")

        assert libeasel.sq.esl_sq_IsDigital(new._sq)
        return new

    cpdef TextSequence copy(self):
        """copy(self)\n--

        Duplicate the text sequence, and return the copy.

        """
        assert self._sq != NULL
        assert libeasel.sq.esl_sq_IsText(self._sq)

        cdef int          status
        cdef TextSequence new    = TextSequence.__new__(TextSequence)

        with nogil:
            new._sq = libeasel.sq.esl_sq_Create()
            if new._sq == NULL:
                raise AllocationError("ESL_SQ")
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
            UserWarning: When the sequence contains unknown characters.

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

    """

    def __cinit__(self, Alphabet alphabet, *args, **kwargs):
        self.alphabet = alphabet

    def __init__(self,
        Alphabet alphabet,
        bytes    name=None,
        bytes    description=None,
        bytes    accession=None,
        char[:]  sequence=None,
        bytes    source=None,
    ):
        """__init__(self, alphabet, name=None, description=None, accession=None, sequence=None, source=None)\n--

        Create a new digital-mode sequence with the given attributes.

        .. versionadded:: 0.1.4

        """
        cdef int     status
        cdef int64_t n

        # create an empty digital sequence
        self._sq = libeasel.sq.esl_sq_CreateDigital(alphabet._abc)
        if self._sq == NULL:
            raise AllocationError("ESL_SQ")

        # NB: because the easel sequence has sentinel bytes that we hide from
        #     the user, we cannot just copy the sequence here or use the libeasel
        #     internals; instead, if a sequence is given, we need to emulate
        #     the `esl_sq_CreateDigitalFrom` but copy the sequence with different
        #     offsets.
        if sequence is not None:
            # we can release the GIL while copying memory
            with nogil:
                # grow the sequence buffer so it can hold `n` residues
                n = sequence.shape[0]
                status = libeasel.sq.esl_sq_GrowTo(self._sq, n)
                if status != libeasel.eslOK:
                    raise UnexpectedError(status, "esl_sq_GrowTo")
                # update the digital sequence buffer
                self._sq.dsq[0] = self._sq.dsq[n+1] = libeasel.eslDSQ_SENTINEL
                memcpy(<void*> &self._sq.dsq[1], <const void*> &sequence[0], n)
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

        assert libeasel.sq.esl_sq_IsDigital(self._sq)
        assert self._sq.name != NULL
        assert self._sq.desc != NULL
        assert self._sq.acc != NULL

    @property
    def sequence(self):
        """`VectorU8`: The raw sequence digits, as a byte vector.

        Note:
            The internal ``ESL_SQ`` object allocates a buffer of size
            :math:`n+2` (where :math:`n` is the number of residues in the
            sequence), with the first and the last element of the buffer
            being sentinel values. This vector does not expose the sentinel
            values, only the :math:`n` elements of the buffer in between.

        .. versionchanged:: v0.4.0
           Property is now a `VectorU8` instead of a memoryview.

        """
        assert self._sq != NULL

        cdef VectorU8 seq = VectorU8.__new__(VectorU8)
        seq._n = seq._shape[0] = self._sq.n
        seq._data = &self._sq.dsq[1]
        seq._owner = self

        return seq

    cpdef DigitalSequence copy(self):
        """copy(self)\n--

        Duplicate the digital sequence, and return the copy.

        """
        assert self._sq != NULL
        assert libeasel.sq.esl_sq_IsDigital(self._sq)

        cdef int             status
        cdef ESL_ALPHABET*   abc    = self.alphabet._abc
        cdef DigitalSequence new    = DigitalSequence.__new__(DigitalSequence, self.alphabet)

        with nogil:
            new._sq = libeasel.sq.esl_sq_CreateDigital(abc)
            if new._sq == NULL:
                raise AllocationError("ESL_SQ")

            status = libeasel.sq.esl_sq_Copy(self._sq, new._sq)
            if status != libeasel.eslOK:
                raise UnexpectedError(status, "esl_sq_Copy")

        assert libeasel.sq.esl_sq_IsDigital(new._sq)
        return new

    cpdef TextSequence textize(self):
        """textize(self)\n--

        Convert the digital sequence to a text sequence.

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
                raise AllocationError("ESL_SQ")

            status = libeasel.sq.esl_sq_Copy(self._sq, new._sq)
            if status != libeasel.eslOK:
                raise UnexpectedError(status, "esl_sq_Copy")

        assert libeasel.sq.esl_sq_IsText(new._sq)
        return new

    cpdef DigitalSequence reverse_complement(self, bint inplace=False):
        """Build the reverse complement of the sequence.

        Arguments:
            inplace (`bool`): Whether or not to copy the sequence before
                computing its reverse complement. With `False` (the default),
                the method will return a copy of the sequence that has been
                reverse-complemented. With `True`, it will reverse-complement
                inplace and return `None`.

        Raises:
            ValueError: When the alphabet of the `DigitalSequence` does
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


# --- Sequence File ----------------------------------------------------------

cdef class SequenceFile:
    """A wrapper around a sequence file, containing unaligned sequences.

    This class supports reading sequences stored in different formats, such
    as FASTA, GenBank or EMBL. The format of each file can be automatically
    detected, but it is also possible to pass an explicit format specifier
    when the `SequenceFile` is instantiated.

    .. versionadded:: 0.2.0
       The ``alphabet`` attribute.

    """

    _formats = {
        "fasta": libeasel.sqio.eslSQFILE_FASTA,
        "embl": libeasel.sqio.eslSQFILE_EMBL,
        "genbank": libeasel.sqio.eslSQFILE_GENBANK,
        "ddbj": libeasel.sqio.eslSQFILE_DDBJ,
        "uniprot": libeasel.sqio.eslSQFILE_UNIPROT,
        "ncbi": libeasel.sqio.eslSQFILE_NCBI,
        "daemon": libeasel.sqio.eslSQFILE_DAEMON,
        "hmmpgmd": libeasel.sqio.eslSQFILE_DAEMON,
        "fmindex": libeasel.sqio.eslSQFILE_FMINDEX,
    }


    # --- Class methods ------------------------------------------------------

    @classmethod
    def parse(cls, bytes buffer, str format):
        """parse(cls, buffer, format)\n--

        Parse a sequence from a binary ``buffer`` using the given ``format``.

        """
        cdef Sequence seq = TextSequence.__new__(TextSequence)
        seq._sq = libeasel.sq.esl_sq_Create()
        if not seq._sq:
            raise AllocationError("ESL_SQ")
        seq._sq.abc = NULL
        return cls.parseinto(seq, buffer, format)

    @classmethod
    def parseinto(cls, Sequence seq, bytes buffer, str format):
        """parseinto(cls, seq, buffer, format)\n--

        Parse a sequence from a binary ``buffer`` into ``seq``.

        """
        assert seq._sq != NULL

        cdef int fmt = libeasel.sqio.eslSQFILE_UNKNOWN
        if format is not None:
            format_ = format.lower()
            if format_ not in cls._formats:
                raise ValueError("Invalid sequence format: {!r}".format(format))
            fmt = cls._formats[format_]

        cdef int status = libeasel.sqio.esl_sqio_Parse(buffer, len(buffer), seq._sq, fmt)
        if status == libeasel.eslEFORMAT:
            raise AllocationError("ESL_SQFILE")
        elif status == libeasel.eslOK:
            return seq
        else:
            raise UnexpectedError(status, "esl_sqio_Parse")


    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self.alphabet = None
        self._sqfp = NULL

    def __init__(self, str file, str format=None, bint ignore_gaps=False):
        """__init__(self, file, format=None, ignore_gaps=False)\n--

        Create a new sequence file parser wrapping the given ``file``.

        Arguments:
            file (`str`): The path to a file containing sequences in one of
                the supported file formats.
            format (`str`, optional): The format of the file, or `None` to
                autodetect. Supported values are: ``fasta``, ``embl``,
                ``genbank``, ``ddbj``, ``uniprot``, ``ncbi``, ``daemon``,
                ``hmmpgmd``, ``fmindex``.
            ignore_gaps (`bool`): When set to `True`, allow ignoring gap
                characters ('-') when they are present in ungapped formats
                such as ``fasta``. With `False`, stick to the default Easel
                behaviour.

        .. versionchanged:: 0.4.4
           Added the ``ignore_gaps`` parameter.

        """
        # get format from string passed as input
        cdef int fmt = libeasel.sqio.eslSQFILE_UNKNOWN
        if format is not None:
            format_ = format.lower()
            if format_ not in self._formats:
                raise ValueError("Invalid sequence format: {!r}".format(format))
            fmt = self._formats[format_]

        # open the given filename
        cdef bytes fspath = os.fsencode(file)
        cdef int status = libeasel.sqio.esl_sqfile_Open(fspath, fmt, NULL, &self._sqfp)
        if status == libeasel.eslENOTFOUND:
            raise FileNotFoundError(2, "No such file or directory: {!r}".format(file))
        elif status == libeasel.eslEMEM:
            raise AllocationError("ESL_SQFILE")
        elif status == libeasel.eslEFORMAT:
            if format is None:
                raise ValueError("Could not determine format of file: {!r}".format(file))
            else:
                raise EOFError("Sequence file is empty")
        elif status != libeasel.eslOK:
            raise UnexpectedError(status, "esl_sqfile_Open")

        # HACK(@althonos): allow ignoring the gap character if explicitly
        #                  requested, which is normally not allowed by
        #                  Easel for ungapped formats (althonos/pyhmmer#7).
        if ignore_gaps:
            libeasel.sqio.esl_sqio_Ignore(self._sqfp, b"-")

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
        seq = self.read()
        if seq is None:
            raise StopIteration()
        return seq


    # --- Read methods -------------------------------------------------------

    cpdef Sequence read(self, bint skip_info=False, bint skip_sequence=False):
        """read(self, skip_info=False, skip_sequence=False)\n--

        Read the next sequence from the file.

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
        if self.alphabet is None:
            seq = TextSequence()
        else:
            seq = DigitalSequence(self.alphabet)
        return self.readinto(seq, skip_info=skip_info, skip_sequence=skip_sequence)

    cpdef Sequence readinto(self, Sequence seq, bint skip_info=False, bint skip_sequence=False):
        """readinto(self, seq, skip_info=False, skip_sequence=False)\n--

        Read the next sequence from the file, using ``seq`` to store data.

        Arguments:
            seq (`~pyhmmer.easel.Sequence`): A sequence object to use to store
                the next entry in the file. If this sequence was used before,
                it must be properly reset (using the `Sequence.clear` method)
                before using it again with `readinto`.
            skip_info (`bool`): Pass `True` to disable reading the sequence
                *metadata*, and only read the sequence *letters*. Defaults to
                False`.
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

            >>> with SequenceFile("vendor/hmmer/testsuite/ecori.fa") as sf:
            ...     seq = TextSequence()
            ...     while sf.readinto(seq) is not None:
            ...         # ... process seq here ... #
            ...         seq.clear()

        """
        assert seq._sq != NULL

        cdef int (*funcread)(ESL_SQFILE *sqfp, ESL_SQ *sq) nogil
        cdef str   funcname
        cdef int   status

        if not skip_info and not skip_sequence:
            funcname = "esl_sqio_Read"
            funcread = libeasel.sqio.esl_sqio_Read
        elif not skip_info:
            funcname = "esl_sqio_ReadInfo"
            funcread = libeasel.sqio.esl_sqio_ReadInfo
        elif not skip_sequence:
            funcname = "esl_sqio_ReadSequence"
            funcread = libeasel.sqio.esl_sqio_ReadSequence
        else:
            raise ValueError("Cannot skip reading both sequence and metadata.")

        if self._sqfp == NULL:
            raise ValueError("I/O operation on closed file.")
        else:
            status = funcread(self._sqfp, seq._sq)

        if status == libeasel.eslOK:
            return seq
        elif status == libeasel.eslEOF:
            return None
        elif status == libeasel.eslEFORMAT:
            msg = <bytes> libeasel.sqio.esl_sqfile_GetErrorBuf(self._sqfp)
            raise ValueError("Could not parse file: {}".format(msg.decode()))
        else:
            raise UnexpectedError(status, funcname)


    # --- Fetch methods ------------------------------------------------------

    # cpdef Sequence fetch(self, bytes key, bint skip_info=False, bint skip_sequence=False):
    #     cdef Sequence seq = TextSequence()
    #     return self.fetchinto(seq, key, skip_info=skip_info, skip_sequence=skip_sequence)
    #
    # cpdef Sequence fetchinto(self, Sequence seq, bytes key, bint skip_info=False, bint skip_sequence=False):
    #     raise NotImplementedError("TODO SequenceFile.fetchinto")


    # --- Utils --------------------------------------------------------------

    cpdef void close(self):
        """close(self)\n--

        Close the file and free the resources used by the parser.

        """
        libeasel.sqio.esl_sqfile_Close(self._sqfp)
        self._sqfp = NULL

    cpdef Alphabet guess_alphabet(self):
        """guess_alphabet(self)\n--

        Guess the alphabet of an open `SequenceFile`.

        This method tries to guess the alphabet of a sequence file by
        inspecting the first sequence in the file. It returns the alphabet,
        or `None` if the file alphabet cannot be reliably guessed.

        Raises:
            `EOFError`: if the file is empty.
            `OSError`: if a parse error occurred.
            `ValueError`: if this methods is called after the file was closed.

        """
        cdef int ty
        cdef int status
        cdef Alphabet alphabet

        if self._sqfp == NULL:
            raise ValueError("I/O operation on closed file.")

        status = libeasel.sqio.esl_sqfile_GuessAlphabet(self._sqfp, &ty)
        if status == libeasel.eslOK:
            alphabet = Alphabet.__new__(Alphabet)
            alphabet._init_default(ty)
            return alphabet
        elif status == libeasel.eslENOALPHABET:
            return None
        elif status == libeasel.eslENODATA:
            raise EOFError("Sequence file appears to be empty.")
        elif status == libeasel.eslEFORMAT:
            msg = <bytes> libeasel.sqio.esl_sqfile_GetErrorBuf(self._sqfp)
            raise ValueError("Could not parse file: {}".format(msg.decode()))
        else:
            raise UnexpectedError(status, "esl_sqfile_GuessAlphabet")

    cpdef Alphabet set_digital(self, Alphabet alphabet):
        """set_digital(self, alphabet)\n--

        Set the `SequenceFile` to read in digital mode with ``alphabet``.

        This method can be called even after the first sequences have been
        read; it only affects subsequent sequences in the file.

        Returns:
            `~pyhmmer.easel.Alphabet`: The alphabet it was given. Useful to
            wrap `~SequenceFile.guess_alphabet` calls and get the resulting
            alphabet.

        .. versionchanged:: 0.4.0
            Returns the `Alphabet` given as argument instead of `None`.

        """
        if self._sqfp == NULL:
            raise ValueError("I/O operation on closed file.")

        cdef int status = libeasel.sqio.esl_sqfile_SetDigital(self._sqfp, alphabet._abc)
        if status == libeasel.eslOK:
            self.alphabet = alphabet
            return self.alphabet
        else:
            raise UnexpectedError(status, "esl_sqfile_SetDigital")


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

    def __init__(self, str file):
        """__init__(self, file)\n--

        Create a new SSI file reader for the file at the given location.

        Arguments:
            file (`str`): The path to a sequence/subsequence index file to
                read.

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
        """file_info(self, fd)\n--

        Retrieve the `~pyhmmer.easel.SSIReader.FileInfo` of the descriptor.

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
        """find_name(self, key)\n--

        Retrieve the `~pyhmmer.easel.SSIReader.Entry` for the given name.

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
        """close(self)\n--

        Close the SSI file reader.

        """
        libeasel.ssi.esl_ssi_Close(self._ssi)
        self._ssi = NULL


cdef class SSIWriter:
    """A writer for sequence/subsequence index files.
    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._newssi = NULL

    def __init__(self, str file, bint exclusive = False):
        """__init__(self, file)\n--

        Create a new SSI file write for the file at the given location.

        Arguments:
            file (`str`): The path to a sequence/subsequence index file to
                write.
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

    cpdef uint16_t add_file(self, str filename, int format = 0):
        """add_file(self, filename, format=0)\n--

        Add a new file to the index.

        Arguments:
            filename (str): The name of the file to register.
            format (int): A format code to associate with the file, or *0*.

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
            raise UnexpectedError(status, "esl_newssi_AddFile")

    cpdef void add_key(
        self,
        bytes key,
        uint16_t fd,
        off_t record_offset,
        off_t data_offset = 0,
        int64_t record_length = 0
    ):
        """add_key(self, key, fd, record_offset, data_offset=0, record_length=0)\n--

        Add a new entry to the index with the given ``key``.

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
            raise UnexpectedError(status, "esl_newssi_AddKey")

    cpdef void add_alias(self, bytes alias, bytes key):
        """add_alias(self, alias, key)\n--

        Make ``alias`` an alias of ``key`` in the index.

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
            raise UnexpectedError(status, "esl_newssi_AddAlias")

    def close(self):
        """close(self)\n--

        Close the SSI file writer.

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


# --- Module init code -------------------------------------------------------

include "exceptions.pxi"
