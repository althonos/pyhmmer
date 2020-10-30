# coding: utf-8
# cython: language_level=3, linetrace=True

from libc.stdint cimport uint32_t

cimport libeasel
cimport libeasel.alphabet
cimport libeasel.sq


cdef class Alphabet:

    # --- Default constructors -----------------------------------------------

    cdef _init_default(self, int ty):
        self._abc = libeasel.alphabet.esl_alphabet_Create(ty)
        if not self._abc:
            raise MemoryError("could not allocate ESL_ALPHABET")

    @classmethod
    def amino(cls):
        """Create a default Aminoacid alphabet.
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

    # def __init__(self, str alphabet, int K, int Kp):
    #     buffer = alphabet.encode('ascii')
    #     self._alphabet = libeasel.alphabet.esl_alphabet_CreateCustom(<char*> buffer, K, Kp)
    #     if not self._alphabet:
    #         raise MemoryError("could not allocate ESL_ALPHABET")

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._abc = NULL

    def __dealloc__(self):
        libeasel.alphabet.esl_alphabet_Destroy(self._abc)

    def __repr__(self):
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



cdef class Sequence:

    def __cinit__(self):
        self._sq = NULL

    def __init__(self, name, seq, description=None, accession=None, secondary_structure=None):
        # FIXME
        self._sq = libeasel.sq.esl_sq_CreateFrom(name, seq, NULL, NULL, NULL)
        if not self._sq:
            raise MemoryError("could not allocate ESL_SQ")

    def __dealloc__(self):
        libeasel.sq.esl_sq_Destroy(self._sq)

    def __eq__(self, Sequence other):
        return libeasel.sq.esl_sq_Compare(self._sq, other._sq) == libeasel.eslOK

    # --- Properties ---------------------------------------------------------

    @property
    def accession(self):
        """`str` or `None`: The accession of the sequence, if any.
        """
        accession = <bytes> self._sq.acc
        return accession or None

    @accession.setter
    def accession(self, accession):
        if accession is None:
            accession = b""
        cdef int status = libeasel.sq.esl_sq_SetAccession(self._sq, <const char*> accession)
        if status == libeasel.eslEMEM:
            raise MemoryError("could not allocate string")
        elif status != libeasel.eslOK:
            raise RuntimeError("unexpected error in esl_sq_SetAccession")

    @property
    def name(self):
        """`str` or `None`: The name of the sequence, if any.
        """
        name = <bytes> self._sq.name
        return name or None

    @name.setter
    def name(self, name):
        if name is None:
            name = b""
        cdef int status = libeasel.sq.esl_sq_SetName(self._sq, <const char*> name)
        if status == libeasel.eslEMEM:
            raise MemoryError("could not allocate string")
        elif status != libeasel.eslOK:
            raise RuntimeError("unexpected error in esl_sq_SetName")

    @property
    def description(self):
        """`str` or `None`: The description of the sequence, if any.
        """
        desc = <bytes> self._sq.desc
        return desc or None

    @description.setter
    def description(self, desc):
        if desc is None:
            desc = b""
        cdef int status = libeasel.sq.esl_sq_SetDesc(self._sq, <const char*> desc)
        if status == libeasel.eslEMEM:
            raise MemoryError("could not allocate string")
        elif status != libeasel.eslOK:
            raise RuntimeError("unexpected error in esl_sq_SetDesc")

    @property
    def source(self):
        """`str` or `None`: The source of the sequence, if any.
        """
        desc = <bytes> self._sq.source
        return desc or None

    @source.setter
    def source(self, src):
        if src is None:
            src = b""
        cdef int status = libeasel.sq.esl_sq_SetDesc(self._sq, <const char*> src)
        if status == libeasel.eslEMEM:
            raise MemoryError("could not allocate string")
        elif status != libeasel.eslOK:
            raise RuntimeError("unexpected error in esl_sq_SetDesc")

    # --- Methods ------------------------------------------------------------

    def checksum(self):
        """Calculate a 32-bit checksum for a sequence.
        """
        cdef uint32_t checksum = 0
        cdef int status = libeasel.sq.esl_sq_Checksum(self._sq, &checksum)
        if status == libeasel.eslOK:
            return checksum
        else:
            raise RuntimeError("unexpected error in esl_sq_Checksum")
