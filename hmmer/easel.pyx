cimport libeasel
cimport libeasel.alphabet

from libeasel.alphabet cimport ESL_ALPHABET


cdef class Alphabet:

    cdef _init_default(self, int ty):
        self._alphabet = libeasel.alphabet.esl_alphabet_Create(ty)
        if not self._alphabet:
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

    def __cinit__(self):
        self._alphabet = NULL

    def __dealloc__(self):
        libeasel.alphabet.esl_alphabet_Destroy(self._alphabet)

    def __repr__(self):
        if self._alphabet.type == libeasel.alphabet.eslRNA:
            return "Alphabet.rna()"
        elif self._alphabet.type == libeasel.alphabet.eslDNA:
            return "Alphabet.dna()"
        elif self._alphabet.type == libeasel.alphabet.eslAMINO:
            return "Alphabet.amino()"
        else:
            return "Alphabet({!r}, K={!r}, Kp={!r})".format(
                self._alphabet.sym.decode('ascii'),
                self._alphabet.K,
                self._alphabet.Kp
            )
