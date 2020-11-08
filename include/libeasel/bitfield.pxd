from libc.stdint cimport uint64_t


cdef extern from "esl_bitfield.h" nogil:

    ctypedef struct ESL_BITFIELD:
        uint64_t* b
        int nb

    void esl_bitfield_Set(ESL_BITFIELD *b, int i)
    void esl_bitfield_Clear(ESL_BITFIELD *b, int i)
    void esl_bitfield_Toggle(ESL_BITFIELD *b, int i)
    bint esl_bitfield_IsSet(const ESL_BITFIELD *b, int i)

    ESL_BITFIELD* esl_bitfield_Create(int nb)
    int esl_bitfield_Count(const ESL_BITFIELD *b)
    void esl_bitfield_Destroy(ESL_BITFIELD *b)
