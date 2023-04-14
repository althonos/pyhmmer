from libeasel cimport ESL_DSQ
from libhmmer.impl_sse.p7_oprofile cimport P7_OPROFILE


cdef extern from "impl_sse/impl_sse.h" nogil:

    const size_t p7O_EXTRA_SB

    cdef void impl_Init()
    cdef int p7_SSVFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc)
