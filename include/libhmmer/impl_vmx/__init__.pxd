from libeasel cimport ESL_DSQ
from libhmmer.impl_vmx.p7_oprofile cimport P7_OPROFILE
from libhmmer.impl_vmx.p7_omx cimport P7_OMX


cdef extern from "impl_vmx/impl_vmx.h" nogil:

    const size_t p7O_EXTRA_SB

    cdef void impl_Init()
    cdef int p7_MSVFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc)
