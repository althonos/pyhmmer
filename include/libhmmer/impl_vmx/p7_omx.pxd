from libhmmer.impl_vmx.p7_oprofile cimport P7_OPROFILE

cdef extern from "impl_vmx/impl_vmx.h" nogil:

    ctypedef struct P7_OM_BLOCK:
        int count
        int listSize
        P7_OPROFILE** list

    ctypedef p7_omx_s P7_OMX
    cdef struct p7_omx_s:
        pass

    P7_OMX* p7_omx_Create(int allocM, int allocL, int allocXL)
    int     p7_omx_GrowTo(P7_OMX *ox, int allocM, int allocL, int allocXL)
    int     p7_omx_Reuse(P7_OMX *ox)
    void    p7_omx_Destroy(P7_OMX *ox)
