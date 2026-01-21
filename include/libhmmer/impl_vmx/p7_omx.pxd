cdef extern from "impl_vmx/impl_vmx.h" nogil:

    ctypedef p7_omx_s P7_OMX
    cdef struct p7_omx_s:
        pass

    P7_OMX* p7_omx_Create(int allocM, int allocL, int allocXL)
    int     p7_omx_GrowTo(P7_OMX *ox, int allocM, int allocL, int allocXL)
    int     p7_omx_Reuse(P7_OMX *ox)
    void    p7_omx_Destroy(P7_OMX *ox)
