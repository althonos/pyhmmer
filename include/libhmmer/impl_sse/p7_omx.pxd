from libc.stdio cimport FILE


cdef extern from "emmintrin.h":

    ctypedef struct __m128i:
        pass

    ctypedef struct __m128:
        pass


cdef extern from "impl_sse/impl_sse.h" nogil:

    ctypedef p7_omx_s P7_OMX
    ctypedef struct p7_omx_s:
        int       M
        int       L

        __m128  **dpf
        __m128i **dpw
        __m128i **dpb
        void     *dp_mem
        int       allocR
        int       validR
        int       allocQ4
        int       allocQ8
        int       allocQ16
        size_t    ncells

        float    *xmx
        void     *x_mem
        int       allocXR
        float     totscale
        bint       has_own_scales

        bint     debugging
        FILE   *dfp


    P7_OMX* p7_omx_Create(int allocM, int allocL, int allocXL)
    int     p7_omx_GrowTo(P7_OMX *ox, int allocM, int allocL, int allocXL)
    int     p7_omx_Reuse(P7_OMX *ox)
    void    p7_omx_Destroy(P7_OMX *ox)
