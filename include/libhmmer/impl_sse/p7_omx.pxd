from libc.stdio cimport FILE

from libhmmer.impl_sse.p7_oprofile cimport P7_OPROFILE


cdef extern from "emmintrin.h":

    ctypedef struct __m128i:
        pass

    ctypedef struct __m128:
        pass


cdef extern from "impl_sse/impl_sse.h" nogil:

    ctypedef struct P7_OM_BLOCK:
        int count
        int listSize
        P7_OPROFILE** list

    P7_OM_BLOCK *p7_oprofile_CreateBlock(int size)
    void p7_oprofile_DestroyBlock(P7_OM_BLOCK *block)

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

        bint     debugging;
        FILE   *dfp;


    P7_OMX* p7_omx_Create(int allocM, int allocL, int allocXL)
    int     p7_omx_GrowTo(P7_OMX *ox, int allocM, int allocL, int allocXL)
    int     p7_omx_Reuse(P7_OMX *ox)
    void    p7_omx_Destroy(P7_OMX *ox)
