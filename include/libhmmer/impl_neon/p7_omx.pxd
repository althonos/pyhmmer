from libc.stdio cimport FILE

from libhmmer.impl_neon.p7_oprofile cimport P7_OPROFILE


cdef extern from "<arm_neon.h>":

    ctypedef struct float32x4_t:
        pass

    ctypedef struct uint8x16_t:
        pass
    
    ctypedef struct int16x8_t:
        pass
        

cdef extern from "impl_neon/impl_neon.h" nogil:

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

        float32x4_t **dpf
        int16x8_t   **dpw
        uint8x16_t  **dpb
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
