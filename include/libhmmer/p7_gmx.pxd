from libc.stdint cimport uint64_t
from libc.stdio cimport FILE


cdef extern from "hmmer.h" nogil:

    cdef size_t p7G_NSCELLS = 3
    cdef enum p7g_scells_e:
        p7G_M = 0
        p7G_I = 1
        p7G_D = 2

    cdef size_t p7G_NXCELLS = 5
    cdef enum p7g_xcells_e:
        p7G_E  = 0
        p7G_N  = 1
        p7G_J  = 2
        p7G_B  = 3
        p7G_C  = 4

    ctypedef p7_gmx_s P7_GMX
    cdef struct p7_gmx_s:
        int  M
        int  L
        int      allocR
        int      validR
        int      allocW
        uint64_t ncells
        float **dp
        float  *xmx
        float  *dp_mem

    P7_GMX *p7_gmx_Create (int allocM, int allocL)
    int     p7_gmx_GrowTo (P7_GMX *gx, int allocM, int allocL)
    size_t  p7_gmx_Sizeof (P7_GMX *gx)
    int     p7_gmx_Reuse  (P7_GMX *gx)
    void    p7_gmx_Destroy(P7_GMX *gx)
    int     p7_gmx_Compare(P7_GMX *gx1, P7_GMX *gx2, float tolerance)
    int     p7_gmx_Dump(FILE *fp, P7_GMX *gx, int flags)
    int     p7_gmx_DumpWindow(FILE *fp, P7_GMX *gx, int istart, int iend, int kstart, int kend, int show_specials)
