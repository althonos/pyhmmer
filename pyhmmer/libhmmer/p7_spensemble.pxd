cdef extern from "hmmer.h" nogil:

    ctypedef struct p7_spcoord_s:
        int idx
        int i
        int j
        int k
        int m
        float prob

    ctypedef p7_spensemble_s P7_SPENSEMBLE
    cdef struct p7_spensemble_s:
        int nsamples
        p7_spcoord_s* sp
        int nalloc
        int n

        int* workspace
        int* assignment
        int nc

        int* epc
        int epc_alloc

        p7_spcoord_s* sigc
        int nsigc
        int nsigc_alloc


    P7_SPENSEMBLE *p7_spensemble_Create(int init_n, int init_epc, int init_sigc)
    int            p7_spensemble_Reuse(P7_SPENSEMBLE *sp)
    int            p7_spensemble_Add(P7_SPENSEMBLE *sp, int sampleidx, int i, int j, int k, int m)
    int            p7_spensemble_Cluster(P7_SPENSEMBLE *sp,
    				          float min_overlap, int of_smaller, int max_diagdiff,
    				          float min_posterior, float min_endpointp,
    				          int *ret_nclusters)
    int     p7_spensemble_GetClusterCoords(P7_SPENSEMBLE *sp, int which,
    					        int *ret_i, int *ret_j, int *ret_k, int *ret_m, float *ret_p)
    void    p7_spensemble_Destroy(P7_SPENSEMBLE *sp)
