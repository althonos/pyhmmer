from posix.types cimport off_t

from libeasel.alphabet cimport ESL_ALPHABET


cdef extern from "p7_config.h" nogil:
    DEF p7_MAXABET = 20


cdef extern from "hmmer.h" nogil:
    DEF p7_NEVPARAM = 6
    DEF p7_NCUTOFFS = 6
    DEF p7_NOFFSETS = 3

    DEF p7P_NXSTATES = 4
    cdef enum p7p_xstates_e:
        p7P_E = 0
        p7P_N = 1
        p7P_J = 2
        p7P_C = 3

    DEF p7P_NXTRANS = 2
    cdef enum p7p_xtransitions_e:
        p7P_LOOP = 0
        p7P_MOVE = 1

    DEF p7P_NTRANS = 8
    cdef enum p7p_tsc_e:
        p7P_MM = 0
        p7P_IM = 1
        p7P_DM = 2
        p7P_BM = 3
        p7P_MD = 4
        p7P_DD = 5
        p7P_MI = 6
        p7P_II = 7

    DEF p7P_NR = 2
    cdef enum p7p_rsc_e:
        p7P_MSC = 0
        p7P_ISC = 1

    cdef enum:
        p7_NO_MODE   = 0
        p7_LOCAL     = 1
        p7_GLOCAL    = 2
        p7_UNILOCAL  = 3
        p7_UNIGLOCAL = 4

    ctypedef p7_profile_s P7_PROFILE
    cdef struct p7_profile_s:
        float* tsc
        float** rsc
        float[p7P_NXSTATES][p7P_NXTRANS] xsc

        int mode
        int L
        int allocM
        int M
        int max_length
        float nj

        char* name
        char* acc
        char* desc
        char* rf
        char* mm
        char* cs
        char* consensus
        float[p7_NEVPARAM] evparam
        float[p7_NCUTOFFS] cutoff
        float[p7_MAXABET] float

        off_t[p7_NOFFSETS] offs
        off_t roff
        off_t eoff

        const ESL_ALPHABET* abc


    P7_PROFILE *p7_profile_Create(int M, const ESL_ALPHABET *abc)
    P7_PROFILE *p7_profile_Clone(const P7_PROFILE *gm)
    int         p7_profile_Copy(const P7_PROFILE *src, P7_PROFILE *dst)
    int         p7_profile_SetNullEmissions(P7_PROFILE *gm)
    int         p7_profile_Reuse(P7_PROFILE *gm)
    size_t      p7_profile_Sizeof(P7_PROFILE *gm)
    void        p7_profile_Destroy(P7_PROFILE *gm)
    int         p7_profile_IsLocal(const P7_PROFILE *gm)
    int         p7_profile_IsMultihit(const P7_PROFILE *gm)
    int         p7_profile_GetT(const P7_PROFILE *gm, char st1, int k1, char st2, int k2, float *ret_tsc);
    int         p7_profile_Validate(const P7_PROFILE *gm, char *errbuf, float tol)
    int         p7_profile_Compare(P7_PROFILE *gm1, P7_PROFILE *gm2, float tol)
