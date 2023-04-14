from libc.stdio cimport FILE

from libeasel.alphabet cimport ESL_ALPHABET
from libeasel.msa cimport ESL_MSA
from libeasel.sq cimport ESL_DSQ
from libhmmer.p7_profile cimport P7_PROFILE
from libhmmer.p7_hmm cimport P7_HMM
from libhmmer.p7_gmx cimport P7_GMX


cdef extern from "hmmer.h" nogil:

    const size_t p7T_NSTATETYPES
    cdef enum p7t_statetype_e:
        p7T_BOGUS =  0
        p7T_M     =  1
        p7T_D     =  2
        p7T_I     =  3
        p7T_S     =  4
        p7T_N     =  5
        p7T_B     =  6
        p7T_E     =  7
        p7T_C     =  8
        p7T_T     =  9
        p7T_J     = 10
        p7T_X     = 11


    ctypedef p7_trace_s P7_TRACE
    cdef struct p7_trace_s:
        int    N
        int    nalloc
        char  *st
        int   *k
        int   *i
        float *pp
        int    M
        int    L
        int   ndom
        int  *tfrom
        int  *tto
        int  *sqfrom
        int  *sqto
        int  *hmmfrom
        int  *hmmto
        int   ndomalloc


    P7_TRACE *p7_trace_Create()
    P7_TRACE *p7_trace_CreateWithPP()
    int  p7_trace_Reuse(P7_TRACE *tr)
    int  p7_trace_Grow(P7_TRACE *tr)
    int  p7_trace_GrowIndex(P7_TRACE *tr)
    int  p7_trace_GrowTo(P7_TRACE *tr, int N)
    int  p7_trace_GrowIndexTo(P7_TRACE *tr, int ndom)
    void p7_trace_Destroy(P7_TRACE *tr)
    void p7_trace_DestroyArray(P7_TRACE **tr, int N)

    int  p7_trace_GetDomainCount   (const P7_TRACE *tr, int *ret_ndom)
    int  p7_trace_GetStateUseCounts(const P7_TRACE *tr, int *counts)
    int  p7_trace_GetDomainCoords  (const P7_TRACE *tr, int which, int *ret_i1, int *ret_i2, int *ret_k1, int *ret_k2)

    int   p7_trace_Validate(const P7_TRACE *tr, const ESL_ALPHABET *abc, const ESL_DSQ *dsq, char *errbuf)
    int   p7_trace_Dump(FILE *fp, const P7_TRACE *tr, const P7_PROFILE *gm, const ESL_DSQ *dsq)
    int   p7_trace_Compare(P7_TRACE *tr1, P7_TRACE *tr2, float pptol)
    int   p7_trace_Score(P7_TRACE *tr, ESL_DSQ *dsq, P7_PROFILE *gm, float *ret_sc)
    int   p7_trace_SetPP(P7_TRACE *tr, const P7_GMX *pp)
    float p7_trace_GetExpectedAccuracy(const P7_TRACE *tr)

    int  p7_trace_Append(P7_TRACE *tr, char st, int k, int i)
    int  p7_trace_AppendWithPP(P7_TRACE *tr, char st, int k, int i, float pp)
    int  p7_trace_Reverse(P7_TRACE *tr)
    int  p7_trace_Index(P7_TRACE *tr)

    int  p7_trace_FauxFromMSA(ESL_MSA *msa, int *matassign, int optflags, P7_TRACE **tr)
    int  p7_trace_Doctor(P7_TRACE *tr, int *opt_ndi, int *opt_nid)

    int  p7_trace_Count(P7_HMM *hmm, ESL_DSQ *dsq, float wt, P7_TRACE *tr)
