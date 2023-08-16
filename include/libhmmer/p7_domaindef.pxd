from libc.stdio cimport FILE

from libeasel.random cimport ESL_RANDOMNESS
from libeasel.sq cimport ESL_SQ
from libhmmer.p7_domain cimport P7_DOMAIN
from libhmmer.p7_spensemble cimport P7_SPENSEMBLE
from libhmmer.p7_trace cimport P7_TRACE
from libhmmer.p7_alidisplay cimport P7_ALIDISPLAY

if HMMER_IMPL == "VMX":
    from libhmmer.impl_vmx.p7_omx cimport P7_OMX
    from libhmmer.impl_vmx.p7_oprofile cimport P7_OPROFILE
elif HMMER_IMPL == "SSE":
    from libhmmer.impl_sse.p7_omx cimport P7_OMX
    from libhmmer.impl_sse.p7_oprofile cimport P7_OPROFILE
elif HMMER_IMPL == "NEON":
    from libhmmer.impl_neon.p7_omx cimport P7_OMX
    from libhmmer.impl_neon.p7_oprofile cimport P7_OPROFILE


cdef extern from "hmmer.h" nogil:

    cdef struct p7_domaindef_s:

        float* mocc
        float* btot
        float* etot
        int L
        int Lalloc

        float* n2sc

        ESL_RANDOMNESS* r
        int do_reseeding
        P7_SPENSEMBLE* sp
        P7_TRACE* tr
        P7_TRACE *gtr

        float rt1
        float rt2
        float rt3

        int nsamples
        float min_overlap
        int of_smaller
        int max_diagdiff
        float min_posterior
        float min_endpointp

        P7_DOMAIN* dcl
        int ndom
        int nalloc

        float nexpected
        int nregions
        int nclustered
        int noverlaps
        int nenvelopes
    ctypedef p7_domaindef_s P7_DOMAINDEF


    P7_DOMAINDEF *p7_domaindef_Create (ESL_RANDOMNESS *r)
    int           p7_domaindef_Fetch  (P7_DOMAINDEF *ddef, int which, int *opt_i, int *opt_j, float *opt_sc, P7_ALIDISPLAY **opt_ad)
    int           p7_domaindef_GrowTo (P7_DOMAINDEF *ddef, int L)
    int           p7_domaindef_Reuse  (P7_DOMAINDEF *ddef)
    int           p7_domaindef_DumpPosteriors(FILE *ofp, P7_DOMAINDEF *ddef)
    void          p7_domaindef_Destroy(P7_DOMAINDEF *ddef)

    # int p7_domaindef_ByViterbi            (P7_PROFILE *gm, const ESL_SQ *sq, const ESL_SQ *ntsq, P7_GMX *gx1, P7_GMX *gx2, P7_DOMAINDEF *ddef)
    # int p7_domaindef_ByPosteriorHeuristics(const ESL_SQ *sq, const ESL_SQ *ntsq, P7_OPROFILE *om, P7_OMX *oxf, P7_OMX *oxb, P7_OMX *fwd, P7_OMX *bck,
    # 				                                  P7_DOMAINDEF *ddef, P7_BG *bg, int long_target,
    # 				                                  P7_BG *bg_tmp, float *scores_arr, float *fwd_emissions_arr);
