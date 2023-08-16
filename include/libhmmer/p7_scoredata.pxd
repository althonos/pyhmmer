from libhmmer.p7_profile cimport P7_PROFILE

if HMMER_IMPL == "VMX":
    from libhmmer.impl_vmx.p7_oprofile cimport P7_OPROFILE
elif HMMER_IMPL == "SSE":
    from libhmmer.impl_sse.p7_oprofile cimport P7_OPROFILE
elif HMMER_IMPL == "NEON":
    from libhmmer.impl_neon.p7_oprofile cimport P7_OPROFILE

cdef extern from "hmmer.h" nogil:

    cdef enum p7_scoredatatype_e:
        p7_sd_std = 0
        p7_sd_fm = 1

    cdef struct p7_scoredata_s:
        p7_scoredatatype_e type
        int M
        float* prefix_lengths
        float* suffix_lengths
        float* fwd_scores
        float** fwd_transitions
        float** opt_ext_fwd
        float** opt_ext_rev
    ctypedef p7_scoredata_s P7_SCOREDATA

    cdef void p7_hmm_ScoreDataDestroy(P7_SCOREDATA* data)
    cdef P7_SCOREDATA* p7_hmm_ScoreDataCreate(P7_OPROFILE* om, P7_PROFILE* gm)
    cdef P7_SCOREDATA* p7_hmm_ScoreDataClone(P7_SCOREDATA* src, int Kp)
    cdef int p7_hmm_ScoreDataComputeRest(P7_OPROFILE* om, P7_SCOREDATA* data)
